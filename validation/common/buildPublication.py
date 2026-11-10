#!/usr/bin/env python3
"""Typeset the retained verification reports as a vector-figure technical volume.

Requires matplotlib, reportlab and pypdf. Numerical records are never modified.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import html
import json
from pathlib import Path
import re
import textwrap
from urllib.parse import urlsplit

import matplotlib
matplotlib.use('Agg')
from matplotlib import mathtext
from matplotlib.font_manager import FontProperties
from pypdf import PdfReader, PdfWriter, Transformation
from pypdf.generic import ContentStream, DictionaryObject, NameObject
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen.canvas import Canvas
from reportlab.platypus import (BaseDocTemplate, Flowable, Frame, KeepTogether,
    PageBreak, PageTemplate, Paragraph, Spacer, Table, TableStyle)
from reportlab.platypus.tableofcontents import TableOfContents

CORE = Path(__file__).resolve().parents[2]
VALIDATION = CORE / 'validation'
WIDTH = 180 * mm
INK = colors.HexColor('#202020')
MUTED = colors.HexColor('#555555')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def clean(text):
    return text.translate(str.maketrans({'\u2010': '-', '\u2011': '-', '\u2012': '-',
        '\u2013': '-', '\u2014': '-', '\u2212': '-', '\u00a0': ' '}))


class VectorPDF(Flowable):
    def __init__(self, path, width, placements):
        super().__init__()
        self.path = Path(path)
        page = PdfReader(self.path).pages[0]
        self.original = (float(page.mediabox.width), float(page.mediabox.height))
        self.width = width
        self.height = width * self.original[1] / self.original[0]
        self.placements = placements

    def draw(self):
        x, y = self.canv.absolutePosition(0, 0)
        self.placements.append({'page': self.canv.getPageNumber(), 'x': x, 'y': y,
            'width': self.width, 'height': self.height, 'pdf': str(self.path)})


class NumberedCanvas(Canvas):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault('initialFontName', 'STIX')
        super().__init__(*args, **kwargs)
        self._case_label = 'Scope and contents'

    def showPage(self):
        self.saveState()
        self.setStrokeColor(colors.HexColor('#999999'))
        self.setLineWidth(.35)
        self.line(15 * mm, 283 * mm, 195 * mm, 283 * mm)
        self.setFont('STIX', 8)
        self.setFillColor(MUTED)
        self.drawString(15 * mm, 286 * mm, 'FUNDEM  |  VERIFICATION REPORT SERIES')
        self.drawRightString(195 * mm, 286 * mm, self._case_label)
        self.drawString(15 * mm, 11 * mm, 'Historical evidence; publication-format edition, 31 October 2026')
        self.drawRightString(195 * mm, 11 * mm, str(self._pageNumber))
        self.restoreState()
        Canvas.showPage(self)


class ReportDoc(BaseDocTemplate):
    def __init__(self, path, placements, **kwargs):
        super().__init__(str(path), pagesize=A4, leftMargin=15 * mm,
            rightMargin=15 * mm, topMargin=20 * mm, bottomMargin=19 * mm, **kwargs)
        self.placements = placements
        self.addPageTemplates(PageTemplate(id='paper', frames=[Frame(15 * mm, 19 * mm,
            WIDTH, 258 * mm, leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0)]))

    def beforeDocument(self):
        self.placements.clear()

    def afterFlowable(self, flowable):
        if hasattr(flowable, 'case_key'):
            self.canv._case_label = flowable.case_id
            self.canv.bookmarkPage(flowable.case_key)
            title = flowable.case_id + ': ' + flowable.navigation_title
            self.canv.addOutlineEntry(title, flowable.case_key, 0)
            self.notify('TOCEntry', (0, title, self.page, flowable.case_key))


class Typesetter:
    def __init__(self, assets, placements):
        self.assets = assets
        self.placements = placements
        fonts = Path(matplotlib.get_data_path()) / 'fonts/ttf'
        for name, filename in [('STIX', 'STIXGeneral.ttf'), ('STIX-Bold', 'STIXGeneralBol.ttf'),
                               ('STIX-Italic', 'STIXGeneralItalic.ttf'), ('STIX-BoldItalic', 'STIXGeneralBolIta.ttf'),
                               ('Code', 'DejaVuSansMono.ttf')]:
            pdfmetrics.registerFont(TTFont(name, str(fonts / filename)))
        pdfmetrics.registerFontFamily('STIX', normal='STIX', bold='STIX-Bold',
            italic='STIX-Italic', boldItalic='STIX-BoldItalic')
        self.styles = {
            'body': ParagraphStyle('body', fontName='STIX', fontSize=10.5, leading=14,
                textColor=INK, spaceAfter=6, splitLongWords=True, autoLeading='max'),
            'title': ParagraphStyle('title', fontName='STIX-Bold', fontSize=18, leading=22,
                spaceBefore=4, spaceAfter=12, keepWithNext=True),
            'h2': ParagraphStyle('h2', fontName='STIX-Bold', fontSize=12, leading=15,
                spaceBefore=12, spaceAfter=5, keepWithNext=True),
            'h3': ParagraphStyle('h3', fontName='STIX-Bold', fontSize=10.6, leading=14,
                spaceBefore=8, spaceAfter=4, keepWithNext=True),
            'caption': ParagraphStyle('caption', fontName='STIX', fontSize=9, leading=11.5,
                spaceBefore=4, spaceAfter=9, textColor=INK, autoLeading='max'),
            'code': ParagraphStyle('code', fontName='Code', fontSize=7.3, leading=10,
                spaceBefore=5, spaceAfter=8, leftIndent=5, rightIndent=5, backColor=colors.HexColor('#f5f5f5')),
            'small': ParagraphStyle('small', fontName='STIX', fontSize=9, leading=12,
                spaceAfter=6, textColor=MUTED),
            'appendix': ParagraphStyle('appendix', fontName='STIX', fontSize=9, leading=12,
                spaceAfter=5, textColor=INK, autoLeading='max'),
        }
        self.math_cache = {}
        self.table_audit = []
        matplotlib.rcParams.update({'mathtext.fontset': 'stix', 'font.family': 'STIXGeneral',
            'pdf.fonttype': 42, 'svg.fonttype': 'none'})

    def equation(self, expression, size=10.5, vector=False):
        expression = expression.strip().replace(r'\tfrac', r'\frac').replace(r'\dfrac', r'\frac')
        expression = re.sub(r'\\tag\{[^}]*\}', '', expression).strip()
        key = hashlib.sha256((expression + str(size) + str(vector)).encode()).hexdigest()[:20]
        if key not in self.math_cache:
            filename = self.assets / ('eq-' + key + ('.pdf' if vector else '.png'))
            prop = FontProperties(size=size)
            parsed = mathtext.MathTextParser('path').parse('$' + expression + '$', dpi=72, prop=prop)
            mathtext.math_to_image('$' + expression + '$', filename, prop=prop,
                dpi=600, format='pdf' if vector else 'png', color='#202020')
            self.math_cache[key] = (filename, float(parsed.width), float(parsed.height), float(parsed.depth))
        return self.math_cache[key]

    def inline(self, text, source, size=10.5):
        protected = []
        def keep(value):
            protected.append(value)
            return f'ZZTOKEN{len(protected)-1}ZZ'
        def math_replace(match):
            expression = match.group(1) if match.group(1) is not None else match.group(2)
            path, width, height, depth = self.equation(expression, size=size)
            return keep(f'<img src="{path}" width="{width:.3f}" height="{height:.3f}" valign="{-depth:.3f}"/>')
        text = re.sub(r'<br\s*/?>', lambda m: keep('<br/>'), text)
        superscript = str.maketrans('⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺', '0123456789-+')
        text = re.sub(r'[⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺]+', lambda m: keep('<super>' + m[0].translate(superscript) + '</super>'), text)
        text = re.sub(r'\\\((.*?)\\\)|(?<!\$)\$(?!\$)(.+?)(?<!\$)\$(?!\$)', math_replace, clean(text))
        text = re.sub(r'`([^`]+)`', lambda m: keep('<font face="Code" size="8">' + html.escape(m[1]) + '</font>'), text)
        def link_replace(match):
            label, url = match[1], match[2].strip('<>')
            if not urlsplit(url).scheme:
                target = (source.parent / url.split('#')[0]).resolve()
                fragment = '#' + url.split('#', 1)[1] if '#' in url else ''
                url = ('validation/' + target.relative_to(VALIDATION).as_posix()
                    if target.is_relative_to(VALIDATION) else target.as_uri()) + fragment
            return keep('<link href="' + html.escape(url, quote=True) + '" color="#204a70">' + html.escape(label) + '</link>')
        text = re.sub(r'\[([^\]]+)\]\(([^)]+)\)', link_replace, text)
        text = html.escape(text)
        text = re.sub(r'\*\*(.+?)\*\*', r'<b>\1</b>', text)
        text = re.sub(r'(?<!\*)\*([^*]+)\*(?!\*)', r'<i>\1</i>', text)
        for i in reversed(range(len(protected))):
            text = text.replace(f'ZZTOKEN{i}ZZ', protected[i])
        return text

    def paragraph(self, text, source, style='body'):
        return Paragraph(self.inline(text, source, self.styles[style].fontSize), self.styles[style])

    def table(self, lines, source):
        def cells(line):
            return [value.strip().replace(r'\|', '|') for value in re.split(r'(?<!\\)\|', line.strip().strip('|'))]
        data = [cells(line) for line in lines]
        data = [row for row in data if not all(re.fullmatch(r':?-{2,}:?', value.replace(' ', '')) for value in row)]
        columns = len(data[0])
        if any(len(row) != columns for row in data):
            raise ValueError(f'Invalid Markdown table: {source}')
        if columns == 2:
            widths = [WIDTH * .30, WIDTH * .70]
        else:
            weights = [max(8, min(23, max(len(re.sub(r'[`*]', '', row[i])) for row in data[:8]))) for i in range(columns)]
            widths = [WIDTH * weight / sum(weights) for weight in weights]
        font_size = 8.2 if columns <= 7 else 7.4
        # A TeX header's source length is unrelated to its printed width. Give
        # each numeric token enough physical space before distributing the
        # remaining width; a decimal or exponent must never wrap mid-number.
        numeric = re.compile(r'(?<![A-Za-z0-9_])[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?%?')
        minimum = [max([18.0] + [pdfmetrics.stringWidth(token, 'STIX', font_size) + 9
            for row in data[1:] for token in numeric.findall(clean(row[i]))]) for i in range(columns)]
        if sum(minimum) > WIDTH:
            raise ValueError(f'Table numeric tokens cannot fit the publication width: {source}')
        pending = set(range(columns))
        fixed = {}
        while pending:
            remaining = WIDTH - sum(fixed.values())
            total_weight = sum(widths[i] for i in pending)
            proposed = {i: remaining * widths[i] / total_weight for i in pending}
            cramped = {i for i in pending if proposed[i] < minimum[i]}
            if not cramped:
                fixed.update(proposed)
                break
            for i in cramped:
                fixed[i] = minimum[i]
            pending.difference_update(cramped)
        widths = [fixed[i] for i in range(columns)]
        base = ParagraphStyle('cell', parent=self.styles['body'], fontSize=font_size,
            leading=font_size + 2.3, spaceAfter=0, spaceBefore=0)
        head = ParagraphStyle('headcell', parent=base, fontName='STIX-Bold')
        rendered = [[Paragraph(self.inline(cell, source, font_size), head if i == 0 else base) for cell in row] for i, row in enumerate(data)]
        numeric_cells = 0
        for i, row in enumerate(data[1:], 1):
            for j, cell in enumerate(row):
                value = clean(cell).strip('`* ')
                if numeric.fullmatch(value):
                    _, height = rendered[i][j].wrap(widths[j] - 8, 10000)
                    if height > base.leading + .01:
                        raise ValueError(f'Numeric table cell would wrap: {source}: row {i + 1}, column {j + 1}: {value}')
                    numeric_cells += 1
        self.table_audit.append({'source': str(source.relative_to(CORE)), 'columns': columns,
            'rows': len(data), 'column_widths_pt': widths, 'unwrapped_numeric_cells_checked': numeric_cells})
        table = Table(rendered, colWidths=widths, repeatRows=1, hAlign='LEFT')
        table.setStyle(TableStyle([('VALIGN', (0,0),(-1,-1),'TOP'), ('LEFTPADDING',(0,0),(-1,-1),4),
            ('RIGHTPADDING',(0,0),(-1,-1),4), ('TOPPADDING',(0,0),(-1,-1),4),
            ('BOTTOMPADDING',(0,0),(-1,-1),4), ('LINEABOVE',(0,0),(-1,0),.7,INK),
            ('LINEBELOW',(0,0),(-1,0),.4,INK), ('LINEBELOW',(0,-1),(-1,-1),.7,INK)]))
        return table

    def markdown(self, source):
        lines = source.read_text(encoding='utf-8').splitlines()
        story = []
        index = 0
        in_appendix = False
        while index < len(lines):
            line = lines[index].strip()
            if not line or line.startswith('# '):
                index += 1; continue
            if line.startswith('```'):
                block = []
                index += 1
                while index < len(lines) and not lines[index].startswith('```'):
                    block.extend(textwrap.wrap(clean(lines[index]), width=106, subsequent_indent='    ',
                        replace_whitespace=False, drop_whitespace=False, break_long_words=False, break_on_hyphens=False) or [''])
                    index += 1
                story.append(Paragraph('<br/>'.join(html.escape(x).replace(' ', '&#160;') for x in block), self.styles['code']))
                index += 1; continue
            if line.startswith('$$') or line == r'\[':
                block = [line[2:]] if line.startswith('$$') else []
                index += 1
                if block and block[0].endswith('$$'):
                    block[0] = block[0][:-2]
                else:
                    while index < len(lines) and lines[index].strip() not in ('$$', r'\]'):
                        block.append(lines[index]); index += 1
                    index += 1
                expression = ' '.join(block).strip()
                number = re.search(r'\\tag\{([^}]*)\}', expression)
                path, width, height, _ = self.equation(expression, size=11, vector=True)
                image = VectorPDF(path, min(width + .3, WIDTH - 38), self.placements)
                image.hAlign = 'CENTER'
                row = Table([[image, Paragraph('(' + number[1] + ')' if number else '', self.styles['body'])]],
                    colWidths=[WIDTH - 30, 30])
                row.setStyle(TableStyle([('VALIGN',(0,0),(-1,-1),'MIDDLE'),('ALIGN',(0,0),(0,0),'CENTER'),
                    ('LEFTPADDING',(0,0),(-1,-1),0),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5)]))
                story.append(row); continue
            if line.startswith('|'):
                block = []
                while index < len(lines) and lines[index].strip().startswith('|'):
                    block.append(lines[index]); index += 1
                story.extend([self.table(block, source), Spacer(1, 7)])
                continue
            image = re.fullmatch(r'!\[[^\]]*\]\(([^)]+)\)', line)
            if image:
                path = (source.parent / image[1]).with_suffix('.pdf')
                if not path.is_file():
                    raise FileNotFoundError(path)
                item = VectorPDF(path, WIDTH, self.placements)
                item.keepWithNext = True
                story.append(item); index += 1; continue
            heading = re.match(r'^(#{2,})\s+(.*)', line)
            if heading:
                if len(heading[1]) == 2 and heading[2].startswith('Appendix '):
                    in_appendix = True
                story.append(self.paragraph(heading[2], source, 'h2' if len(heading[1]) == 2 else 'h3'))
                index += 1; continue
            bullet = re.match(r'^[-*]\s+(.*)', line)
            block = [bullet[1] if bullet else line]
            index += 1
            while index < len(lines) and lines[index].strip() and not re.match(r'^(#{1,6}\s|\||```|\$\$|!\[|[-*]\s)', lines[index].strip()):
                block.append(lines[index].strip()); index += 1
            text = ' '.join(block)
            style = 'caption' if re.match(r'^(?:\*\*)?(?:Figure|Table)\s+\d+', text) else ('appendix' if in_appendix else 'body')
            paragraph = self.paragraph(('- ' if bullet else '') + text, source, style)
            if style == 'caption' and text.lstrip('*').startswith('Table'):
                paragraph.keepWithNext = True
            story.append(paragraph)
        grouped = []
        index = 0
        while index < len(story):
            item = story[index]
            if isinstance(item, Paragraph) and item.style.name == 'h2' and item.getPlainText() == 'References':
                block = [item]
                index += 1
                while index < len(story) and not (isinstance(story[index], Paragraph) and story[index].style.name == 'h2'):
                    block.append(story[index]); index += 1
                grouped.append(KeepTogether(block))
            else:
                grouped.append(item); index += 1
        return grouped


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path, help='new PDF path outside the core source tree')
    args = parser.parse_args()
    output = args.output.expanduser().resolve()
    if output.is_relative_to(CORE) or output.exists():
        parser.error('choose a new PDF path outside the source tree')
    output.parent.mkdir(parents=True, exist_ok=True)
    assets = output.parent / (output.stem + '-assets')
    assets.mkdir(exist_ok=False)
    placements = []
    typesetter = Typesetter(assets, placements)
    cases = []
    for family, prefix in [('LAMMPSComparison', 'L'), ('PFCComparison', 'P')]:
        rows = list(csv.DictReader((VALIDATION / family / 'HISTORICAL-INDEX.csv').open(encoding='utf-8-sig')))
        for i, record in enumerate(rows, 1):
            cases.append((f'{prefix}{i:02d}', record['case'], VALIDATION / family / 'cases' / record['case'] / 'REPORT.md'))
    cases.append(('O01', 'ondem-bumpy-wall-shear', VALIDATION / 'ONDEMComparison/shearFlow/README.md'))
    inputs = {str(path.relative_to(CORE)): digest(path) for _, _, path in cases}
    front = VALIDATION / 'PUBLICATION.md'
    story = [Spacer(1, 22 * mm), Paragraph('Verification of FunDEM<br/>Contact Models',
        ParagraphStyle('cover', parent=typesetter.styles['title'], fontSize=28, leading=33)),
        typesetter.paragraph('Reference benchmarks, analytical comparisons and documented limitations', front, 'h2'),
        Spacer(1, 12 * mm), typesetter.paragraph('PUBLICATION-FORMAT TECHNICAL REPORT SERIES', front, 'small'),
        typesetter.paragraph('Edition: 31 October 2026<br/>Historical LAMMPS/PFC campaign: 28-29 October 2026', front),
        Spacer(1, 12 * mm), typesetter.paragraph('35 LAMMPS reference cases: 33 historical PASS and 2 retained FAIL.<br/>3 PFC analytical comparisons: 3 historical PASS; no local PFC execution.<br/>1 ON-DEM adaptation: unresolved transient discrepancies; full benchmark not run.', front),
        Spacer(1, 10 * mm), typesetter.paragraph('The volume reformats retained evidence. It does not certify the current solver revision, establish experimental validation, or claim peer review. Exact records and qualified failures remain available with each case.', front),
        PageBreak(), typesetter.paragraph('Scope and notation', front, 'title'),
        typesetter.paragraph('The comparisons distinguish analytical checks, published saved-state comparisons and external benchmark adaptations. Agreement of CPU and GPU outputs measures implementation consistency; it is not an independent proof of physical-model equivalence. Two time steps alone do not establish an asymptotic convergence order.', front),
        typesetter.paragraph('Reference identity, source revision, contact-law mapping, observable definition, tolerance and physical check time are specified in each report. Tables, equations and figures are numbered locally within a case. Sparse checkpoint symbols do not represent complete trajectories. No new solver runs were performed for this presentation edition.', front),
        typesetter.paragraph('Absolute error is |q - q_ref|. Relative error divides this magnitude by |q_ref| only for a nonzero reference. A normalized acceptance ratio divides the stated error by its case-specific tolerance; the acceptance inequality is reproduced from the case evidence. PFC restitution uses an absolute coefficient error below 0.01. The same-code YAML epsilon is not promoted to a cross-software criterion.', front),
        typesetter.paragraph('Figures use a 180 mm final width, embedded vector text and graphics, and separate 600 dpi raster exports. Sources and reproducible generators accompany the reports. Relative evidence hyperlinks resolve in the accompanying validation/ directory; external references use their cited source URLs.', front),
        PageBreak(), typesetter.paragraph('Contents', front, 'title')]
    toc = TableOfContents()
    toc.levelStyles = [ParagraphStyle('toc', parent=typesetter.styles['body'], fontSize=8.5, leading=10.5, spaceAfter=0)]
    story.append(toc)
    for key, case_id, source in cases:
        story.append(PageBreak())
        title = source.read_text(encoding='utf-8').splitlines()[0].lstrip('# ')
        story.append(typesetter.paragraph(f'CASE {key}  |  {case_id}', source, 'small'))
        heading = typesetter.paragraph(title, source, 'title')
        heading.case_key, heading.case_id = key, case_id
        heading.navigation_title = title.replace(r'\(', '').replace(r'\)', '').replace(r'\beta', 'beta')
        story.append(heading)
        story.extend(typesetter.markdown(source))
    base = assets / 'typeset-base.pdf'
    doc = ReportDoc(base, placements, title='Verification of FunDEM Contact Models',
        author='FunDEM', subject='Historical reference comparisons and limitations')
    doc.multiBuild(story, canvasmaker=NumberedCanvas)
    reader = PdfReader(base)
    writer = PdfWriter(clone_from=reader)
    # ReportLab shares resource dictionaries. Isolate them before page overlays
    # so one page does not accumulate another page's vector font aliases.
    for page in writer.pages:
        resources = DictionaryObject(page['/Resources'].items())
        for key in ('/Font', '/XObject', '/ExtGState', '/ColorSpace', '/Pattern', '/Shading'):
            if key in resources and isinstance(resources[key], DictionaryObject):
                resources[NameObject(key)] = DictionaryObject(resources[key].items())
        page[NameObject('/Resources')] = resources
    for item in placements:
        page = PdfReader(item['pdf']).pages[0]
        scale = item['width'] / float(page.mediabox.width)
        writer.pages[item['page'] - 1].merge_transformed_page(page,
            Transformation().scale(scale).translate(item['x'], item['y']))
    # ReportLab declares an unused base Helvetica font in its initial text state.
    # Map that initial state to the embedded body font and drop the base resource.
    for page in writer.pages:
        fonts = page['/Resources'].get('/Font', {}).get_object()
        fallback = next((key for key, ref in fonts.items()
            if ref.get_object().get('/Subtype') == '/TrueType'
            and str(ref.get_object().get('/BaseFont')).endswith('STIXGeneral-Regular')), None)
        base_fonts = {key for key, ref in fonts.items() if ref.get_object().get('/BaseFont') == '/Helvetica'}
        if base_fonts:
            if fallback is None:
                raise RuntimeError('Embedded body font is unavailable for the initial PDF text state')
            content = ContentStream(page.get_contents(), writer)
            for operands, operator in content.operations:
                if operator == b'Tf' and operands[0] in base_fonts:
                    operands[0] = NameObject(fallback)
            page.replace_contents(content)
            for key in base_fonts:
                del fonts[key]
        page.compress_content_streams()
    writer.add_metadata({'/Title': 'Verification of FunDEM Contact Models', '/Author': 'FunDEM',
        '/Subject': 'Publication-format historical verification reports; qualified outcomes retained'})
    with output.open('wb') as stream:
        writer.write(stream)
    changed = [name for name, value in inputs.items() if digest(CORE / name) != value]
    if changed:
        raise RuntimeError('Report sources changed during layout: ' + ', '.join(changed))
    audit = {'pdf': str(output), 'sha256': digest(output), 'pages': len(reader.pages),
        'reports': len(cases), 'vector_placements': placements, 'source_reports_sha256': inputs,
        'source_changes_during_build': changed, 'fonts': 'STIX General, DejaVu Sans Mono; embedded TrueType',
        'page_size': 'A4', 'text_width_mm': 180, 'body_font_pt': 10.5, 'appendix_font_pt': 9,
        'table_layout_checks': typesetter.table_audit,
        'scope': 'Typesetting and provenance only; visual review is required before publication delivery.'}
    (output.parent / 'publication-build-audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    print(json.dumps({k: audit[k] for k in ('pdf', 'sha256', 'pages', 'reports')}, indent=2))


if __name__ == '__main__':
    main()
