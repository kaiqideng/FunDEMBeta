"""Build API HTML and repository attachments outside the source tree."""

import argparse
import hashlib
import html
from html.parser import HTMLParser
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from urllib.parse import quote, unquote, urlsplit, urlunsplit
import xml.etree.ElementTree as ET


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def page_marker(relative):
    return "fundem_api_" + hashlib.sha256(relative.encode("utf-8")).hexdigest()[:20]


def source_files(root):
    return {
        p.relative_to(root).as_posix(): digest(p)
        for p in root.rglob("*")
        if p.is_file()
        and not any(part == ".git" or part == "build" or part.startswith("build-") for part in p.relative_to(root).parts[:-1])
        and not p.is_relative_to(root / "docs" / "api")
    }


class Links(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.anchors = set()
        self.targets = []
        self.headings = {}
        self.heading = None
        self.slug_counts = {}

    def handle_starttag(self, tag, attributes):
        attrs = dict(attributes)
        if re.fullmatch(r"h[1-6]", tag):
            self.heading = {"text": [], "anchor": attrs.get("id")}
        if self.heading is not None and attrs.get("id"):
            self.heading["anchor"] = attrs["id"]
        self.anchors.update(attrs[key] for key in ("id", "name") if attrs.get(key))
        self.targets.extend(attrs[key] for key in ("href", "src") if attrs.get(key))

    def handle_data(self, text):
        if self.heading is not None:
            self.heading["text"].append(text)

    def handle_endtag(self, tag):
        if re.fullmatch(r"h[1-6]", tag) and self.heading is not None:
            if self.heading["anchor"]:
                title = " ".join("".join(self.heading["text"]).split()).lower()
                title = re.sub(r"<[^>]+>", "", title)
                slug = re.sub(r"[^\w -]", "", title).replace(" ", "-")
                count = self.slug_counts.get(slug, 0)
                self.slug_counts[slug] = count + 1
                self.headings[slug + (f"-{count}" if count else "")] = self.heading["anchor"]
            self.heading = None


def repair_fragments(folder, attributes):
    pages = {}
    for path in folder.rglob("*.html"):
        document = Links()
        document.feed(path.read_text(encoding="utf-8"))
        pages[path.resolve()] = document
    repaired = 0
    for page in pages:
        def replace(match):
            nonlocal repaired
            url = urlsplit(html.unescape(match.group(3)))
            if url.scheme or url.netloc or not url.fragment:
                return match.group()
            target = (page.parent / unquote(url.path)).resolve() if url.path else page
            document = pages.get(target)
            fragment = unquote(url.fragment)
            if not document or fragment in document.anchors or fragment not in document.headings:
                return match.group()
            repaired += 1
            target_url = urlunsplit(("", "", url.path, url.query, quote(document.headings[fragment], safe="-_.:")))
            return match.group(1) + "=" + match.group(2) + html.escape(target_url, quote=True) + match.group(2)
        text = page.read_text(encoding="utf-8")
        updated = attributes.sub(replace, text)
        if updated != text:
            page.write_text(updated, encoding="utf-8")
    return repaired


def check_links(folder):
    pages = {}
    for path in folder.rglob("*.html"):
        parser = Links()
        parser.feed(path.read_text(encoding="utf-8"))
        pages[path.resolve()] = parser
    errors = []
    count = 0
    for path, parser in pages.items():
        for target in parser.targets:
            url = urlsplit(target)
            if url.scheme or url.netloc or not target:
                continue
            count += 1
            resolved = (path.parent / unquote(url.path)).resolve() if url.path else path
            if not resolved.is_relative_to(folder.resolve()) or not resolved.is_file():
                errors.append({"page": str(path.relative_to(folder)), "target": target, "error": "missing file"})
            elif url.fragment and resolved.suffix.lower() == ".html":
                document = pages.get(resolved)
                if document and unquote(url.fragment) not in document.anchors:
                    errors.append({"page": str(path.relative_to(folder)), "target": target, "error": "missing anchor"})
    return {"html_pages": len(pages), "checked_local_links": count, "errors": errors}


def main():
    root = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=root.parent / "build-api", help="new directory outside the source tree")
    args = parser.parse_args()
    output = args.output.expanduser().resolve()
    if output.is_relative_to(root) or root.is_relative_to(output):
        parser.error("output must be outside the source tree and cannot contain it")
    if output.exists():
        parser.error("output already exists; choose a fresh directory")
    if not shutil.which("doxygen"):
        parser.error("doxygen is not available on PATH")
    output.mkdir(parents=True)
    before = source_files(root)
    config = output / "Doxyfile"
    config.write_text(
        (root / "Doxyfile").read_text(encoding="utf-8")
        + f'\nOUTPUT_DIRECTORY = "{output}"\nGENERATE_TAGFILE = "{output / "tags.xml"}"\nSOURCE_BROWSER = NO\n',
        encoding="utf-8",
    )
    environment = dict(os.environ, FUNDEM_API_SOURCE_ROOT=str(root))
    log = output / "doxygen.log"
    with log.open("w", encoding="utf-8") as stream:
        result = subprocess.run(["doxygen", str(config)], cwd=root, env=environment, stdout=stream, stderr=subprocess.STDOUT)
    if result.returncode:
        print(f"Doxygen failed; see {log}", file=sys.stderr)
        return result.returncode

    folder = output / "html"
    tags = ET.parse(output / "tags.xml").getroot()
    markers = {page_marker(name): root / name for name in before if name.endswith(".md")}
    pages = {}
    for anchor in tags.iter("docanchor"):
        if anchor.text in markers:
            pages[markers[anchor.text]] = anchor.attrib["file"]
    contexts = {folder / page: source for source, page in pages.items()}
    for compound in tags.findall("compound"):
        if compound.attrib.get("kind") != "file":
            continue
        name = compound.findtext("name", "")
        source = root / name
        if not source.is_file():
            source = root / compound.findtext("path", "") / Path(name).name
        if source.is_file() and source.resolve().is_relative_to(root):
            contexts[folder / compound.findtext("filename")] = source

    attachments = {}
    rewritten = 0
    attributes = re.compile(r'''\b(href|src)=(['"])(.*?)\2''')
    for page, source in contexts.items():
        if not page.is_file():
            continue

        def replace(match):
            nonlocal rewritten
            url = urlsplit(html.unescape(match.group(3)))
            if url.scheme or url.netloc or not url.path:
                return match.group()
            generated = (page.parent / unquote(url.path)).resolve()
            if generated.is_file() and generated.is_relative_to(folder):
                return match.group()
            target = (source.parent / unquote(url.path)).resolve()
            if not target.is_relative_to(root) or not target.is_file():
                return match.group()
            if target in pages:
                destination = folder / pages[target]
            else:
                relative = target.relative_to(root).as_posix()
                destination = folder / "assets" / "source" / relative
                if relative not in attachments:
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(target, destination)
                    attachments[relative] = digest(destination)
            path = quote(os.path.relpath(destination, page.parent).replace(os.sep, "/"), safe="/")
            target_url = urlunsplit(("", "", path, url.query, url.fragment))
            rewritten += 1
            return match.group(1) + "=" + match.group(2) + html.escape(target_url, quote=True) + match.group(2)

        text = page.read_text(encoding="utf-8")
        updated = attributes.sub(replace, text)
        if updated != text:
            page.write_text(updated, encoding="utf-8")

    repaired_fragments = repair_fragments(folder, attributes)
    links = check_links(folder)
    after = source_files(root)
    changed = sorted(name for name in before.keys() | after.keys() if before.get(name) != after.get(name))
    audit = {
        "source": str(root), "output": str(output),
        "source_files_before": before, "source_changes_during_build": changed,
        "rendered_markdown_pages": len(pages), "rewritten_links": rewritten,
        "repaired_heading_fragments": repaired_fragments,
        "copied_attachments_sha256": attachments, "local_links": links,
        "all_checks_passed": not changed and not links["errors"],
    }
    (output / "api-build-audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    if not audit["all_checks_passed"]:
        print(f"API checks failed; see {output / 'api-build-audit.json'}", file=sys.stderr)
        return 1
    print(f"API HTML: {folder / 'index.html'}")
    print(f"Checked {links['checked_local_links']} local links; copied {len(attachments)} repository attachments.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
