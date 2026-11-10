"""Adapt Markdown math and links for Doxygen without editing sources.

Used only by Doxygen's input filter. Keep line counts unchanged so diagnostics
still refer to the original Markdown. Fenced and inline code are left intact.
"""

from pathlib import Path
import hashlib
import os
import re
import sys
from urllib.parse import urlsplit


TOKENS = re.compile(
    r"(?P<code>`+).*?(?P=code)|\\[()\[\]]|\\.|(?P<display_math>\$\$)|"
    r"(?P<inline_math>\$(?![\s$])(?:\\.|(?!\]\()[^$\\`\r\n])*?(?<!\s)\$(?![$\d]))|"
    r"\]\((?P<target>[^)\n]+)\)"
)
MATH = {r"\(": r"\f$", r"\)": r"\f$", r"\[": r"\f[", r"\]": r"\f]"}
MEDIA = {".gif", ".png", ".jpg", ".jpeg", ".svg", ".mp4", ".webm"}


def adapt_markdown(text, source):
    """Translate supported math delimiters and flatten existing local media URLs."""
    lines = []
    fence = None
    for line in text.splitlines(keepends=True):
        marker = re.match(r"^ {0,3}(`{3,}|~{3,})", line)
        if marker:
            value = marker.group(1)
            if fence is None:
                fence = value
            elif value[0] == fence[0] and len(value) >= len(fence) and not line[marker.end():].strip():
                fence = None
            lines.append((line, True))
        else:
            lines.append((line, fence is not None))
    delimiters = [
        (index, match.start())
        for index, (line, protected) in enumerate(lines)
        if not protected
        for match in TOKENS.finditer(line)
        if match.group("display_math")
    ]
    # Leave an unmatched delimiter unchanged rather than opening an unterminated Doxygen formula.
    display_markers = {
        position: r"\f[" if index % 2 == 0 else r"\f]"
        for index, position in enumerate(delimiters[:len(delimiters) // 2 * 2])
    }
    display_active = False

    def replace(match, line_index):
        nonlocal display_active
        token = match.group()
        if match.group("display_math"):
            replacement = display_markers.get((line_index, match.start()), token)
            if replacement != token:
                display_active = not display_active
            return replacement
        if match.group("inline_math"):
            return token if display_active else r"\f$" + token[1:-1] + r"\f$"
        if token in MATH:
            return MATH[token]
        target = match.group("target")
        if target:
            url = urlsplit(target)
            path = source.parent / url.path
            if not os.environ.get("FUNDEM_API_SOURCE_ROOT") and not url.scheme and not url.netloc and path.suffix.lower() in MEDIA and path.is_file():
                return "](" + url._replace(path=path.name).geturl() + ")"
        return token

    result = [
        line if protected else TOKENS.sub(lambda match: replace(match, index), line)
        for index, (line, protected) in enumerate(lines)
    ]
    source_root = os.environ.get("FUNDEM_API_SOURCE_ROOT")
    if source_root:
        relative = source.resolve().relative_to(Path(source_root).resolve()).as_posix()
        anchor = "fundem_api_" + hashlib.sha256(relative.encode("utf-8")).hexdigest()[:20]
        # The first non-fenced blank line preserves original diagnostic line numbers.
        fenced = None
        seen_content = False
        for index, line in enumerate(result):
            marker = re.match(r"^ {0,3}(`{3,}|~{3,})", line)
            if marker:
                value = marker.group(1)
                if fenced is None:
                    fenced = value
                elif value[0] == fenced[0] and len(value) >= len(fenced) and not line[marker.end():].strip():
                    fenced = None
            elif not fenced and seen_content and not line.strip():
                ending = "\r\n" if line.endswith("\r\n") else "\n" if line.endswith("\n") else ""
                result[index] = "\\anchor " + anchor + ending
                break
            seen_content = seen_content or bool(line.strip())
        else:
            raise ValueError(f"No blank line is available for an API page marker: {source}")
    return "".join(result)


if __name__ == "__main__":
    path = Path(sys.argv[1])
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stdout.write(adapt_markdown(path.read_text(encoding="utf-8"), path))
