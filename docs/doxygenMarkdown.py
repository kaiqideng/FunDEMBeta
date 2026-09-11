"""Adapt Markdown math and local media links for Doxygen without editing sources.

Used only by Doxygen's input filter. Keep line counts unchanged so diagnostics
still refer to the original Markdown. Fenced and inline code are left intact.
"""

from pathlib import Path
import re
import sys
from urllib.parse import urlsplit


TOKENS = re.compile(r"(`+).*?\1|\\[()\[\]]|\]\(([^)\n]+)\)")
MATH = {r"\(": r"\f$", r"\)": r"\f$", r"\[": r"\f[", r"\]": r"\f]"}
MEDIA = {".gif", ".png", ".jpg", ".jpeg", ".svg", ".mp4", ".webm"}


def adapt_markdown(text, source):
    """Translate supported math delimiters and flatten existing local media URLs."""
    def replace(match):
        token = match.group()
        if token in MATH:
            return MATH[token]
        target = match.group(2)
        if target:
            url = urlsplit(target)
            path = source.parent / url.path
            if not url.scheme and not url.netloc and path.suffix.lower() in MEDIA and path.is_file():
                return "](" + url._replace(path=path.name).geturl() + ")"
        return token

    result = []
    fence = None
    for line in text.splitlines(keepends=True):
        marker = re.match(r"^ {0,3}(`{3,}|~{3,})", line)
        if marker:
            value = marker.group(1)
            if fence is None:
                fence = value
            elif value[0] == fence[0] and len(value) >= len(fence) and not line[marker.end():].strip():
                fence = None
            result.append(line)
        else:
            result.append(line if fence else TOKENS.sub(replace, line))
    return "".join(result)


if __name__ == "__main__":
    path = Path(sys.argv[1])
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stdout.write(adapt_markdown(path.read_text(encoding="utf-8"), path))
