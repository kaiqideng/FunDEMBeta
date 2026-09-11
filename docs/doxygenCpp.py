"""Expand storage aliases only in Doxygen's C++ parser input.

Doxygen 1.9.8 drops template arguments while resolving aliases of the
hostAoSDeviceSoA partial specialization. It consequently merges overloads with
different particle types. Spelling out the same underlying type preserves those
overloads. Discover aliases from module headers so new storage types are handled
without maintaining a separate list of affected functions or signatures.

Alias declarations, comments, strings, and line counts remain intact. Keep
FILTER_SOURCE_FILES disabled so source listings show the actual C++ files.
"""

from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parent.parent
NON_CODE = re.compile(
    r'(?:(?:u8|u|U|L)?R)"(?P<delimiter>[^ ()\\\t\r\n]{0,16})'
    r'\(.*?\)(?P=delimiter)"'
    r'|//(?:\\\r?\n|[^\n])*|/\*.*?\*/'
    r'|"(?:\\[\s\S]|[^"\\])*"'
    r"|'(?:\\[\s\S]|[^'\\\n])*'",
    re.DOTALL,
)
STORAGE_ALIAS = re.compile(
    r'\busing\s+(?P<name>\w+)\s*=\s*'
    r'(?P<type>hostAoSDeviceSoA\s*<[^;]+>)\s*;'
)
IDENTIFIER = re.compile(r'\b[A-Za-z_]\w*\b')


def code_only(text):
    """Mask non-code without changing character offsets or source line numbers."""
    return NON_CODE.sub(lambda match: re.sub(r'[^\n]', ' ', match.group()), text)


def storage_aliases(root):
    """Read concrete container aliases from the project's module headers."""
    aliases = {}
    for header in sorted(root.glob('*/*.h')):
        code = code_only(header.read_text(encoding='utf-8'))
        for match in STORAGE_ALIAS.finditer(code):
            name = match.group('name')
            target = ' '.join(match.group('type').split())
            if name in aliases and aliases[name] != target:
                raise ValueError(f'Ambiguous storage alias {name} in {header}')
            aliases[name] = target
    return aliases


def expand_storage_aliases(text, aliases):
    """Expand type uses while preserving their original alias declarations."""
    code = code_only(text)
    declarations = list(STORAGE_ALIAS.finditer(code))
    declaration_index = 0
    pieces = []
    previous = 0
    for token in IDENTIFIER.finditer(code):
        while declaration_index < len(declarations) and declarations[declaration_index].end() <= token.start():
            declaration_index += 1
        if declaration_index < len(declarations) and declarations[declaration_index].start() <= token.start():
            continue
        target = aliases.get(token.group())
        if target is not None:
            pieces.extend((text[previous:token.start()], target))
            previous = token.end()
    pieces.append(text[previous:])
    return ''.join(pieces)


if __name__ == '__main__':
    source = Path(sys.argv[1])
    sys.stdout.reconfigure(encoding='utf-8')
    sys.stdout.write(expand_storage_aliases(source.read_text(encoding='utf-8'), storage_aliases(ROOT)))
