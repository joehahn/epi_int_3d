#!/usr/bin/env python3
"""Builds paper/epi_int_3d.bib, the paper's own bibliography, from the author's bibliography.

The author's paper/jmh_bibliography.bib is the master copy; it stays local and is never committed.
This script copies out the entries that paper/epi_int_3d.tex cites, keeping the comment lines just
above each entry (the SciX-check notes) and dropping the personal keywords field. It stops if a
cited key is missing. Rerun it whenever the paper's citations change, and commit the result with
the paper.

Usage:  uv run python tools/make_bib.py
"""

import re
from pathlib import Path

PAPER = Path(__file__).resolve().parent.parent/'paper'
TEX, MASTER, OUT = PAPER/'epi_int_3d.tex', PAPER/'jmh_bibliography.bib', PAPER/'epi_int_3d.bib'
DROP = ('keywords',)


def cited_keys(tex):
    """Keys of every \\cite-like command, in order of first citation; % comments are skipped."""
    text = re.sub(r'(?<!\\)%.*', '', tex)
    keys = []
    for group in re.findall(r'\\cite[a-zA-Z]*\*?(?:\[[^\]]*\])*\{([^}]*)\}', text):
        for key in group.split(','):
            if key.strip() and key.strip() not in keys:
                keys.append(key.strip())
    return keys


def entries(bib):
    """{key (lower case): text} of each @entry, with the comment lines just above it."""
    found = {}
    for m in re.finditer(r'@\w+\s*\{\s*([^,\s]+)\s*,', bib):
        depth, end = 0, m.end() - 1
        for end in range(m.start(), len(bib)):
            depth += {'{': 1, '}': -1}.get(bib[end], 0)
            if depth == 0 and bib[end] == '}':
                break
        before = bib[:m.start()].rstrip('\n').split('\n')
        notes = []
        while before and before[-1].startswith('%'):
            notes.insert(0, before.pop())
        found[m.group(1).lower()] = '\n'.join(notes + [bib[m.start():end + 1]])
    return found


def strip_fields(entry, drop=DROP):
    """Removes the named fields, whose values may span lines."""
    for field in drop:
        m = re.search(r'\n\s*' + field + r'\s*=\s*', entry, re.IGNORECASE)
        if not m:
            continue
        depth, i = 0, m.end()
        while i < len(entry):
            c = entry[i]
            depth += {'{': 1, '}': -1}.get(c, 0)
            if depth == 0 and c in ',\n' and i > m.end():
                break
            i += 1
        entry = entry[:m.start()] + entry[i + 1 if entry[i] == ',' else i:]
    return entry


def main():
    keys = cited_keys(TEX.read_text())
    master = entries(MASTER.read_text(errors='replace'))
    missing = [k for k in keys if k.lower() not in master]
    if missing:
        raise SystemExit(f'not in {MASTER.name}: {", ".join(missing)}')
    head = ('% epi_int_3d.bib: the entries cited by epi_int_3d.tex, made by tools/make_bib.py from the\n'
            '% author\'s bibliography (not in the repo). Do not edit; rerun the script.\n')
    body = [strip_fields(master[k.lower()]) for k in sorted(keys, key=str.lower)]
    OUT.write_text(head + '\n' + '\n\n'.join(body) + '\n')
    print(f'{OUT.relative_to(PAPER.parent)}: {len(keys)} entries')


if __name__ == '__main__':
    main()
