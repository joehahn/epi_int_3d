"""Mechanical checks of paper/epi_int_3d.tex, run before the author sees a revision.

Checks: every \\ref has a \\label and no label is defined twice; every \\cite key is in paper/epi_int_3d.bib;
every \\includegraphics file exists; the equation and figure reference style of joes_papers/writing_style.md
("Equation (n)" and "Figure n" start a sentence, "Eqn. (n)" and "Fig. n" elsewhere); banned words from
that guide's Vocabulary section; and, with --compile, the TeX warnings tectonic reports (undefined
references, overfull boxes). Exit status 1 when anything is found.

Run:  uv run python tools/lint_paper.py [--compile]
"""
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
PAPER = REPO/'paper'

#(pattern, advice): words the style guide rules out, matched case-insensitively outside comments
BANNED = [(r'\btest particles?\b', '"massless particles" or "tracer particles" (writing_style.md, Vocabulary)'),
          (r'\bTherefore,', '"Consequently," (writing_style.md, Sentences and connectives)')]
CASE_SENSITIVE = {r'\bTherefore,'}


def strip_comments(line):
    """The line without its TeX comment (an unescaped %)."""
    return re.sub(r'(?<!\\)%.*', '', line)


def body_lines(tex):
    """(line number, text) for each line between \\begin{document} and \\end{document}, comments removed."""
    lines = tex.splitlines()
    start = next(i for i, l in enumerate(lines) if l.startswith(r'\begin{document}'))
    return [(i + 1, strip_comments(l)) for i, l in enumerate(lines) if i > start]


def sentence_start(text, pos):
    """Whether pos starts a sentence: preceded, apart from whitespace and '(' or '{', by '.', '?', '!', ':', a
    paragraph break, or a sectioning command's closing brace."""
    before = text[:pos].rstrip(' \t({')
    if before == '' or before.endswith('\n\n'):
        return True
    display = re.search(r'\\end\{(equation|align|eqnarray)\*?\}$', before.rstrip())
    if display:                                                    #after a display: did its math end a sentence?
        before = re.sub(r'(\\label\{[^}]*\}|\\end\{\w+\*?\}|\\nonumber|\\\\|\s)+$', '',
                        before.rstrip()[:display.start()])
        return before.endswith('.')
    if before.endswith(('\\', '~')):                               #"Eqn.\ (" itself, or a tie
        return False
    return before[-1] in '.?!:' or bool(re.search(r'\\(sub)*section\*?\{[^}]*\}\s*$', before)) \
        or before.endswith(r'\item')


def main():
    problems = []
    tex = (PAPER/'epi_int_3d.tex').read_text()
    lines = body_lines(tex)
    text = '\n'.join(l for _, l in lines)
    starts = {}                                                    #character offset of each line in text
    off = 0
    for n, l in lines:
        starts[off] = n
        off += len(l) + 1
    offsets = sorted(starts)

    def line_of(pos):
        import bisect
        return starts[offsets[bisect.bisect_right(offsets, pos) - 1]]

    labels = re.findall(r'\\label\{([^}]+)\}', tex)
    for lab in sorted({l for l in labels if labels.count(l) > 1}):
        problems.append(f'label defined more than once: {lab}')
    for m in re.finditer(r'\\(?:eq)?ref\{([^}]+)\}', text):
        if m.group(1) not in labels:
            problems.append(f'line {line_of(m.start())}: \\ref to undefined label {m.group(1)}')

    bib = (PAPER/'epi_int_3d.bib').read_text()
    keys = set(re.findall(r'@\w+\{([^,]+),', bib))
    for m in re.finditer(r'\\cite[pt]?\*?(?:\[[^]]*\])*\{([^}]+)\}', text):
        for key in m.group(1).split(','):
            if key.strip() not in keys:
                problems.append(f'line {line_of(m.start())}: citation {key.strip()} is not in epi_int_3d.bib')

    for m in re.finditer(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}', text):
        if not (PAPER/m.group(1)).exists():
            problems.append(f'line {line_of(m.start())}: missing figure file {m.group(1)}')

    #reference style: long forms start sentences, short forms elsewhere
    other = {'Equation': 'Eqn.', 'Equations': 'Eqns.', 'Figure': 'Fig.', 'Figures': 'Figs.'}
    other.update({v: k for k, v in other.items()})
    for long, short, ref in [('Equations?', 'Eqns?\\.', r'\(\\ref'), ('Figures?', 'Figs?\\.', r'\\ref')]:
        for m in re.finditer(r'\b(' + long + r')(?:~|\\ |\s)+' + ref, text):
            if not sentence_start(text, m.start()):
                problems.append(f'line {line_of(m.start())}: "{m.group(1)}" mid-sentence; use "{other[m.group(1)]}"')
        for m in re.finditer(r'(' + short + r')(?:\\ |~|\s)*' + ref, text):
            if sentence_start(text, m.start()):
                problems.append(f'line {line_of(m.start())}: sentence starts with "{m.group(1)}"; '
                                f'spell out "{other[m.group(1)]}"')

    for pattern, advice in BANNED:
        for m in re.finditer(pattern, text, 0 if pattern in CASE_SENSITIVE else re.IGNORECASE):
            problems.append(f'line {line_of(m.start())}: "{m.group(0)}": use {advice}')

    if '--compile' in sys.argv:
        run = subprocess.run(['tectonic', 'epi_int_3d.tex'], cwd=PAPER, capture_output=True, text=True)
        seen = set()
        for l in (run.stdout + run.stderr).splitlines():
            if re.search(r'undefined|Overfull|Underfull|error', l, re.IGNORECASE) and 'lineno.sty' not in l \
                    and l not in seen:
                seen.add(l)
                problems.append('tectonic: ' + l.strip())

    for p in problems:
        print(p)
    print(f'{len(problems)} problem(s)' if problems else 'no problems found')
    sys.exit(1 if problems else 0)


if __name__ == '__main__':
    main()
