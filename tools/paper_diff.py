"""The paper's changes since a git revision, marked up for review (latexdiff: deletions red, additions blue).

Writes paper/diff/epi_int_3d_diff.pdf (git-ignored) and prints the pages that carry changes, so the author
reviews those pages rather than the whole paper. Needs latexdiff (brew install latexdiff) and tectonic.

Run:  uv run python tools/paper_diff.py [revision]     (default HEAD: the uncommitted changes)
"""
import re
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
PAPER = REPO/'paper'
OUT = PAPER/'diff'
#latexdiff brackets every change with \DIF{add,del}begin ... \DIF{add,del}end; make each write "b" or "e" and its
#page number (expanded when the page ships out) to epi_int_3d_diff.difpages
PAGE_HOOKS = r'''\newwrite\difpages \immediate\openout\difpages=\jobname.difpages
''' + ''.join(r'''\let\DIF%(k)sOrig\DIF%(k)s \renewcommand{\DIF%(k)s}{\DIF%(k)sOrig\write\difpages{%(c)s\thepage}}
''' % dict(k=k, c=k[3]) for k in ['addbegin', 'addend', 'delbegin', 'delend'])


def changed_pages(log):
    """The pages spanned by the changes, from the begin/end page numbers that PAGE_HOOKS writes."""
    pages, start = set(), None
    for mark in (log.read_text().split() if log.exists() else []):
        if mark[0] == 'b':
            start = int(mark[1:])
        elif start is not None:
            pages.update(range(start, int(mark[1:]) + 1))
            start = None
    return sorted(pages)


def main():
    rev = sys.argv[1] if len(sys.argv) > 1 else 'HEAD'
    if shutil.which('latexdiff') is None:
        sys.exit('latexdiff not found: brew install latexdiff')
    OUT.mkdir(exist_ok=True)
    old = OUT/'old.tex'
    old.write_text(subprocess.run(['git', 'show', f'{rev}:paper/epi_int_3d.tex'], cwd=REPO, check=True,
                                  capture_output=True, text=True).stdout)
    diff = subprocess.run(['latexdiff', '--math-markup=whole', '--type=UNDERLINE', '--graphics-markup=none',
                           str(old), str(PAPER/'epi_int_3d.tex')], check=True, capture_output=True,
                          text=True).stdout
    diff = diff.replace(r'\begin{document}', PAGE_HOOKS + r'\begin{document}', 1)
    tex = OUT/'epi_int_3d_diff.tex'
    tex.write_text(diff)
    for f in ['aastex702.cls', 'aasjournalv7.bst', 'epi_int_3d.bib']:
        shutil.copy(PAPER/f, OUT/f)
    figs = OUT/'figures'
    if figs.is_symlink() or figs.exists():
        figs.unlink()
    figs.symlink_to(PAPER/'figures')
    run = subprocess.run(['tectonic', '--keep-intermediates', tex.name], cwd=OUT, capture_output=True, text=True)
    pdf = OUT/'epi_int_3d_diff.pdf'
    if run.returncode != 0 or not pdf.exists():
        print(run.stderr[-3000:])
        sys.exit('tectonic failed on the marked-up paper')
    print('wrote', pdf)
    pages = changed_pages(OUT/'epi_int_3d_diff.difpages')
    print('pages with changes since', rev + ':', ', '.join(map(str, pages)) if pages else 'none found')


if __name__ == '__main__':
    main()
