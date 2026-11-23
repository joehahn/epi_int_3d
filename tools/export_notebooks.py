#export_notebooks.py
#
#by Joe Hahn and Claude Code, joe.hahn@jmh-datasciences.com, 13 November 2026.
#
#export the executed test notebooks to HTML on the project's GitHub Pages site, docs/notebooks/<name>.html, so the
#paper ledger can link each claim or figure to the notebook test that backs it. The notebooks have no markdown
#(style/testing.md T3): a test is the code cell whose first line is '#test N: ...', and the export puts an anchor
#<a id="testN"></a> just before that cell, which the 'notebook' fields of notes/chain.json link to
#(bending_waves#test6), through tools/ledger_cards.py. The notebooks are exported as they are, outputs included, and not re-executed:
#execute them first (uv run jupyter nbconvert --to notebook --execute --inplace <notebook>).
#
#from the repo's top folder,
#   uv run python tools/export_notebooks.py                 #every notebook below
#   uv run python tools/export_notebooks.py bending_waves   #only those named

import os
import re
import subprocess
import sys

#paths
repo_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
docs_folder = os.path.join(repo_folder, 'docs', 'notebooks')

#the notebooks that back the paper's claims
notebooks = {
    'bending_waves': 'tests/bending_waves/bending_waves.ipynb',
    'vertical_resonance': 'tests/vertical_resonance/vertical_resonance.ipynb',
    'inclined': 'tests/inclined/inclined.ipynb',
    'lindblad_resonance': 'tests/lindblad_resonance/lindblad_resonance.ipynb',
    'density_waves': 'tests/density_waves/density_waves.ipynb',
}

names = sys.argv[1:] or list(notebooks)
os.makedirs(docs_folder, exist_ok=True)
for name in names:
    source = os.path.join(repo_folder, notebooks[name])
    cmd = ['jupyter', 'nbconvert', '--to', 'html', '--output-dir', docs_folder, '--output', name, source]
    subprocess.run(cmd, check=True, capture_output=True)
    #an anchor before each test's cell, found by its first line's comment '#test N:'
    path = os.path.join(docs_folder, name + '.html')
    page = open(path).read()
    for m in reversed(list(re.finditer(r'<pre><span></span><span class="c1">#test (\d+):', page))):
        cell = page.rfind('<div class="jp-Cell ', 0, m.start())
        page = page[:cell] + '<a id="test' + m.group(1) + '"></a>' + page[cell:]
    open(path, 'w').write(page)
    print('wrote', os.path.relpath(os.path.join(docs_folder, name + '.html'), repo_folder))
