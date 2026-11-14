#export_notebooks.py
#
#by Joe Hahn and Claude Code, joe.hahn@jmh-datasciences.com, 13 November 2026.
#
#export the executed test notebooks to HTML on the project's GitHub Pages site, docs/notebooks/<name>.html, so the
#paper ledger can link each claim or figure to the notebook section that tests it (a section's anchor is the
#<a id="testN"></a> at the top of its heading cell; the links come from the 'notebook' fields of notes/chain.json,
#written by tools/ledger_cards.py). The notebooks are exported as they are, outputs included, and not re-executed:
#execute them first (uv run jupyter nbconvert --to notebook --execute --inplace <notebook>).
#
#from the repo's top folder,
#   uv run python tools/export_notebooks.py                 #every notebook below
#   uv run python tools/export_notebooks.py bending_waves   #only those named

import os
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
}

names = sys.argv[1:] or list(notebooks)
os.makedirs(docs_folder, exist_ok=True)
for name in names:
    source = os.path.join(repo_folder, notebooks[name])
    cmd = ['jupyter', 'nbconvert', '--to', 'html', '--output-dir', docs_folder, '--output', name, source]
    subprocess.run(cmd, check=True, capture_output=True)
    print('wrote', os.path.relpath(os.path.join(docs_folder, name + '.html'), repo_folder))
