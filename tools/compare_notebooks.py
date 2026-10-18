#!/usr/bin/env python3
"""Compare epi_int_lite's test notebooks as run on Python 2 and as re-executed on the port.

The Python 2 notebooks, with their original outputs, are in git commit 99a7753; the working
copies are their Python 3 ports, executed on the port's output with
    uv run jupyter nbconvert --to notebook --execute --inplace <notebook>
For each notebook this places every figure beside the same cell's Python 2 figure, writing
notes/phase0/figs/<notebook>_cellNN.png, and prints each cell's printed output from both runs.

Run:  uv run python tools/compare_notebooks.py
"""
import base64
import io
import json
import subprocess
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parent.parent
PY2_COMMIT = '99a7753'
FIGS = REPO/'notes'/'phase0'/'figs'
NOTEBOOKS = ['tests/unperturbed/check_orbits.ipynb', 'tests/spreading/spreading.ipynb',
             'tests/libration/narrow_ring/ring_mass_1.0e-9/libration.ipynb',
             'tests/libration/libration_frequency.ipynb']


def outputs(nb):
    """Map each code cell's index to its printed text and its figures."""
    out = {}
    for i, c in enumerate(nb['cells']):
        if c['cell_type'] != 'code':
            continue
        text, images = [], []
        for o in c.get('outputs', []):
            if 'text' in o:
                text.append(''.join(o['text']).rstrip())
            elif 'data' in o:
                if 'image/png' in o['data']:
                    images.append(base64.b64decode(o['data']['image/png']))
                elif 'text/plain' in o['data']:
                    text.append(''.join(o['data']['text/plain']).rstrip())
        out[i] = ('\n'.join(text), images)
    return out


def side_by_side(png2, png3, title, path):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for ax, png, label in zip(axes, (png2, png3), ('Python 2 (epi_int_lite)', 'Python 3 (this port)')):
        ax.imshow(mpimg.imread(io.BytesIO(png), format='png'))
        ax.set_axis_off()
        ax.set_title(label, fontsize=11)
    fig.suptitle(title, fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def main():
    FIGS.mkdir(parents=True, exist_ok=True)
    for f in NOTEBOOKS:
        py2 = json.loads(subprocess.run(['git', 'show', f'{PY2_COMMIT}:{f}'], cwd=REPO,
                                        capture_output=True, text=True, check=True).stdout)
        py3 = json.load(open(REPO/f))
        o2, o3 = outputs(py2), outputs(py3)
        name = Path(f).stem
        print(f'\n######## {f}')
        for i in sorted(o3):
            text2, images2 = o2.get(i, ('', []))
            text3, images3 = o3[i]
            if text2 or text3:
                print(f'--- cell {i}\n  Python 2: ' + text2.replace('\n', '\n            ') +
                      '\n  Python 3: ' + text3.replace('\n', '\n            '))
            for k, (png2, png3) in enumerate(zip(images2, images3)):
                suffix = f'_{k}' if len(images3) > 1 else ''
                side_by_side(png2, png3, f'{Path(f).name}, cell {i}',
                             FIGS/f'{name}_cell{i:02d}{suffix}.png')
            if len(images2) != len(images3):
                print(f'  cell {i}: {len(images2)} Python 2 figures, {len(images3)} Python 3 figures')
    print(f'\nfigures written to {FIGS.relative_to(REPO)}/')


if __name__ == '__main__':
    main()
