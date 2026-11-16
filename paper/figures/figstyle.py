#figstyle.py
#
#by Joe Hahn and Claude Code, joe.hahn@jmh-datasciences.com, 16 November 2026.
#
#The style the paper's figures are drawn in, as one context: matplotlib's defaults first, then seaborn's darkgrid at
#the given font scale, in DejaVu Sans. A figure drawn inside it is the same whatever its caller has set, so a notebook
#that draws a paper figure with that figure's own function shows exactly the paper's (style/testing.md T28).

import contextlib

import matplotlib.pyplot as plt
import seaborn as sns

RC = {'font.family': ['DejaVu Sans'], 'mathtext.fontset': 'dejavusans'}


@contextlib.contextmanager
def style(font_scale=1.5):
    with plt.style.context('default'), sns.axes_style('darkgrid'), \
            sns.plotting_context('notebook', font_scale=font_scale), plt.rc_context(RC):
        yield sns.color_palette('deep').as_hex()
