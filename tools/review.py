#review.py
#
#by Joe Hahn and Claude Code, joe.hahn@jmh-datasciences.com, 8 November 2026.
#
#send a problem, code, a test or the paper to outside AI models via OpenRouter, and save their raw
#replies in reviews/. Adapted from agent-banking-watch's scripts/review.py. Phase 3b of PLAN.md.
#
#usage, from the repo's top folder:
#   uv run python tools/review.py --derive reviews/problems/vertical_motion.md
#   uv run python tools/review.py --code fast/bookkeeping.py --context notes/equations.md
#   uv run python tools/review.py --code fast/bookkeeping.py --diff HEAD~3 --context notes/equations.md
#   uv run python tools/review.py --test tests/test_vertical.py --context tests/vertical_analysis.py fast/bookkeeping.py
#   uv run python tools/review.py --paper paper/epi_int_3d.tex --context notes/equations.md --models gpt,gemini,grok
#add --dry-run to see the models, the prompt's size and the most it could cost, without sending anything;
#--focus "question" adds a question every reviewer must answer first; --models picks reviewers by name.
#
#the four modes, one per activity in the ledger:
#   --derive    a blind re-derivation: the problem statement only, never Claude's answer
#   --code      the code (or its diff since --diff) against the equations it implements (--context)
#   --test      would this test fail if the code were wrong? (code under test in --context)
#   --paper     every claim against the result behind it, plus a skeptical referee's read
#each reply is saved unedited as reviews/<date>-<activity>-<reviewer>.md and its cost appended to
#reviews/costs.jsonl. Claude then records a verdict on each finding (reviews/<date>-<activity>-verdicts.json),
#and tools/tiebreak.py sends the disputed ones to a model not involved in them.
#The OpenRouter key is read from .env (OPENROUTER_API_KEY=...), which is git-ignored.

import argparse
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from openrouter import (repo_folder, reviews_folder, today, models, max_output_tokens, read_key, model_prices,
    file_text, diff_text, ask, save_reply, log_cost)

#default reviewers per activity; Claude is left out, since Claude wrote what the others review
default_reviewers = {'derive': 'gpt,gemini', 'code': 'gpt,gemini', 'test': 'gpt,gemini', 'paper': 'gpt,gemini'}

#style guides sent along with each activity
style_files = {
    'derive': ['style/derivations.md'],
    'code': ['style/code.md'],
    'test': ['style/testing.md', 'style/dataviz.md'],
    'paper': ['joes_papers/writing_style.md', 'style/dataviz.md'],
}

#what every reviewer of code, tests or the paper is told about the project
project_note = '''The project is epi_int_3d, an N-body code for narrow planetary rings orbiting an oblate planet
(mass M_p, radius R_p, second zonal harmonic J_2), written by Joe Hahn with Claude Code (Anthropic) as the
single author of code, derivations and drafts. It extends epi_int_lite (Hahn et al. 2025, ApJ 986, 220), whose
particles move in the planet's equatorial plane, into the vertical direction, to simulate linear spiral
bending waves (and later linear density waves) launched by a satellite. Each particle is advanced by a drift
step (motion under the planet's gravity alone, done exactly through epicyclic orbit elements) and kick steps
(ring gravity, viscosity, the satellite). Code units: G = 1, M_p = 1, the ring's inner edge at r = 1.
Notation follows Goldreich & Tremaine: Omega_a, kappa_a, nu_a are the angular, epicyclic and vertical
frequencies at semimajor axis a; Omega (no subscript) is the longitude of ascending node.'''

#the findings format shared by the code, test and paper reviews, so that the verdicts can cite them
findings_format = '''Return numbered findings, most important first. For each finding give:
- **Type:** {types}
- **Severity:** high | medium | low
- **Where:** the file and line, function or section, and a short exact quote
- **Finding:** what is wrong or could be better, in one to three sentences
- **Suggestion:** the concrete change, or "none"
Flag a breach of the author's style guide (given after the material) as Type: style, citing its rule
number (e.g. C12, T7, P31). Don't rewrite the material in your own style. If you find nothing wrong of a
kind, say so; don't pad.'''

#instructions for each activity
instructions = {}
instructions['derive'] = '''You are a theoretical astrophysicist working on planetary ring dynamics.
Derive the result asked for below from scratch. Show every step, state every approximation and the
order of the terms it neglects, and check your result in its limits. Write the derivation as a
draft appendix of a journal paper (ApJ), following the author's derivation style guide, which is
given after the problem. Use plain text with LaTeX math ($...$), and number the displayed equations.
End with a short list of the final results.'''
instructions['code'] = '''You are reviewing scientific code before it is used in simulations for a journal paper.
''' + project_note + '''

The code to review comes first, then the equations or derivations it implements (the context). Check:
1. Does the code implement those equations exactly: every term, factor, sign and index?
2. Units and normalization: are G, M_p, R_p and J_2 used consistently with code units?
3. Numerical issues: loss of precision, division by zero or a branch cut (atan2, sqrt of a negative
   number), angle wrapping, a quantity that should be conserved but is not.
4. Bugs of any other kind: wrong array index, off-by-one, aliasing, a variable used before it is set.
Do not suggest refactoring, type hints, docstrings, classes or input validation: the author's style
guide forbids them. For each claimed bug, give a concrete input that would expose it if you can.

''' + findings_format.replace('{types}', 'bug | equation-mismatch | units | numerics | style | other')
instructions['test'] = '''You are reviewing a test of scientific code before its results go into a journal paper.
''' + project_note + '''

The test (a pytest file or a Jupyter notebook, given as text with its printed outputs) comes first, then
the code it tests and any theory it compares with (the context). The test's purpose is to confirm that
the simulation agrees with theory, that its output agrees with expectations, or that constants of the
motion are conserved. Answer, as your first findings:
1. Name the most plausible bug in the code under test that this test would NOT catch, and say what
   check would catch it.
2. Is the tolerance, or the agreement shown in a plot, tight enough to catch a factor-of-2 error, a sign
   error, or a missing term of the expected order? Say what the tolerance actually constrains.
3. Is the theory the test compares with the right theory for this setup, and are its assumptions met?
Then anything else: a comparison that is circular (the test checks the code against itself), an output
that is printed but never checked, a plot that does not make its point.

''' + findings_format.replace('{types}', 'blind-spot | tolerance | theory | circular | plot | style | other')
instructions['paper'] = '''You are a skeptical referee for The Astrophysical Journal, reviewing a draft paper.
''' + project_note + '''

The paper's LaTeX source comes first, then the notes and equation tracker its results rest on (the
context). Check every claim against the result behind it: does each number, equation and statement of
agreement with theory follow from the derivation, simulation or citation it rests on? Flag claims that
are unsupported, overstated, or contradicted by the context. Then add three short sections: "Strongest
parts", "Weakest parts", and "What a skeptical referee would push back on". You cannot check citations
against the literature here; flag a citation only if it looks wrong for the claim it supports.

''' + findings_format.replace('{types}', 'unsupported | overstated | error | unclear | missing | style | other')

#the question the author adds with --focus
focus_note = 'The author asks you specifically, and wants it answered as your first finding: '


#assemble the prompt for one activity
def make_prompt(activity, targets, context_files, rev, focus):
    prompt = instructions[activity]
    if (focus is not None):
        prompt += '\n\n' + focus_note + focus
    if (activity == 'derive'):
        prompt += '\n\n# The problem\n\n'
        with open(os.path.join(repo_folder, targets[0])) as fp:
            prompt += fp.read()
    else:
        prompt += '\n\n# The material to review\n'
        if (rev is not None):
            prompt += diff_text(targets, rev)
        for path in targets:
            prompt += file_text(path)
        if (len(context_files) > 0):
            prompt += '\n\n# Context\n'
            for path in context_files:
                prompt += file_text(path)
    for style_file in style_files[activity]:
        with open(os.path.join(repo_folder, style_file)) as fp:
            prompt += '\n\n# The author\'s style guide (' + style_file + ')\n\n' + fp.read()
    return prompt


#parse the command line
parser = argparse.ArgumentParser()
mode = parser.add_mutually_exclusive_group(required=True)
mode.add_argument('--derive', type=str, help='problem statement for a blind derivation')
mode.add_argument('--code', type=str, nargs='+', help='code files to review')
mode.add_argument('--test', type=str, nargs='+', help='test files or notebooks to review')
mode.add_argument('--paper', type=str, nargs='+', help='the paper\'s LaTeX source')
parser.add_argument('--context', type=str, nargs='*', default=[], help='equations, derivations or code the reviewer needs')
parser.add_argument('--diff', type=str, default=None, help='also send the changes since this git revision')
parser.add_argument('--focus', type=str, default=None, help='a question every reviewer answers first')
parser.add_argument('--models', type=str, default=None, help='comma-separated names from models')
parser.add_argument('--dry-run', action='store_true', help='show models, prompt size and cost estimate only')
args = parser.parse_args()
for activity in ['derive', 'code', 'test', 'paper']:
    if (getattr(args, activity) is not None):
        targets = getattr(args, activity)
        break
if (activity == 'derive'):
    targets = [targets]
chosen_names = args.models
if (chosen_names is None):
    chosen_names = default_reviewers[activity]
chosen_names = chosen_names.split(',')
inputs = ' '.join(targets)
if (args.diff is not None):
    inputs += ' (diff since ' + args.diff + ')'

#prompt, models and a cost estimate, assuming the longest possible reply
prompt = make_prompt(activity, targets, args.context, args.diff, args.focus)
prices = model_prices(chosen_names)
prompt_tokens_estimate = len(prompt)//4
print('activity =', activity, '   input =', inputs)
print('context =', ' '.join(args.context))
print('prompt length (characters) =', len(prompt), '   (about', prompt_tokens_estimate, 'tokens)')
total_cost_max = 0.0
for name in chosen_names:
    price_in, price_out = prices[name]
    cost_max = prompt_tokens_estimate*price_in + max_output_tokens*price_out
    total_cost_max += cost_max
    print(name, '=', models[name], '   most it could cost ($) =', round(cost_max, 2))
print('most the run could cost ($) =', round(total_cost_max, 2))
if (args.dry_run):
    sys.exit(0)

#send the prompt to each model in turn
key = read_key('OPENROUTER_API_KEY')
os.makedirs(reviews_folder, exist_ok=True)
total_cost = 0.0
for name in chosen_names:
    model = models[name]
    print('asking', name, '=', model, '...')
    text, usage = ask(model, prompt, key)
    if (text is None):
        continue
    path = save_reply(activity, name, model, inputs, text, usage)
    row = log_cost(activity, name, model, inputs, usage, os.path.relpath(path, repo_folder))
    total_cost += row['usd']
    print('   saved', os.path.relpath(path, repo_folder), '   cost ($) =', row['usd'], '   tokens out =', row['tokens_out'])
print('total cost ($) =', round(total_cost, 4))
