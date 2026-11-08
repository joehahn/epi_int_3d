#tiebreak.py
#
#by Joe Hahn and Claude Code, joe.hahn@jmh-datasciences.com, 8 November 2026.
#
#tie-break the review findings that Claude rejected or only partly accepted. Claude wrote the code, tests
#and paper and then judged the reviews of them, so these are the verdicts most open to self-bias. Each
#disputed finding goes to a model that is neither Claude nor the reviewer who made it (GPT's to Gemini,
#everyone else's to GPT), which rules from the material. Adapted from agent-banking-watch's
#scripts/tiebreak.py. Phase 3b of PLAN.md.
#
#usage, from the repo's top folder:
#   uv run python tools/tiebreak.py reviews/2026-11-08-code-verdicts.json --dry-run
#   uv run python tools/tiebreak.py reviews/2026-11-08-code-verdicts.json
#
#the verdicts file is written by Claude after reading the reviews:
#   {"activity": "code", "inputs": ["fast/bookkeeping.py"], "context": ["notes/equations.md"],
#    "findings": [{"id": "gpt-3", "reviewer": "gpt", "reply": "reviews/2026-11-08-code-gpt.md",
#                  "severity": "high", "where": "...", "finding": "...",
#                  "verdict": "agree | partly | disagree", "comment": "Claude's reason", "fix": "what was changed"}]}
#tiebreak.py adds each ruling to its finding in that file, writes a readable reviews/<date>-tiebreak-<activity>.md,
#and logs its cost in reviews/costs.jsonl under activity "tiebreak". Editorial decisions stay the author's.

import argparse
import json
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from openrouter import (repo_folder, reviews_folder, today, models, max_output_tokens, read_key, model_prices,
    file_text, ask, log_cost)

#who judges a finding: never Claude, never the reviewer who made it
judge_of = {'gpt': 'gemini', 'gemini': 'gpt', 'grok': 'gpt'}

#instructions to the judge
judge_instructions = '''You are breaking ties in the review of a scientific project, epi_int_3d, an N-body code for
planetary rings (Joe Hahn, with Claude Code from Anthropic as the single author of its code, derivations
and drafts). Outside models reviewed the material below, and Claude then judged their findings. These are
only the findings Claude rejected or accepted in part, because an author judging critiques of its own work
may be biased toward it.

For each item decide, from the material given, who is right:
- "reviewer": the finding is right and should be accepted in full;
- "claude": Claude's verdict is right;
- "split": both have a point; say exactly what to change;
- "need more evidence": the material can't settle it; say what would (a test, a derivation, a run).
Judge on physics and correctness. Style findings are judged against the author's style guide, not your
own taste. Be independent: some findings may have come from a model like you; judge them like any other.

Return only a JSON array, one object per item, in order:
[{"id": "<item id>", "ruling": "reviewer|claude|split|need more evidence",
  "reason": "<one to three sentences>", "fix": "<the change to make, or empty>"}]'''


#the findings Claude disputed
def disputed(verdicts):
    items = []
    for f in verdicts['findings']:
        if (f['verdict'] in ['partly', 'disagree']):
            items.append(f)
    return items


#the prompt for one judge: instructions, the disputed items, then the material they are about
def judge_prompt(verdicts, items):
    payload = []
    for f in items:
        payload.append({'id': f['id'], 'reviewer': f['reviewer'], 'severity': f.get('severity'),
            'where': f.get('where'), 'finding': f['finding'], 'claude_verdict': f['verdict'],
            'claude_comment': f['comment']})
    prompt = judge_instructions + '\n\n# The disputed findings\n\n' + json.dumps(payload, indent=1)
    prompt += '\n\n# The material reviewed\n'
    for path in verdicts['inputs']:
        prompt += file_text(path)
    if (len(verdicts.get('context', [])) > 0):
        prompt += '\n\n# Context\n'
        for path in verdicts['context']:
            prompt += file_text(path)
    return prompt


#the JSON array in a judge's reply
def parse_rulings(text):
    start = text.find('[')
    end = text.rfind(']')
    if ((start < 0) or (end < start)):
        return []
    return json.loads(text[start:end + 1])


#parse the command line
parser = argparse.ArgumentParser()
parser.add_argument('verdicts', type=str, help='reviews/<date>-<activity>-verdicts.json')
parser.add_argument('--dry-run', action='store_true', help='list the disputed findings and the cost estimate only')
args = parser.parse_args()
verdicts_file = os.path.join(repo_folder, args.verdicts)
with open(verdicts_file) as fp:
    verdicts = json.load(fp)
activity = verdicts['activity']

#group the disputed findings by judge
items = disputed(verdicts)
by_judge = {}
for f in items:
    judge = judge_of[f['reviewer']]
    if (judge not in by_judge):
        by_judge[judge] = []
    by_judge[judge].append(f)
print('activity =', activity, '   disputed findings =', len(items), 'of', len(verdicts['findings']))
for f in items:
    print('  ', f['id'], '(' + f['verdict'] + ', judge ' + judge_of[f['reviewer']] + '):', f['finding'][:80])
if (len(items) == 0):
    sys.exit(0)

#a cost estimate, assuming the longest possible reply
prices = model_prices(list(by_judge.keys()))
prompts = {}
total_cost_max = 0.0
for judge in by_judge:
    prompts[judge] = judge_prompt(verdicts, by_judge[judge])
    price_in, price_out = prices[judge]
    cost_max = (len(prompts[judge])//4)*price_in + max_output_tokens*price_out
    total_cost_max += cost_max
    print(judge, '=', models[judge], '   items =', len(by_judge[judge]), '   most it could cost ($) =', round(cost_max, 2))
print('most the run could cost ($) =', round(total_cost_max, 2))
if (args.dry_run):
    sys.exit(0)

#ask each judge, and add its rulings to the verdicts file
key = read_key('OPENROUTER_API_KEY')
rulings = {}
total_cost = 0.0
for judge in by_judge:
    print('asking', judge, '=', models[judge], '...')
    text, usage = ask(models[judge], prompts[judge], key)
    if (text is None):
        continue
    for r in parse_rulings(text):
        r['judge'] = judge
        rulings[r['id']] = r
    row = log_cost('tiebreak', judge, models[judge], args.verdicts, usage, args.verdicts)
    total_cost += row['usd']
    print('   cost ($) =', row['usd'])
for f in verdicts['findings']:
    if (f['id'] in rulings):
        f['tiebreak'] = rulings[f['id']]
with open(verdicts_file, 'w') as fp:
    json.dump(verdicts, fp, indent=1)

#a readable record
record_file = os.path.join(reviews_folder, today + '-tiebreak-' + activity + '.md')
with open(record_file, 'w') as fp:
    fp.write('# Tie-break: ' + activity + '\n\n' + today + ' · ' + str(len(items)) + ' disputed findings · ')
    fp.write('verdicts `' + args.verdicts + '` · cost $' + str(round(total_cost, 4)) + '\n\n')
    for f in items:
        r = f.get('tiebreak', {})
        fp.write('## ' + f['id'] + ': ' + r.get('ruling', 'no ruling') + ' (judge ' + r.get('judge', '?') + ')\n\n')
        fp.write('- **Finding (' + f['reviewer'] + '):** ' + f['finding'] + '\n')
        fp.write('- **Claude (' + f['verdict'] + '):** ' + f['comment'] + '\n')
        fp.write('- **Reason:** ' + r.get('reason', '') + '\n- **Fix:** ' + r.get('fix', '') + '\n\n')
print('wrote', os.path.relpath(record_file, repo_folder), '   total cost ($) =', round(total_cost, 4))
