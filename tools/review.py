#review.py
#
#by Joe Hahn and Claude Code, joe.hahn@jmh-datasciences.com, 8 November 2026.
#
#send a problem or a draft to outside AI models via OpenRouter, and save their raw replies in reviews/.
#Adapted from agent-banking-watch's scripts/review.py. Phase 3b of PLAN.md.
#
#usage, from the repo's top folder:
#   uv run python tools/review.py --derive reviews/problems/vertical_motion.md --dry-run
#   uv run python tools/review.py --derive reviews/problems/vertical_motion.md
#   uv run python tools/review.py --derive reviews/problems/vertical_motion.md --models gpt,gemini
#--derive sends a problem statement to each model, which derives the result from scratch (a blind
#re-derivation). --dry-run shows the models, the prompt's size and a cost estimate, and sends nothing.
#The OpenRouter key is read from .env (OPENROUTER_API_KEY=...), which is git-ignored. Each reply is
#saved unedited as reviews/<date>-<activity>-<model>.md, and its cost appended to reviews/costs.jsonl.

import argparse
import datetime
import json
import os
import sys
import urllib.error
import urllib.request

#paths and settings
repo_folder = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
reviews_folder = os.path.join(repo_folder, 'reviews')
openrouter_url = 'https://openrouter.ai/api/v1'
max_output_tokens = 96000
max_reasoning_tokens = 32000    #capped separately, so the reply's text always has room
today = datetime.date.today().isoformat()

#the models, chosen 2026-11-08 from OpenRouter's list: each company's most capable model for math
models = {
    'gpt': 'openai/gpt-6-astra',
    'gemini': 'google/gemini-4-argon',
    'claude': 'anthropic/claude-opus-5.5',
}

#style guides sent along with a derivation problem, so that each model writes in the author's style
derivation_style_files = ['style/derivations.md']

#instructions that precede a derivation problem
derive_instructions = '''You are a theoretical astrophysicist working on planetary ring dynamics.
Derive the result asked for below from scratch. Show every step, state every approximation and the
order of the terms it neglects, and check your result in its limits. Write the derivation as a
draft appendix of a journal paper (ApJ), following the author's derivation style guide, which is
given after the problem. Use plain text with LaTeX math ($...$), and number the displayed equations.
End with a short list of the final results.'''


#read the OpenRouter key from .env
def read_key(name):
    key = os.environ.get(name)
    env_file = os.path.join(repo_folder, '.env')
    if ((key is None) and os.path.exists(env_file)):
        with open(env_file) as fp:
            for line in fp:
                if (line.startswith(name + '=')):
                    key = line.split('=', 1)[1].strip()
    if ((key is None) or (key == 'paste-your-key-here')):
        sys.exit(name + ' missing: put it in .env')
    return key


#send a request to url and return its json reply
def http_json(url, body=None, key=None, timeout=900):
    headers = {'Content-Type': 'application/json'}
    if (key is not None):
        headers['Authorization'] = 'Bearer ' + key
        headers['X-Title'] = 'epi_int_3d review'
    data = None
    if (body is not None):
        data = json.dumps(body).encode()
    request = urllib.request.Request(url, data=data, headers=headers)
    with urllib.request.urlopen(request, timeout=timeout) as reply:
        return json.load(reply)


#check that every chosen model is on OpenRouter, and get its prices in $ per token
def model_prices(chosen_names):
    listed = {}
    for m in http_json(openrouter_url + '/models')['data']:
        listed[m['id']] = m
    prices = {}
    for name in chosen_names:
        model = models[name]
        if (model not in listed):
            sys.exit(model + ' is not on OpenRouter; update models in tools/review.py')
        pricing = listed[model]['pricing']
        prices[name] = (float(pricing['prompt']), float(pricing['completion']))
    return prices


#assemble the prompt for a derivation problem
def derive_prompt(problem_file):
    with open(os.path.join(repo_folder, problem_file)) as fp:
        problem = fp.read()
    prompt = derive_instructions + '\n\n# The problem\n\n' + problem
    for style_file in derivation_style_files:
        with open(os.path.join(repo_folder, style_file)) as fp:
            prompt += '\n\n# The author\'s style guide (' + style_file + ')\n\n' + fp.read()
    return prompt


#send one prompt to one model, return the reply's text and OpenRouter's usage record; an empty or
#truncated reply is still returned, flagged, so that its cost is logged
def ask(model, prompt, key):
    body = {'model': model, 'messages': [{'role': 'user', 'content': prompt}],
        'max_tokens': max_output_tokens, 'reasoning': {'max_tokens': max_reasoning_tokens},
        'usage': {'include': True}}
    try:
        reply = http_json(openrouter_url + '/chat/completions', body=body, key=key)
    except urllib.error.HTTPError as error:
        print('OpenRouter HTTP error', error.code, 'for', model, ':', error.read().decode()[:400])
        return None, None
    choice = reply['choices'][0]
    text = choice['message'].get('content')
    if (text is None):
        text = ''
    usage = reply.get('usage', {})
    finish_reason = choice.get('finish_reason')
    if ((finish_reason != 'stop') or (text == '')):
        print('   warning: finish_reason =', finish_reason, '   reply length (characters) =', len(text))
        text = '*(warning: the reply ended with finish_reason = ' + str(finish_reason) + ')*\n\n' + text
    return text, usage


#save a reply unedited, without overwriting an earlier one from the same day
def save_reply(activity, name, model, problem_file, text, usage):
    base = os.path.join(reviews_folder, today + '-' + activity + '-' + name)
    path = base + '.md'
    number = 2
    while (os.path.exists(path)):
        path = base + '-' + str(number) + '.md'
        number += 1
    header = '# ' + activity + ': ' + model + '\n\n'
    header += 'Problem: `' + problem_file + '` · ' + today
    header += ' · tokens in ' + str(usage.get('prompt_tokens')) + ', out ' + str(usage.get('completion_tokens'))
    header += ' · cost $' + str(round(usage.get('cost', 0.0), 4)) + '\n\n---\n\n'
    with open(path, 'w') as fp:
        fp.write(header + text + '\n')
    return path


#append one line per reply to reviews/costs.jsonl; activity tags the cost for the ledger
def log_cost(activity, name, model, problem_file, usage):
    row = {'date': today, 'activity': activity, 'reviewer': name, 'model': model, 'input': problem_file,
        'usd': round(usage.get('cost', 0.0), 4), 'tokens_in': usage.get('prompt_tokens'),
        'tokens_out': usage.get('completion_tokens')}
    with open(os.path.join(reviews_folder, 'costs.jsonl'), 'a') as fp:
        fp.write(json.dumps(row) + '\n')
    return row


#parse the command line
parser = argparse.ArgumentParser()
parser.add_argument('--derive', type=str, required=True, help='problem statement for a blind derivation')
parser.add_argument('--models', type=str, default='gpt,gemini,claude', help='comma-separated names from models')
parser.add_argument('--dry-run', action='store_true', help='show models, prompt size and cost estimate only')
args = parser.parse_args()
activity = 'derive'
chosen_names = args.models.split(',')

#prompt, models and a cost estimate, assuming the longest possible reply
prompt = derive_prompt(args.derive)
prices = model_prices(chosen_names)
prompt_tokens_estimate = len(prompt)//4
print('problem =', args.derive)
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
    path = save_reply(activity, name, model, args.derive, text, usage)
    row = log_cost(activity, name, model, args.derive, usage)
    total_cost += row['usd']
    print('   saved', os.path.relpath(path, repo_folder), '   cost ($) =', row['usd'], '   tokens out =', row['tokens_out'])
print('total cost ($) =', round(total_cost, 4))
