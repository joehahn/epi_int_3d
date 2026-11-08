#openrouter.py
#
#by Joe Hahn and Claude Code, joe.hahn@jmh-datasciences.com, 8 November 2026.
#
#what tools/review.py and tools/tiebreak.py share: the models, the OpenRouter calls, files as text for a
#prompt, and saving replies and their costs in reviews/.

import datetime
import json
import os
import subprocess
import sys
import urllib.error
import urllib.request

#paths and settings
repo_folder = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
reviews_folder = os.path.join(repo_folder, 'reviews')
openrouter_url = 'https://openrouter.ai/api/v1'
max_output_tokens = 96000
max_reasoning_tokens = 32000    #capped separately, so the reply's text always has room
max_output_characters = 3000    #longest notebook cell output sent to a reviewer
today = datetime.date.today().isoformat()

#the models, chosen 2026-11-08 from OpenRouter's list: each company's most capable model. The 2026-11-08
#bake-off made Claude the author of derivations and GPT the main math reviewer (PLAN.md, Phase 3b).
models = {
    'gpt': 'openai/gpt-6-astra',
    'gemini': 'google/gemini-4-argon',
    'claude': 'anthropic/claude-opus-5.5',
    'grok': 'x-ai/grok-4.7',
}

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
def http_json(url, body=None, key=None, timeout=1800):
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


#a notebook as text: each cell's source and its printed outputs, with plots replaced by a marker
def notebook_text(path):
    with open(path) as fp:
        notebook = json.load(fp)
    text = ''
    for n, cell in enumerate(notebook['cells']):
        text += '\n# ---- cell ' + str(n) + ' (' + cell['cell_type'] + ')\n' + ''.join(cell['source']) + '\n'
        for output in cell.get('outputs', []):
            printed = ''
            if ('text' in output):
                printed = ''.join(output['text'])
            elif ('data' in output):
                if ('image/png' in output['data']):
                    printed = '[plot]'
                elif ('text/plain' in output['data']):
                    printed = ''.join(output['data']['text/plain'])
            if (len(printed) > max_output_characters):
                printed = printed[:max_output_characters] + '\n[output truncated]'
            if (printed != ''):
                text += '# output:\n' + printed + '\n'
    return text


#one file as text, labeled with its path
def file_text(path):
    full_path = os.path.join(repo_folder, path)
    if (path.endswith('.ipynb')):
        body = notebook_text(full_path)
    else:
        with open(full_path) as fp:
            body = fp.read()
    return '\n## File: ' + path + '\n\n```\n' + body + '\n```\n'


#the changes to files since git revision rev
def diff_text(paths, rev):
    command = ['git', 'diff', rev, '--'] + paths
    diff = subprocess.run(command, cwd=repo_folder, capture_output=True, text=True, check=True).stdout
    return '\n## Changes since ' + rev + ' (git diff)\n\n```diff\n' + diff + '\n```\n'


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
def save_reply(activity, name, model, inputs, text, usage):
    base = os.path.join(reviews_folder, today + '-' + activity + '-' + name)
    path = base + '.md'
    number = 2
    while (os.path.exists(path)):
        path = base + '-' + str(number) + '.md'
        number += 1
    header = '# ' + activity + ': ' + model + '\n\n'
    header += 'Input: `' + inputs + '` · ' + today
    header += ' · tokens in ' + str(usage.get('prompt_tokens')) + ', out ' + str(usage.get('completion_tokens'))
    header += ' · cost $' + str(round(usage.get('cost', 0.0), 4)) + '\n\n---\n\n'
    with open(path, 'w') as fp:
        fp.write(header + text + '\n')
    return path


#append one line per reply to reviews/costs.jsonl; activity tags the cost for the ledger
def log_cost(activity, name, model, inputs, usage, reply_file):
    row = {'date': today, 'activity': activity, 'reviewer': name, 'model': model, 'input': inputs,
        'reply': reply_file, 'usd': round(usage.get('cost', 0.0), 4), 'tokens_in': usage.get('prompt_tokens'),
        'tokens_out': usage.get('completion_tokens')}
    with open(os.path.join(reviews_folder, 'costs.jsonl'), 'a') as fp:
        fp.write(json.dumps(row) + '\n')
    return row
