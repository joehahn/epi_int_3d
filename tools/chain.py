#chain.py
#
#by Joe Hahn and Claude Code, joe.hahn@jmh-datasciences.com, 9 November 2026.
#
#track the dependency chain, claim -> equation -> sympy check -> code -> test -> figure, numbers -> paper, and
#mark everything downstream of a change as stale. Phase 3b of PLAN.md.
#
#The chain is notes/chain.json, written by hand (by Claude): a list of nodes, each with an id, a kind,
#an anchor (the piece of a file it is), the nodes it needs, and optionally a check command that confirms
#it. When a node is confirmed, the hash of its anchor and the hashes of the nodes it needs are saved in
#notes/chain_state.json. Later a node is
#   edited   if its own anchor changed since it was confirmed,
#   stale    if a node it needs changed, or is itself edited or stale,
#   new      if it was never confirmed,
#   missing  if its anchor can't be found,
#   ok       otherwise.
#Anchors hash only a piece of a file, so an edit elsewhere in it doesn't mark the node stale:
#   path                    the whole file (a notebook's code cells only, not its outputs)
#   path::name              a Python function, or Class.method, found with ast
#   path#ID                 .md: the table row whose first cell is ID (its first two cells only, so
#                           ticking the Hand column doesn't count as an edit); .py: the file's preamble
#                           plus each '# --- ... ---' (or '#---') section whose header names ID (H1-H3
#                           names H2);
#                           .tex: the equation or figure environment holding \label{ID}, or the
#                           (sub)section that \label{ID} follows
#   path@phrase             the paragraph (text between blank lines) that contains phrase
#Whitespace is collapsed before hashing, so refilling a paragraph isn't an edit.
#
#usage, from the repo's top folder:
#   uv run python tools/chain.py                     nodes that aren't ok, and what to do about them
#   uv run python tools/chain.py --all               every node
#   uv run python tools/chain.py --tree eq:V8        what eq:V8 needs, and what needs it
#   uv run python tools/chain.py --run               run the check commands of nodes that aren't ok,
#                                                    and confirm those that pass
#   uv run python tools/chain.py --confirm code:vertical_drift paper:drift_rotation
#                                                    confirm nodes checked by reading them
#   uv run python tools/chain.py --csv nodes.csv     also write every node and its status, for tools/ledger.py

import argparse
import ast
import csv
import datetime
import hashlib
import json
import os
import re
import subprocess
import sys

#paths
repo_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
chain_file = os.path.join(repo_folder, 'notes', 'chain.json')
state_file = os.path.join(repo_folder, 'notes', 'chain_state.json')
today = datetime.date.today().isoformat()

#what confirms each kind of node, when it has no check command
how_to_confirm = {
    'claim': 'reread the claim against what it rests on',
    'equation': 'rederive, or check against its source',
    'sympy': 'run the sympy check',
    'code': 'check the code against its equations (or a --code review)',
    'test': 'run the test',
    'notebook': 'rerun the notebook and look at its outputs',
    'figure': 'rerun the figure script and look at the figure',
    'paper': 'reread the passage against what it rests on',
    'numbers': 'rerun the simulation (tools/run_simulations.py), then tools/paper_numbers.py; reread the claims that quote them',
}


#collapse whitespace, then a short sha256
def hash_text(text):
    text = ' '.join(text.split())
    return hashlib.sha256(text.encode()).hexdigest()[:12]


#the text of a whole file; a notebook's code cells only
def whole_file(path):
    with open(path) as fp:
        text = fp.read()
    if (path.endswith('.ipynb')):
        cells = json.loads(text)['cells']
        sources = [''.join(cell['source']) for cell in cells if (cell['cell_type'] == 'code')]
        text = '\n'.join(sources)
    return text


#a Python function, or Class.method
def python_function(path, name):
    text = whole_file(path)
    tree = ast.parse(text)
    parts = name.split('.')
    body = tree.body
    found = None
    for part in parts:
        found = None
        for item in body:
            if (isinstance(item, (ast.FunctionDef, ast.ClassDef)) and (item.name == part)):
                found = item
        if (found is None):
            return None
        body = found.body
    return ast.get_source_segment(text, found)


#does a header like '# --- V2: ...; V3: ... ---' or '# --- H1-H3: ... ---' name ID?
def header_names(header, ID):
    if (re.search(r'\b' + re.escape(ID) + r'\b', header) is not None):
        return True
    match = re.match(r'([A-Za-z]+)(\d+)$', ID)
    if (match is None):
        return False
    letters = match.group(1)
    number = int(match.group(2))
    for low, high in re.findall(r'\b' + letters + r'(\d+)\s*-\s*' + letters + r'?(\d+)\b', header):
        if ((int(low) <= number) and (number <= int(high))):
            return True
    return False


#the sections of a sympy script that check ID, plus the preamble they all use
def python_sections(path, ID):
    lines = whole_file(path).split('\n')
    preamble = []
    sections = []
    current = None
    for line in lines:
        if (re.match(r'#\s*---', line) is not None):
            if (header_names(line, ID)):
                current = [line]
                sections.append(current)
            else:
                current = 'skip'
        elif (current is None):
            preamble.append(line)
        elif (current != 'skip'):
            current.append(line)
    if (len(sections) == 0):
        return None
    pieces = ['\n'.join(preamble)]
    for section in sections:
        pieces.append('\n'.join(section))
    return '\n'.join(pieces)


#a markdown table row whose first cell is ID: its first two cells
def markdown_row(path, ID):
    for line in whole_file(path).split('\n'):
        cells = [cell.strip() for cell in re.split(r'(?<!\\)\|', line)]     #an escaped \| is inside a cell
        if ((len(cells) > 3) and (cells[1] == ID)):
            return cells[1] + ' | ' + cells[2]
    return None


#the LaTeX environment holding \label{ID}, or the (sub)section that the label follows
def latex_label(path, ID):
    text = whole_file(path)
    position = text.find('\\label{' + ID + '}')
    if (position < 0):
        return None
    before = text[:position]
    section_starts = [m.start() for m in re.finditer(r'\\(sub)*section\*?\{', before)]
    #a section label: from that \section to the next \section or \subsection
    if ((len(section_starts) > 0) and (before[section_starts[-1]:].count('\n') <= 1)):
        start = section_starts[-1]
        after = re.search(r'\\(sub)*section\*?\{|\\appendix|\\end\{document\}', text[position:])
        end = len(text)
        if (after is not None):
            end = position + after.start()
        return text[start:end]
    #an equation, eqnarray or figure label: its environment
    begins = [m for m in re.finditer(r'\\begin\{(equation|eqnarray|figure|table|align)\*?\}', before)]
    if (len(begins) == 0):
        return None
    begin = begins[-1]
    environment = begin.group(1)
    end = re.search(r'\\end\{' + environment + r'\*?\}', text[position:])
    if (end is None):
        return None
    return text[begin.start():position + end.end()]


#the paragraph that contains phrase
def paragraph(path, phrase):
    for piece in re.split(r'\n\s*\n', whole_file(path)):
        if (phrase in ' '.join(piece.split())):
            return piece
    return None


#the text of an anchor, or None if it can't be found
def anchor_text(anchor):
    match = re.match(r'([^:#@]+)(::|#|@)?(.*)$', anchor)
    path = os.path.join(repo_folder, match.group(1))
    kind = match.group(2)
    target = match.group(3)
    if (not os.path.exists(path)):
        return None
    if (kind is None):
        return whole_file(path)
    if (kind == '::'):
        return python_function(path, target)
    if (kind == '@'):
        return paragraph(path, target)
    if (path.endswith('.md')):
        return markdown_row(path, target)
    if (path.endswith('.py')):
        return python_sections(path, target)
    if (path.endswith('.tex')):
        return latex_label(path, target)
    return None


#hash of every node's anchor, None where it is missing
def current_hashes(nodes):
    hashes = {}
    for node in nodes:
        text = anchor_text(node['anchor'])
        hashes[node['id']] = None
        if (text is not None):
            hashes[node['id']] = hash_text(text)
    return hashes


#nodes in dependency order, needs first
def ordered(nodes):
    by_id = {node['id']: node for node in nodes}
    done = []
    seen = set()
    def visit(ID):
        if (ID in seen):
            return
        seen.add(ID)
        for need in by_id[ID]['needs']:
            visit(need)
        done.append(by_id[ID])
    for node in nodes:
        visit(node['id'])
    return done


#each node's status and the reason for it
def statuses(nodes, state, hashes):
    status = {}
    reason = {}
    for node in ordered(nodes):
        ID = node['id']
        saved = state.get(ID)
        if (hashes[ID] is None):
            status[ID] = 'missing'
            reason[ID] = 'anchor not found: ' + node['anchor']
            continue
        if (saved is None):
            status[ID] = 'new'
            reason[ID] = 'never confirmed'
            continue
        changed = [need for need in node['needs'] if (saved['needs'].get(need) != hashes[need])]
        waiting = [need for need in node['needs'] if ((need not in changed) and (status[need] != 'ok'))]
        if (saved['hash'] != hashes[ID]):
            status[ID] = 'edited'
            reason[ID] = 'its own text changed since ' + saved['date']
        elif (len(changed) > 0):
            status[ID] = 'stale'
            reason[ID] = 'changed since ' + saved['date'] + ': ' + ', '.join(changed)
        elif (len(waiting) > 0):
            status[ID] = 'stale'
            reason[ID] = 'waits on ' + ', '.join(waiting)
        else:
            status[ID] = 'ok'
            reason[ID] = 'confirmed ' + saved['date']
    return status, reason


#record a node as confirmed now
def confirm(node, state, hashes):
    needs = {need: hashes[need] for need in node['needs']}
    state[node['id']] = {'hash': hashes[node['id']], 'needs': needs, 'date': today}


#what a node needs, and what needs it, as indented trees
def print_tree(nodes, ID, status):
    by_id = {node['id']: node for node in nodes}
    needed_by = {node['id']: [] for node in nodes}
    for node in nodes:
        for need in node['needs']:
            needed_by[need].append(node['id'])
    def walk(ID, links, depth):
        print('   '*depth + ID + '   (' + status[ID] + ')   ' + by_id[ID]['anchor'])
        for other in links[ID]:
            walk(other, links, depth + 1)
    needs = {node['id']: node['needs'] for node in nodes}
    print('what', ID, 'needs:')
    walk(ID, needs, 1)
    print('what needs', ID + ':')
    walk(ID, needed_by, 1)


#parse the command line
parser = argparse.ArgumentParser()
parser.add_argument('--all', action='store_true', help='list every node, not only those that are not ok')
parser.add_argument('--tree', type=str, default=None, help='show what a node needs and what needs it')
parser.add_argument('--run', action='store_true', help='run the checks of nodes that are not ok; confirm those that pass')
parser.add_argument('--confirm', type=str, nargs='+', default=None, help='confirm these nodes')
parser.add_argument('--csv', type=str, default=None, help='write every node and its status to this file')
args = parser.parse_args()

#read the chain and its state, and hash every anchor
with open(chain_file) as fp:
    nodes = json.load(fp)['nodes']
state = {}
if (os.path.exists(state_file)):
    with open(state_file) as fp:
        state = json.load(fp)
by_id = {node['id']: node for node in nodes}
for node in nodes:
    for need in node['needs']:
        if (need not in by_id):
            print('unknown need', need, 'in', node['id'])
            sys.exit(1)
hashes = current_hashes(nodes)
status, reason = statuses(nodes, state, hashes)

#show one node's tree
if (args.tree is not None):
    print_tree(nodes, args.tree, status)
    sys.exit(0)

#confirm nodes by hand, needs first; a node whose needs aren't ok is left alone
if (args.confirm is not None):
    for node in ordered(nodes):
        if (node['id'] in args.confirm):
            blocked = [need for need in node['needs'] if (status[need] != 'ok')]
            if (len(blocked) > 0):
                print('not confirmed', node['id'], '   first confirm:', ', '.join(blocked))
            elif (hashes[node['id']] is None):
                print('not confirmed', node['id'], '   anchor not found')
            else:
                confirm(node, state, hashes)
                status[node['id']] = 'ok'
                print('confirmed', node['id'])

#run each distinct check command once, needs first, and confirm the nodes whose check passed
if (args.run):
    results = {}
    for node in ordered(nodes):
        ID = node['id']
        command = node.get('check')
        if ((status[ID] == 'ok') or (command is None) or (hashes[ID] is None)):
            continue
        blocked = [need for need in node['needs'] if (status[need] != 'ok')]
        if (len(blocked) > 0):
            print('skipped', ID, '   first confirm:', ', '.join(blocked))
            continue
        if (command not in results):
            print('running:', command)
            result = subprocess.run(command, shell=True, cwd=repo_folder, capture_output=True, text=True)
            results[command] = (result.returncode == 0)
            if (result.returncode != 0):
                print(result.stdout[-2000:] + result.stderr[-2000:])
        if (results[command]):
            confirm(node, state, hashes)
            status[ID] = 'ok'
            print('   passed, confirmed', ID)
        else:
            print('   FAILED', ID)

#save the state
if ((args.confirm is not None) or args.run):
    with open(state_file, 'w') as fp:
        json.dump(state, fp, indent=1, sort_keys=True)
    status, reason = statuses(nodes, state, hashes)

#report
counts = {}
for node in nodes:
    counts[status[node['id']]] = counts.get(status[node['id']], 0) + 1
print('nodes =', len(nodes), '  ', '   '.join(key + ' = ' + str(counts[key]) for key in sorted(counts)))
for node in ordered(nodes):
    ID = node['id']
    if ((status[ID] == 'ok') and (not args.all)):
        continue
    todo = ''
    if (status[ID] != 'ok'):
        todo = '   -> ' + (node.get('check') or how_to_confirm[node['kind']])
    print(status[ID].ljust(8), ID.ljust(28), reason[ID] + todo)

#every node and its status, for the ledger
if (args.csv is not None):
    needed_by = {node['id']: [] for node in nodes}
    for node in nodes:
        for need in node['needs']:
            needed_by[need].append(node['id'])
    with open(args.csv, 'w', newline='') as fp:
        writer = csv.writer(fp)
        writer.writerow(['id', 'kind', 'status', 'reason', 'anchor', 'needs', 'needed_by', 'note', 'check', 'confirmed'])
        for node in ordered(nodes):
            ID = node['id']
            confirmed = ''
            if (ID in state):
                confirmed = state[ID]['date']
            writer.writerow([ID, node['kind'], status[ID], reason[ID], node['anchor'], ' '.join(node['needs']),
                ' '.join(needed_by[ID]), node.get('note', ''), node.get('check', ''), confirmed])
