#ledger.py
#
#by Joe Hahn and Claude Code, joe.hahn@jmh-datasciences.com, 10 November 2026.
#
#build the paper ledger, the page where the author reads the paper card by card, one card per section or
#subsection, and leaves notes for Claude to act on (https://claude.ai/artifact/KpgQvmghBVn7c4nTmFWxLZ).
#Each card holds:
#   the section's text, from tools/ledger_cards.py
#   its issues, each with Claude's recommendation: paper items of the dependency chain (notes/chain.json)
#       that aren't current, problems that tools/lint_paper.py finds, and findings of the latest outside
#       review of the paper (reviews/*-paper-verdicts.json) with Claude's verdict
#   its facts: how many chain items it holds and whether they are current, and how many of the equations
#       behind it the author has checked by hand (notes/equations.md)
#and the summary and costs: outside reviews (reviews/costs.jsonl), Claude Code at API prices and the
#author's hours per phase (tools/session_time.py), and what the author paid (the subscription).
#
#writes one JSON file per database document in the git-ignored notes/ledger/db/, which Claude uploads to
#the page: cards/<id>.json, meta/summary.json and meta/costs.json. The figures are uploaded to the page
#once, and again when they change; notes/ledger/assets.json records their ids.
#
#usage, from the repo's top folder:
#   uv run python tools/ledger.py

import csv
import datetime
import glob
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ledger_cards import cards

#paths
repo_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ledger_folder = os.path.join(repo_folder, 'notes', 'ledger')
db_folder = os.path.join(ledger_folder, 'db')
assets_file = os.path.join(ledger_folder, 'assets.json')
figures_folder = os.path.join(repo_folder, 'paper', 'figures')

#what Joe pays for Claude Code: a monthly subscription, assumed billed from the project's start
#(the first row of notes/phases.csv), one payment at the start of each month begun
subscription_usd_per_month = 100.0


#the Hand column of notes/equations.md, by equation id: yes, no or n/a
def hand_checks():
    checks = {}
    for line in open(os.path.join(repo_folder, 'notes', 'equations.md')):
        cells = [cell.strip() for cell in re.split(r'(?<!\\)\|', line)]     #an escaped \| is inside a cell
        if ((len(cells) > 7) and (cells[1] != 'ID') and (not cells[1].startswith('-'))):
            hand = cells[6]
            if ('✅' in hand):
                hand = 'yes'
            elif ('☐' in hand):
                hand = 'no'
            checks[cells[1]] = hand
    return checks


#every node of the chain with its status, from chain.py
def chain_nodes():
    path = os.path.join(ledger_folder, 'nodes.tmp.csv')
    command = [sys.executable, os.path.join(repo_folder, 'tools', 'chain.py'), '--csv', path]
    subprocess.run(command, cwd=repo_folder, check=True, capture_output=True)
    with open(path) as fp:
        nodes = {row['id']: row for row in csv.DictReader(fp)}
    os.remove(path)
    return nodes


#the paper's lint problems: (line, text)
def lint_problems():
    command = [sys.executable, os.path.join(repo_folder, 'tools', 'lint_paper.py')]
    result = subprocess.run(command, cwd=repo_folder, capture_output=True, text=True)
    return [(int(m.group(1)), m.group(2)) for m in re.finditer(r'line (\d+): (.*)', result.stdout)]


#the findings of the latest outside review of the paper, if any
def paper_findings():
    paths = sorted(glob.glob(os.path.join(repo_folder, 'reviews', '*-paper-verdicts.json')))
    if (len(paths) == 0):
        return '', []
    with open(paths[-1]) as fp:
        verdicts = json.load(fp)
    return os.path.basename(paths[-1])[:10], verdicts['findings']


#a short id for an issue that stays the same while its text does, so the author's Use or Skip sticks
def issue_id(text):
    return hashlib.sha256(text.encode()).hexdigest()[:8]


#which card a node's anchor in the paper lies in: by its \label, or by the phrase of its paragraph
def card_of(anchor, paper_cards):
    if (not anchor.startswith('paper/epi_int_3d.tex')):
        return None
    if ('#' in anchor):
        label = anchor.split('#')[1]
        for card in paper_cards:
            if (label in card['labels']):
                return card['id']
    if ('@' in anchor):
        phrase = ' '.join(anchor.split('@')[1].split())
        for card in paper_cards:
            if (phrase in card['tex']):
                return card['id']
    return None


#the cards, with their issues and facts
def write_cards():
    assets = {}
    if os.path.exists(assets_file):
        with open(assets_file) as fp:
            assets = json.load(fp)
    paper_cards = cards(assets)
    nodes = chain_nodes()
    checks = hand_checks()
    problems = lint_problems()
    review_date, findings = paper_findings()
    by_id = {card['id']: card for card in paper_cards}
    for card in paper_cards:
        card['issues'] = []
        card['items'] = []
        card['equations'] = []
    #the chain's paper items, and the equations of notes/equations.md behind each
    for node in nodes.values():
        card_id = card_of(node['anchor'], paper_cards)
        if (card_id is None):
            continue
        card = by_id[card_id]
        card['items'] += [node['id']]
        for need in node['needs'].split():
            if need.startswith('sympy:') and (need.split(':')[1] in checks):
                card['equations'] += [need.split(':')[1]]
        if (node['status'] != 'ok'):
            text = node['id'] + ' is ' + node['status'] + ': ' + node['reason'] + '.'
            card['issues'] += [{'id': issue_id(text), 'kind': 'chain', 'text': text, 'rec': True,
                'detail': 'Recheck it against what changed, then confirm it (tools/chain.py).'}]
    #lint problems, by line
    for line, text in problems:
        for card in paper_cards:
            if (card['lines'][0] <= line < card['lines'][1]):
                card['issues'] += [{'id': issue_id(text), 'kind': 'lint', 'text': text + ' (line ' + str(line) + ')',
                    'rec': True, 'detail': ''}]
    #the latest paper review's findings, by the section or appendix they name
    for finding in findings:
        where = finding.get('where', '')
        found = re.search(r'(?:Section|Sec\.|§)\s*(\d(?:\.\d)?)|Appendix\s+([A-D])', where)
        card_id = None
        if found and found.group(1):
            card_id = 's' + found.group(1).replace('.', '-')
        elif found:
            card_id = 'app' + found.group(2)
        if (card_id not in by_id):
            card_id = 'abstract'
        rec = finding['verdict'] in ['agree', 'partly']
        ruling = finding.get('tiebreak', {}).get('ruling', '')
        detail = 'Claude ' + finding['verdict'] + 's: ' + finding.get('comment', '')
        if ruling:
            detail += ' Tiebreak (' + finding['tiebreak'].get('judge', '') + '): ' + ruling + '.'
        text = finding['reviewer'].upper() + ' (' + review_date + '): ' + finding['finding']
        by_id[card_id]['issues'] += [{'id': issue_id(text), 'kind': 'review', 'text': text, 'rec': rec,
            'detail': detail, 'fix': finding.get('fix', '')}]
    #facts, and write each card
    folder = os.path.join(db_folder, 'cards')
    os.makedirs(folder, exist_ok=True)
    for card in paper_cards:
        equations = sorted(set(card['equations']))
        card['facts'] = {'items': len(card['items']),
            'items_ok': sum(nodes[i]['status'] == 'ok' for i in card['items']),
            'equations': len(equations), 'equations_hand': sum(checks.get(e) == 'yes' for e in equations)}
        card['state'] = 'attention' if card['issues'] else 'ok'
        document = {key: card[key] for key in ['id', 'order', 'number', 'title', 'section', 'level', 'appendix',
            'html', 'words', 'issues', 'facts', 'state']}
        with open(os.path.join(folder, card['id'] + '.json'), 'w') as fp:
            json.dump(document, fp, indent=1)
    #figures that need uploading: new, or changed since their upload
    for card in paper_cards:
        for png in re.findall(r'data-png="([^"]*)"', card['html']):
            digest = hashlib.sha256(open(os.path.join(figures_folder, png), 'rb').read()).hexdigest()[:12]
            if (assets.get(png, {}).get('hash') != digest):
                print('upload figure', png, 'hash', digest)
    print('cards =', len(paper_cards), ' needing attention =', sum(c['state'] == 'attention' for c in paper_cards))
    return paper_cards, nodes, checks


#the phase of a date, from notes/phases.csv
def phase_of(date, phases):
    phase = phases[0][1]
    for start, name in phases:
        if (start[:10] <= date):
            phase = name
    return phase


#what the project cost, per phase: outside reviews, Claude Code at API prices, the author's hours; and
#what the author paid
def costs():
    with open(os.path.join(repo_folder, 'notes', 'phases.csv')) as fp:
        phases = [(row['start'], row['phase']) for row in csv.DictReader(fp)]
    rows = {}
    def row(phase):
        if (phase not in rows):
            rows[phase] = {'phase': phase, 'cc_usd': 0.0, 'cc_hours': 0.0, 'joe_hours': 0.0, 'outside_usd': 0.0,
                'open': False}
        return rows[phase]
    command = [sys.executable, os.path.join(repo_folder, 'tools', 'session_time.py'), '--csv',
        os.path.join(ledger_folder, 'phases.tmp.csv')]
    subprocess.run(command, cwd=repo_folder, check=True, capture_output=True)
    with open(os.path.join(ledger_folder, 'phases.tmp.csv')) as fp:
        for line in csv.DictReader(fp):
            r = row(line['phase'])
            r['cc_usd'] = round(float(line['cost_usd']), 2)
            r['cc_hours'] = round(float(line['cc_hours']), 2)
            r['joe_hours'] = round(float(line['human_hours']) + float(line['off_hours']), 2)
            r['open'] = (line['open'] == '1')
    os.remove(os.path.join(ledger_folder, 'phases.tmp.csv'))
    reviews = []
    for line in open(os.path.join(repo_folder, 'reviews', 'costs.jsonl')):
        cost = json.loads(line)
        r = row(phase_of(cost['date'], phases))
        r['outside_usd'] = round(r['outside_usd'] + cost['usd'], 4)
        reviews += [{'date': cost['date'], 'activity': cost['activity'], 'reviewer': cost['reviewer'],
            'usd': round(cost['usd'], 4)}]
    #the subscription, one payment per month begun since the project started
    start = datetime.date.fromisoformat(phases[0][0][:10])
    today = datetime.date.today()
    months = (today.year - start.year)*12 + (today.month - start.month)
    if (today.day >= start.day):
        months = months + 1
    subscription = months*subscription_usd_per_month
    phase_rows = list(rows.values())
    outside = sum(r['outside_usd'] for r in phase_rows)
    return {'phases': phase_rows, 'reviews': reviews,
        'subscription': {'usd': subscription, 'months': months, 'per_month': subscription_usd_per_month,
            'since': start.isoformat()},
        'paid_usd': round(subscription + outside, 2), 'outside_usd': round(outside, 2),
        'cc_usd': round(sum(r['cc_usd'] for r in phase_rows), 2),
        'cc_hours': round(sum(r['cc_hours'] for r in phase_rows), 1),
        'joe_hours': round(sum(r['joe_hours'] for r in phase_rows), 1),
        'open_session': any(r['open'] for r in phase_rows)}


#the summary at the top of the page, and the costs
def write_meta(paper_cards, nodes, checks):
    cost = costs()
    commit = subprocess.run(['git', 'log', '-1', '--format=%h %cs'], cwd=repo_folder, capture_output=True,
        text=True).stdout.split()
    hand_needed = [ID for ID, hand in checks.items() if (hand in ['yes', 'no'])]
    summary = {'asOf': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds'),
        'commit': commit[0] if commit else '', 'commitDate': commit[1] if len(commit) > 1 else '',
        'cards': len(paper_cards), 'attention': sum(c['state'] == 'attention' for c in paper_cards),
        'issues': sum(len(c['issues']) for c in paper_cards),
        'chain_ok': sum(n['status'] == 'ok' for n in nodes.values()), 'chain_total': len(nodes),
        'hand_yes': sum(checks[ID] == 'yes' for ID in hand_needed), 'hand_needed': len(hand_needed),
        'words': sum(c['words'] for c in paper_cards),
        'paid_usd': cost['paid_usd'], 'cc_usd': cost['cc_usd'], 'outside_usd': cost['outside_usd'],
        'joe_hours': cost['joe_hours'], 'cc_hours': cost['cc_hours']}
    folder = os.path.join(db_folder, 'meta')
    os.makedirs(folder, exist_ok=True)
    for name, document in [('summary', summary), ('costs', cost)]:
        with open(os.path.join(folder, name + '.json'), 'w') as fp:
            json.dump(document, fp, indent=1)
    print('paid = $' + str(summary['paid_usd']), ' chain', summary['chain_ok'], 'of', summary['chain_total'])


#write the database documents
if os.path.exists(db_folder):
    shutil.rmtree(db_folder)
paper_cards, nodes, checks = write_cards()
write_meta(paper_cards, nodes, checks)
print('wrote', os.path.relpath(db_folder, repo_folder))
