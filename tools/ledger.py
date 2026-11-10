#ledger.py
#
#by Joe Hahn and Claude Code, joe.hahn@jmh-datasciences.com, 10 November 2026.
#
#export the project's ledger as three CSV files in notes/ledger/, which the ledger dashboard on claude.ai
#reads (Phase 3b of PLAN.md):
#   nodes.csv      every node of the dependency chain (notes/chain.json), its status, what it needs and
#                  what needs it, with the author's hand check of each equation (notes/equations.md)
#   findings.csv   every finding of the outside reviews, with Claude's verdict and the tiebreak ruling
#                  (reviews/*-verdicts.json)
#   costs.csv      what the project cost: outside reviews (reviews/costs.jsonl), Claude Code per phase
#                  and the author's hours (tools/session_time.py, notes/human_time.md)
#
#usage, from the repo's top folder:
#   uv run python tools/ledger.py
#then Claude uploads the three files to the dashboard.

import csv
import datetime
import glob
import json
import os
import re
import subprocess
import sys

#paths
repo_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ledger_folder = os.path.join(repo_folder, 'notes', 'ledger')
nodes_file = os.path.join(ledger_folder, 'nodes.csv')
findings_file = os.path.join(ledger_folder, 'findings.csv')
costs_file = os.path.join(ledger_folder, 'costs.csv')
phases_table = os.path.join(ledger_folder, 'phases.csv')

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


#the nodes, from chain.py, with the hand checks and each node's file
def write_nodes():
    command = [sys.executable, os.path.join(repo_folder, 'tools', 'chain.py'), '--csv', nodes_file]
    subprocess.run(command, cwd=repo_folder, check=True, capture_output=True)
    with open(nodes_file) as fp:
        rows = list(csv.DictReader(fp))
    checks = hand_checks()
    for row in rows:
        row['hand'] = ''
        if (row['kind'] == 'equation'):
            row['hand'] = checks.get(row['id'].split(':')[1], '')
        row['file'] = row['anchor'].split('::')[0].split('#')[0].split('@')[0]
    columns = list(rows[0].keys())
    with open(nodes_file, 'w', newline='') as fp:
        writer = csv.DictWriter(fp, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    print('nodes =', len(rows))


#the findings of every review, with verdicts and rulings
def write_findings():
    rows = []
    for path in sorted(glob.glob(os.path.join(repo_folder, 'reviews', '*-verdicts.json'))):
        with open(path) as fp:
            verdicts = json.load(fp)
        date = os.path.basename(path)[:10]
        for finding in verdicts['findings']:
            tiebreak = finding.get('tiebreak', {})
            rows += [{'key': date + '-' + verdicts['activity'] + '-' + finding['id'], 'date': date,
                'activity': verdicts['activity'], 'file': ' '.join(verdicts['inputs']), 'id': finding['id'],
                'reviewer': finding['reviewer'], 'severity': finding.get('severity', ''),
                'where': finding.get('where', ''), 'finding': finding['finding'], 'verdict': finding['verdict'],
                'comment': finding.get('comment', ''), 'fix': finding.get('fix', ''),
                'ruling': tiebreak.get('ruling', ''), 'judge': tiebreak.get('judge', ''),
                'ruling_reason': tiebreak.get('reason', ''), 'reply': finding.get('reply', '')}]
    columns = ['key', 'date', 'activity', 'file', 'id', 'reviewer', 'severity', 'where', 'finding', 'verdict',
        'comment', 'fix', 'ruling', 'judge', 'ruling_reason', 'reply']
    with open(findings_file, 'w', newline='') as fp:
        writer = csv.DictWriter(fp, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    print('findings =', len(rows))


#the phase of a date, from notes/phases.csv
def phase_of(date, phases):
    phase = phases[0][1]
    for start, name in phases:
        if (start[:10] <= date):
            phase = name
    return phase


#what the project cost: outside reviews, Claude Code and the author's hours
def write_costs():
    with open(os.path.join(repo_folder, 'notes', 'phases.csv')) as fp:
        phases = [(row['start'], row['phase']) for row in csv.DictReader(fp)]
    rows = []
    #outside reviews, each call
    for line in open(os.path.join(repo_folder, 'reviews', 'costs.jsonl')):
        cost = json.loads(line)
        rows += [{'who': 'outside review', 'phase': phase_of(cost['date'], phases), 'activity': cost['activity'],
            'model': cost['reviewer'], 'usd': cost['usd'], 'hours': '', 'date': cost['date'],
            'note': cost.get('note', '')}]
    #Claude Code and the author, per phase
    command = [sys.executable, os.path.join(repo_folder, 'tools', 'session_time.py'), '--csv', phases_table]
    subprocess.run(command, cwd=repo_folder, check=True, capture_output=True)
    with open(phases_table) as fp:
        for row in csv.DictReader(fp):
            note = ''
            if (row['open'] == '1'):
                note = 'cost of the open session not yet counted'
            rows += [{'who': 'Claude Code', 'phase': row['phase'], 'activity': 'all', 'model': 'claude',
                'usd': row['cost_usd'], 'hours': row['cc_hours'], 'date': '', 'note': note}]
            hours = float(row['human_hours']) + float(row['off_hours'])
            rows += [{'who': 'Joe', 'phase': row['phase'], 'activity': 'all', 'model': '',
                'usd': '', 'hours': round(hours, 3), 'date': '', 'note': 'in session plus off session'}]
    os.remove(phases_table)
    #the subscription, one payment per month begun since the project started
    start = datetime.date.fromisoformat(phases[0][0][:10])
    today = datetime.date.today()
    months = (today.year - start.year)*12 + (today.month - start.month)
    if (today.day >= start.day):
        months = months + 1
    rows += [{'who': 'Joe (paid)', 'phase': 'all', 'activity': 'Claude Code subscription', 'model': 'claude',
        'usd': round(months*subscription_usd_per_month, 2), 'hours': '', 'date': today.isoformat(),
        'note': str(months) + ' month(s) at $' + str(int(subscription_usd_per_month)) + ' from ' + start.isoformat()}]
    columns = ['who', 'phase', 'activity', 'model', 'usd', 'hours', 'date', 'note']
    with open(costs_file, 'w', newline='') as fp:
        writer = csv.DictWriter(fp, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    print('cost rows =', len(rows))


#write the three files
os.makedirs(ledger_folder, exist_ok=True)
write_nodes()
write_findings()
write_costs()
print('wrote', os.path.relpath(ledger_folder, repo_folder))
