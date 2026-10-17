#!/usr/bin/env python3
"""Time, token and cost accounting per project phase, from Claude Code session logs.

Reads the session logs in ~/.claude/projects/<project>/*.jsonl and reports, per phase:
  - CC time: Claude Code's working time per turn, from the logs' turn_duration records.
  - human time in session: from the end of each turn to the next human prompt; gaps
    longer than --away minutes count as away and are excluded.
  - human time off session: hours logged by hand in notes/human_time.md.
  - tokens: API usage summed over each turn's assistant messages.
  - cost: each session's totalCostUSD (written when a session ends), split across phases in
    proportion to CC time. Sessions still running have no cost yet.

Phases are assigned by prompt time, using notes/phases.csv (columns: start, phase), where
each row's phase applies from its start time until the next row's.

Run:  uv run python tools/session_time.py [--away 30] [--by-session]
"""
import argparse
import csv
import json
import re
from collections import defaultdict
from datetime import datetime
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent


def parse_time(s):
    return datetime.fromisoformat(s.replace('Z', '+00:00'))


def default_log_dir():
    # Claude Code names a project's log folder after its path, with '/' and '.' replaced by '-'
    name = re.sub(r'[/.]', '-', str(REPO))
    return Path.home()/'.claude'/'projects'/name


def load_phases(path):
    with open(path) as f:
        rows = [(parse_time(r['start']), r['phase']) for r in csv.DictReader(f)]
    return sorted(rows)


def phase_at(phases, t):
    current = phases[0][1]
    for start, phase in phases:
        if t >= start:
            current = phase
    return current


def is_human_prompt(d):
    # Typed prompts and accepted suggestions start a turn; prompts queued mid-turn do not.
    return (d.get('type') == 'user' and d.get('turnOrigin') == 'human'
            and d.get('promptSource') != 'queued')


def read_session(path):
    """Return the session's turns and its total cost (None if the session is still open)."""
    records = []
    for line in open(path):
        try:
            records.append(json.loads(line))
        except json.JSONDecodeError:
            continue
    cost = None
    turns = []
    turn = None
    seen_messages = set()
    for d in records:
        if d.get('type') == 'cost-state':
            cost = float(d['totalCostUSD'])
        if 'timestamp' not in d:
            continue
        t = parse_time(d['timestamp'])
        if is_human_prompt(d):
            turn = {'start': t, 'end': None, 'cc_s': 0.0,
                    'tokens': defaultdict(int)}
            turns.append(turn)
        elif turn is None:
            continue
        elif d.get('type') == 'assistant':
            msg = d.get('message', {})
            usage = msg.get('usage')
            if usage and msg.get('id') not in seen_messages:
                seen_messages.add(msg.get('id'))
                for k in ('input_tokens', 'output_tokens', 'cache_read_input_tokens',
                          'cache_creation_input_tokens'):
                    turn['tokens'][k] += usage.get(k) or 0
        elif d.get('type') == 'system' and d.get('subtype') == 'turn_duration':
            turn['end'] = t
            turn['cc_s'] += d['durationMs']/1000
    return turns, cost


def read_off_session(path):
    """Hours from the markdown table in notes/human_time.md: | date | phase | task | hours |."""
    hours = defaultdict(float)
    if not path.exists():
        return hours
    for line in open(path):
        cells = [c.strip() for c in line.strip().strip('|').split('|')]
        if len(cells) == 4 and re.fullmatch(r'\d{4}-\d{2}-\d{2}', cells[0]):
            hours[cells[1]] += float(cells[3])
    return hours


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--logs', type=Path, default=default_log_dir())
    p.add_argument('--phases', type=Path, default=REPO/'notes'/'phases.csv')
    p.add_argument('--human-log', type=Path, default=REPO/'notes'/'human_time.md')
    p.add_argument('--away', type=float, default=30, help='minutes; longer gaps are away')
    p.add_argument('--by-session', action='store_true', help='also print each session')
    args = p.parse_args()

    phases = load_phases(args.phases)
    order = [ph for _, ph in phases]
    total = defaultdict(lambda: defaultdict(float))
    open_sessions = 0
    for path in sorted(args.logs.glob('*.jsonl'), key=lambda f: f.stat().st_mtime):
        turns, cost = read_session(path)
        if not turns:
            continue
        cc_by_phase = defaultdict(float)
        for k, turn in enumerate(turns):
            ph = phase_at(phases, turn['start'])
            row = total[ph]
            row['turns'] += 1
            row['cc_h'] += turn['cc_s']/3600
            cc_by_phase[ph] += turn['cc_s']
            for key, n in turn['tokens'].items():
                row[key] += n
            # human time: from this turn's end to the next prompt, unless away
            if turn['end'] and k + 1 < len(turns):
                gap = (turns[k + 1]['start'] - turn['end']).total_seconds()
                if 0 < gap <= args.away*60:
                    row['human_h'] += gap/3600
                else:
                    row['away_gaps'] += 1
        cc_total = sum(cc_by_phase.values())
        if cost is None:
            open_sessions += 1
        elif cc_total > 0:
            for ph, s in cc_by_phase.items():
                total[ph]['cost'] += cost*s/cc_total
        if args.by_session:
            print(f'{path.stem[:8]}  {turns[0]["start"]:%Y-%m-%d}  turns {len(turns):3d}  '
                  f'CC {cc_total/3600:5.2f} h  cost ' + (f'${cost:.2f}' if cost is not None else 'open'))

    off = read_off_session(args.human_log)
    print(f'\n{"phase":<10} {"turns":>5} {"CC h":>6} {"human h":>8} {"off h":>6} '
          f'{"in tok":>8} {"out tok":>8} {"cache rd":>10} {"cache wr":>9} {"cost $":>7}')
    for ph in order + sorted(set(off) - set(order)):
        r = total.get(ph, {})
        if not r and ph not in off:
            continue
        print(f'{ph:<10} {int(r.get("turns", 0)):5d} {r.get("cc_h", 0):6.2f} '
              f'{r.get("human_h", 0):8.2f} {off.get(ph, 0):6.2f} '
              f'{int(r.get("input_tokens", 0)):8d} {int(r.get("output_tokens", 0)):8d} '
              f'{int(r.get("cache_read_input_tokens", 0)):10d} '
              f'{int(r.get("cache_creation_input_tokens", 0)):9d} {r.get("cost", 0):7.2f}')
    if open_sessions:
        print(f'\n{open_sessions} session(s) still open: their cost is not yet in the table.')


if __name__ == '__main__':
    main()
