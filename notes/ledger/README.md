# The paper ledger

A private claude.ai page that shows the paper one section or subsection at a time, each card with its
open issues and a box for the author's notes to Claude: https://claude.ai/artifact/KpgQvmghBVn7c4nTmFWxLZ

- `ledger.html`: the page. Its content lives in the page's database, not in the file:
  - `cards/<id>`: one per section or subsection (`abstract`, `s1`, `s2-1`, …, `appA`, …). Each holds the
    section's text as HTML, with the math left as LaTeX for MathJax, the equations numbered as in the PDF,
    and the number macros filled in; its issues; and its facts (chain items current, equations hand-checked).
    Each issue has an id that stays the same while its text does, `kind` (chain, lint or review), `rec`
    (Claude recommends acting on it) and a `detail`.
  - `meta/summary` and `meta/costs`: the top bar's totals and the Costs view.
  - `notes/<id>`: the author's note on a card (`text`) and his Use/Skip choice on each issue (`fixes`).
  - `replies/<id>`: Claude's reply to that note (`reply`, `at`, `status`: done or question).
- `assets.json`: the ids of the figures uploaded to the page, with each PNG's hash; `tools/ledger.py`
  says when a figure has changed and needs uploading again.
- `notes.json`: the record of the author's notes and Claude's replies, exported each time Claude acts on
  them, so the editing record is in git.
- `db/`: the documents `tools/ledger.py` writes for upload (git-ignored).

To refresh the page after a change to the paper, the chain or the reviews: `uv run python tools/ledger.py`,
then Claude uploads `db/` with one ArtifactData batch (cards, meta) and any figure it names.

To act on notes: the author types them on the page and clicks "Apply changes" (it sends a comment to the
Claude Code session watching the page), or tells Claude Code "apply my ledger notes". Claude then
1. reads `notes/*` and `replies/*`, and takes each note newer than its reply;
2. edits the paper (the paper-revision rule of CLAUDE.md: lint, style guide, render the changed pages),
   and acts on each issue marked Use;
3. runs the knock-on pass and `tools/ledger.py`, uploads the cards, writes a reply on each card it acted on
   (what changed, or a question back), and exports the notes and replies to `notes.json`.
