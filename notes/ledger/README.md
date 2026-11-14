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
  - `notes/<id>`: the author's note on a card (`text`) and his Use/Skip choice on each issue (`fixes`). A note on one of a
    card's paragraphs, subsubsections or figures has the id `<card>~p<k>`, `<card>~sub<k>` or `<card>~fig-<png name>`, with `card` and `part`
    (e.g. "4.1.2 Massless particles", "Figure 3", or a paragraph's section and opening words, since k shifts as
    paragraphs are added) saying where it applies; its reply has the same id.
  - `replies/<id>`: Claude's reply to that note (`reply`, `at`, `status`: done or question). Set `at` to the real current time (`date -u`), never a guess: the page treats a note older than its reply as handled.
- `citations.json`: each cited paper's public PDF (`url`: ADS's scan, or arXiv) and the page, in that PDF, of each
  place the paper cites ("Eqn. (A21) of Shu (1984)"); the cards link the place to `url#page=N`, which opens the whole
  paper in a new tab at that page. A paper with no free copy (BL94) links to its DOI (`nopage`). When the paper cites a
  new place, `tools/ledger.py` says so; find its page in the public PDF (by eye for scans) and add it. The ledger's
  own file store can't serve these: a link that leaves the page can't reach it ("not found").
- Notebook links: a chain node (`notes/chain.json`) with a `notebook` field, e.g. `["bending_waves#test2"]`, puts a
  "notebook §2" link at the end of its claim's paragraph or its figure's caption on the card, to that section of the
  notebook's HTML on the GitHub Pages site (`docs/notebooks/`, written by `tools/export_notebooks.py`; each section's
  heading cell starts with `<a id="testN"></a>`). Export the notebooks again after re-executing them, and push.
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
