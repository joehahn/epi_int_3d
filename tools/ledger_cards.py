#ledger_cards.py
#
#by Joe Hahn and Claude Code, joe.hahn@jmh-datasciences.com, 10 November 2026.
#
#split the paper into the ledger's cards, one per section or subsection, each holding its text as HTML
#for the ledger page (tools/ledger.py). Math stays LaTeX, which the page typesets with MathJax; the
#equations get the numbers the PDF gives them (checked against the compiled paper's labels), the number
#macros of paper/numbers/ are filled in with their values, and citations and references are written out.
#
#used by tools/ledger.py; on its own, from the repo's top folder,
#   uv run python tools/ledger_cards.py      prints each card's number, title and word count

import glob
import html
import json
import os
import re
import subprocess
import tempfile

#paths
repo_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
paper_folder = os.path.join(repo_folder, 'paper')
paper_file = os.path.join(paper_folder, 'epi_int_3d.tex')
citations_file = os.path.join(repo_folder, 'notes', 'ledger', 'citations.json')
chain_file = os.path.join(repo_folder, 'notes', 'chain.json')
#the executed notebooks, exported to HTML on the project's GitHub Pages site by tools/export_notebooks.py
notebook_site = 'https://joehahn.github.io/epi_int_3d/notebooks/'

#a place in a cited paper: "Eqn.\ (A21) of \cite{S84}", "Fig.~5 of \cite{H07}", "Appendix A of \citealt{H25}",
#"Eqns.\ (88) and (91)--(92) of \cite{S84}" (the first number), and \citep[Eqns.\ (97)--(98) of][]{S84}
_place = (r'(?P<kind>Eqns?\.|Figs?\.|Sections?|Appendix)(?:\\ |~|\s)+\(?(?P<number>[A-Z]?\d+(?:\.\d+)*|[A-Z])\)?'
    r'(?:--\(?[A-Z]?\d+\)?)?(?:\s+and\s+\(?[A-Z]?\d+\)?(?:--\(?[A-Z]?\d+\)?)?)?')
LOCATOR = re.compile(_place + r'\s+of\s+\\(?:citealt|citet|cite)\{(?P<key>\w+)\}')
LOCATOR_CITEP = re.compile(r'\\citep\[' + _place + r'\s+of\]\[\]\{(?P<key>\w+)\}')


#the text of a {...} group starting at text[start] == '{', and the position after it
def braced(text, start):
    depth = 0
    for position in range(start, len(text)):
        if (text[position] == '{') and ((position == 0) or (text[position - 1] != '\\')):
            depth += 1
        elif (text[position] == '}') and (text[position - 1] != '\\'):
            depth -= 1
            if (depth == 0):
                return text[start + 1:position], position + 1
    raise ValueError('unbalanced braces after: ' + text[start:start + 60])


#replace each \name{arg} (with optional [..] arguments before it) by what replace(arg, options) returns
def replace_command(text, name, replace):
    pattern = re.compile(r'\\' + name + r'((?:\[[^\]]*\])*)\{')
    while True:
        m = pattern.search(text)
        if (m is None):
            return text
        arg, end = braced(text, m.end() - 1)
        options = re.findall(r'\[([^\]]*)\]', m.group(1))
        text = text[:m.start()] + replace(arg, options) + text[end:]


#the compiled paper's labels (number by key) and bibliography (author and year by key)
def compiled_paper():
    with tempfile.TemporaryDirectory() as folder:
        subprocess.run(['tectonic', '-k', '-o', folder, paper_file], cwd=paper_folder, check=True,
            capture_output=True)
        aux = open(os.path.join(folder, 'epi_int_3d.aux')).read()
        bbl = open(os.path.join(folder, 'epi_int_3d.bbl')).read()
    #a label's number may hold one level of braces: {\hbox {A1}a} is A1a, a lettered line in an appendix
    labels = {key: re.sub(r'\\hbox \{([^}]*)\}', r'\1', value)
        for key, value in re.findall(r'\\newlabel\{([^}]*)\}\{\{((?:[^{}]|\{[^{}]*\})*)\}', aux)}
    references = {}
    for label, key in re.findall(r'\\bibitem\[\{(.*?)\}\]\{([^}]*)\}', bbl, flags=re.S):
        m = re.match(r'.*\((\d{4}[a-z]?)\)(.*)$', ' '.join(label.split()))
        authors = m.group(2).replace('{', '').replace('}', '').replace('\\&', '&').replace('~', ' ')
        references[key] = (authors, m.group(1))
    return labels, references


#the number macros of paper/numbers/*.tex: value by name
def number_macros():
    macros = {}
    for path in sorted(glob.glob(os.path.join(paper_folder, 'numbers', '*.tex'))):
        text = open(path).read()
        for m in re.finditer(r'\\newcommand\{\\(\w+)\}\{', text):
            macros[m.group(1)], end = braced(text, m.end() - 1)
    return macros


#remove TeX comments (an unescaped % to the end of the line)
def strip_comments(text):
    return '\n'.join(re.sub(r'(?<!\\)%.*', '', line) for line in text.split('\n'))


#fill in the number macros
def fill_numbers(text, macros):
    return re.sub(r'\\(\w+)(?![a-zA-Z])', lambda m: macros.get(m.group(1), m.group(0)), text)


#turns the LaTeX of a card into HTML; keeps the equation and figure counters across cards
#an eqnarray's rows: split at its own \\, not at those inside a nested array or cases
def split_rows(body):
    rows, depth, start = [], 0, 0
    for m in re.finditer(r'\\begin\{[^}]*\}|\\end\{[^}]*\}|\\\\', body):
        token = m.group(0)
        if token.startswith('\\begin'):
            depth += 1
        elif token.startswith('\\end'):
            depth -= 1
        elif (depth == 0):
            rows += [body[start:m.start()]]
            start = m.end()
    return rows + [body[start:]]


class Converter:

    def __init__(self, labels, references, macros):
        self.labels = labels
        self.references = references
        self.macros = macros
        self.equation = 0
        self.figure = 0
        self.prefix = ''                 #'' in the body, the appendix letter in an appendix
        self.mismatches = []             #labels whose number here differs from the PDF's
        self.leftovers = set()           #LaTeX commands this converter doesn't know
        self.depth = 0
        self.citations = json.load(open(citations_file))['papers'] if os.path.exists(citations_file) else {}
        self.missing_locators = set()     #cited places that citations.json has no page for
        self.card_number, self.subsub = '', 0   #the card being converted, and its subsubsections so far
        self.pieces = []
        self.sub = None                  #in a subequations environment, the number of lines so far


    #the next equation number, checked against the PDF's when the equation has a label
    def number_equation(self, label):
        if (self.sub is None):
            self.equation += 1
            number = self.prefix + str(self.equation)
        else:                            #6a, 6b, ... on one number
            number = self.prefix + str(self.equation) + 'abcdefghijklmnopqrstuvwxyz'[self.sub]
            self.sub += 1
        if ((label is not None) and (self.labels.get(label) != number)):
            self.mismatches += [(label, number, self.labels.get(label))]
        return number

    #a displayed equation or eqnarray, as LaTeX for MathJax with \tag{n} for its numbers
    def display(self, environment, body):
        body = fill_numbers(body, self.macros)
        if (environment == 'eqnarray'):
            rows = split_rows(body)
            out = []
            for row in rows:
                label = re.search(r'\\label\{([^}]*)\}', row)
                row = re.sub(r'\\label\{[^}]*\}', '', row)
                if ('\\nonumber' in row):
                    row = row.replace('\\nonumber', '\\notag')
                elif (row.strip() != ''):
                    row = row + r' \tag{' + self.number_equation(label and label.group(1)) + '}'
                out += [row.strip()]
            latex = '\\begin{eqnarray}' + ' \\\\\n'.join(out) + '\\end{eqnarray}'
        else:
            label = re.search(r'\\label\{([^}]*)\}', body)
            body = re.sub(r'\\label\{[^}]*\}', '', body).strip()
            if ((environment == 'equation*') or ('\\nonumber' in body)):
                body = body.replace('\\nonumber', '')
                latex = '\\begin{equation*}' + body + '\\end{equation*}'
            else:
                number = self.number_equation(label and label.group(1))
                latex = '\\begin{equation}' + body + r' \tag{' + number + '}\\end{equation}'
        return '<span class="eq">' + html.escape(latex, quote=False) + '</span>'

    #a subsubsection heading, numbered as the paper numbers it (3.1.2)
    def subsubsection(self, title):
        self.subsub += 1
        return ('<h4 class="sub">' + self.card_number + '.' + str(self.subsub) + '\u2002'
            + html.escape(title, quote=False) + '</h4>')

    #a locator and its citation wrapped in a link to that page of the cited paper's PDF in the ledger's file store,
    #or left as it is (and noted in self.missing_locators) when citations.json lacks the paper or the page
    def locator_link(self, m, keep):
        kind = {'eqn': 'eqn', 'fig': 'fig', 'sec': 'section', 'app': 'appendix'}[m.group('kind')[:3].lower()]
        key = kind + ':' + m.group('number')
        paper = self.citations.get(m.group('key'), {})
        page = paper.get('pages', {}).get(key)
        if (page is None) or (not paper.get('url')):
            self.missing_locators.add(m.group('key') + ' ' + key)
            return m.group(0)
        #the whole paper, in a new tab, at that page
        if paper.get('nopage'):
            href, title = paper['url'], 'Opens the paper at its publisher (no free copy); the place is on page ' + str(page)
        else:
            href, title = paper['url'] + '#page=' + str(page), 'Opens the paper at page ' + str(page)
        if key in paper.get('notes', {}):
            title += ', where it is ' + paper['notes'][key]
        return (keep('<a class="ref" href="' + href + '" target="_blank" rel="noopener" title="' + title + '">')
            + m.group(0) + keep('</a>'))

    #citations as the paper prints them
    def cite(self, keys, options, style):
        names = []
        for key in keys.split(','):
            authors, year = self.references[key.strip()]
            if (style == 'citet'):
                names += [authors + ' (' + year + ')']
            else:
                names += [authors + ' ' + year]
        text = '; '.join(names)
        if (style == 'citep'):
            options = [o.replace('\\ ', ' ').replace('---', '\u2014').replace('--', '\u2013').replace('~', '\u00a0')
                for o in options]                #the TeX typography of the running text
            before, after = '', ''
            if (len(options) == 1):
                after = ', ' + options[0]
            elif (len(options) == 2):
                before = options[0] + ' ' if options[0] else ''
                after = ', ' + options[1] if options[1] else ''
            text = '(' + before + text + after + ')'
        return html.escape(text, quote=False)

    #a figure: its image (the PNG beside the PDF) and caption
    def figure_html(self, body, assets):
        self.figure += 1
        name = re.search(r'\\includegraphics(?:\[[^\]]*\])?\{([^}]*)\}', body).group(1)
        png = os.path.basename(name).replace('.pdf', '.png')
        caption, end = braced(body, body.index('\\caption{') + len('\\caption'))
        label = re.search(r'\\label\{([^}]*)\}', body)
        if (label and (self.labels.get(label.group(1)) != str(self.figure))):
            self.mismatches += [(label.group(1), str(self.figure), self.labels.get(label.group(1)))]
        asset = assets.get(png, {}).get('id', '')
        return ('<figure class="fig" data-png="' + png + '" data-asset="' + asset + '">'
            + '<figcaption><b>Figure ' + str(self.figure) + '.</b> ' + self.inline(caption, {}) + '</figcaption></figure>')

    #running text: paragraphs, inline math, citations, references, footnotes, lists
    def inline(self, text, assets):
        top = (self.depth == 0)                          #nested calls (footnotes, captions) share the pieces
        if top:
            self.pieces = []
        self.depth += 1
        pieces = self.pieces                             #protected HTML, put back at the end
        def keep(fragment):
            pieces.append(fragment)
            return '\x00' + str(len(pieces) - 1) + '\x00'
        #figures and displayed math first
        text = re.sub(r'\\begin\{(figure\*?)\}(?:\[[^\]]*\])?(.*?)\\end\{\1\}',
            lambda m: keep(self.figure_html(m.group(2), assets)), text, flags=re.S)
        def subequations(body):            #numbered 6a, 6b, ... under one number
            self.equation += 1
            label = re.match(r'\s*\\label\{([^}]*)\}', body)
            number = self.prefix + str(self.equation)
            if (label and (self.labels.get(label.group(1)) != number)):
                self.mismatches += [(label.group(1), number, self.labels.get(label.group(1)))]
            self.sub = 0
            out = re.sub(r'\\begin\{(equation\*?|eqnarray)\}(.*?)\\end\{\1\}',
                lambda e: keep(self.display(e.group(1), e.group(2))), body, flags=re.S)
            self.sub = None
            return re.sub(r'^\s*\\label\{[^}]*\}', '', out)
        def displayed(m):                  #in the paper's order, so the numbers count up
            if (m.group(1) == 'subequations'):
                return subequations(m.group(2))
            return keep(self.display(m.group(1), m.group(2)))
        text = re.sub(r'\\begin\{(subequations|equation\*?|eqnarray)\}(.*?)\\end\{\1\}', displayed, text, flags=re.S)
        #inline math
        text = re.sub(r'(?<!\\)\$(.+?)(?<!\\)\$',
            lambda m: keep('$' + html.escape(fill_numbers(m.group(1), self.macros), quote=False) + '$'), text, flags=re.S)
        #the number macros in running text (a {} after one keeps the space after it in TeX)
        text = re.sub(r'(\\(?:vres|bw)\w+)\{\}', r'\1', text)
        text = fill_numbers(text, self.macros)
        #a place in a cited paper, "Eqn. (A21) of \\cite{S84}", links to that page of its PDF (notes/ledger/citations.json)
        text = re.sub(LOCATOR, lambda m: self.locator_link(m, keep), text)
        text = re.sub(LOCATOR_CITEP, lambda m: self.locator_link(m, keep), text)
        #commands with arguments
        text = replace_command(text, 'footnote', lambda arg, o: keep('<span class="fn">'
            + self.inline(arg, assets) + '</span>'))
        text = re.sub(r'\\nblink\{([^}]*)\}\{([^}]*)\}', lambda m: keep(' <a class="nb" href="' + html.escape(m.group(1))
            + '" target="_blank" rel="noopener">notebook ' + html.escape(m.group(2)) + '</a>'), text)
        text = re.sub(r'\\href\{([^}]*)\}\{((?:[^{}]|\{[^{}]*\})*)\}', lambda m: keep('<a href="' + html.escape(m.group(1))
            + '" target="_blank" rel="noopener">' + self.inline(m.group(2), assets) + '</a>'), text)
        text = replace_command(text, 'url', lambda arg, o: keep('<a href="' + html.escape(arg) + '" target="_blank" rel="noopener">'
            + html.escape(arg) + '</a>'))
        for style in ['citep', 'citet', 'citealt', 'cite']:
            text = replace_command(text, style, lambda arg, o, style=style:
                keep(self.cite(arg, o, 'citet' if style == 'cite' else style)))
        text = replace_command(text, 'ref', lambda arg, o: keep(html.escape(self.labels.get(arg, '??'))))
        text = replace_command(text, 'label', lambda arg, o: '')
        text = replace_command(text, 'paragraph', lambda arg, o: keep('<b>' + html.escape(arg) + '.</b>'))
        for name in ['code', 'texttt']:
            text = replace_command(text, name, lambda arg, o: keep('<code>'
                + html.escape(arg.replace('\\_', '_'), quote=False) + '</code>'))
        text = replace_command(text, 'emph', lambda arg, o: keep('<em>' + self.inline(arg, assets) + '</em>'))
        text = re.sub(r'\\(centering|noindent|appendix)\b', '', text)
        #subsubsection headings (sections and subsections are the cards themselves)
        text = re.sub(r'\\subsubsection\*?\{([^}]*)\}(\s*\\label\{[^}]*\})?',
            lambda m: keep(self.subsubsection(m.group(1))), text)
        #lists
        text = text.replace('\\begin{enumerate}', keep('<ol>')).replace('\\end{enumerate}', keep('</ol>'))
        text = re.sub(r'\\item\s*', lambda m: keep('<li>'), text)
        #the remaining text, escaped, then TeX's typography
        text = html.escape(text, quote=False)
        for tex, out in [('---', '\u2014'), ('--', '\u2013'), ('``', '\u201c'), ("''", '\u201d'), ('~', '\u00a0'),
                ('\\%', '%'), ('\\_', '_'), ('\\&amp;', '&amp;'), ('\\#', '#'), ('\\$', '$'), ('\\S', '\u00a7'),
                ('\\ldots', '\u2026'), ('\\,', '\u2009'), ('\\ ', ' '), ('.\\@', '.'), ('\\@', '')]:
            text = text.replace(tex, out)
        text = re.sub(r'(?<=\w)\\/', '', text)
        text = re.sub(r'\\raggedright\s*|\\setcounter\{[^}]*\}\{[^}]*\}\s*', '', text)      #layout only
        self.leftovers |= set(re.findall(r'\\[a-zA-Z]+', text))
        #paragraphs
        paragraphs = [p.strip() for p in re.split(r'\n\s*\n', text) if p.strip()]
        if (len(paragraphs) > 1):
            text = ''.join('<p>' + ' '.join(p.split()) + '</p>' for p in paragraphs)
        else:
            text = ' '.join(text.split())
        #put back the protected pieces, innermost last
        self.depth -= 1
        while top and ('\x00' in text):
            text = re.sub(r'\x00(\d+)\x00', lambda m: pieces[int(m.group(1))], text)
        text = re.sub(r'<p>(<h4 class="sub">.*?</h4>)</p>', r'\1', text)
        #a figure is a block of its own: inside <p> the browser splits it into two empty paragraphs
        text = re.sub(r'<p>(<figure\b.*?</figure>)</p>', r'\1', text, flags=re.S)
        return text


#the notebook sections that back each claim or figure: chain.json nodes with a 'notebook' field, e.g.
#["bending_waves#test2"], put as a link at the end of the claim's paragraph (anchor 'paper/...tex@phrase') or of the
#figure's caption (anchor 'paper/...tex#fig:label'), so the author sees them in the text he reads
def notebook_links(body):
    nodes = json.load(open(chain_file))
    nodes = nodes['nodes'] if isinstance(nodes, dict) else nodes
    def marker(sections):
        return ''.join('\\nblink{' + notebook_site + sec.split('#')[0] + '.html' + ('#' + sec.split('#')[1] if '#' in sec else '')
            + '}{' + (('\u00a7' + re.sub(r'\D', '', sec.split('#')[1])) if '#' in sec else sec) + '}' for sec in sections)
    for node in nodes:
        sections = node.get('notebook')
        anchor = node.get('anchor', '')
        if (not sections) or (not anchor.startswith('paper/')):
            continue
        if ('@' in anchor):
            phrase = ' '.join(anchor.split('@', 1)[1].split())
            pieces = re.split(r'(\n\s*\n)', body)
            for k, piece in enumerate(pieces):
                if (phrase in ' '.join(piece.split())):
                    pieces[k] = piece.rstrip() + marker(sections)
                    break
            body = ''.join(pieces)
        elif ('#fig:' in anchor):
            at = body.find('\\label{' + anchor.split('#', 1)[1] + '}')
            cap = body.rfind('\\caption{', 0, at)
            if ((at < 0) or (cap < 0)):
                continue
            _, end = braced(body, cap + len('\\caption'))
            body = body[:end - 1] + marker(sections) + body[end - 1:]
    return body


#the paper's cards: [{id, order, number, title, level, html, words, labels, lines, tex}]
def cards(assets=None):
    assets = assets or {}
    labels, references = compiled_paper()
    macros = number_macros()
    source = open(paper_file).read()
    begin = source.index('\\begin{document}')
    end = source.index('\\end{document}')
    first_line = source[:begin].count('\n') + 1
    body = strip_comments(source[begin:end])
    body = re.sub(r'\\bibliographystyle\{[^}]*\}|\\bibliography\{[^}]*\}', '', body)
    body = notebook_links(body)
    #headings
    heads = []
    for m in re.finditer(r'\\(section|subsection)\*?\{|\\appendix|\\begin\{abstract\}|\\end\{abstract\}', body):
        heads += [m]
    converter = Converter(labels, references, macros)
    out = []
    section = 0
    subsection = 0
    appendix = False
    letters = 'ABCDEFGH'
    for k, m in enumerate(heads):
        start = m.end()
        if (m.group(0) == '\\appendix'):
            appendix = True
            section = 0
            converter.equation = 0
            continue
        if (m.group(0) == '\\end{abstract}'):
            continue
        if (m.group(0) == '\\begin{abstract}'):
            number, title, level, card_id = '', 'Abstract', 1, 'abstract'
        else:
            title, start = braced(body, m.end() - 1)
            if (m.group(1) == 'section'):
                section += 1
                subsection = 0
                level = 1
                number = letters[section - 1] if appendix else str(section)
                converter.prefix = number if appendix else ''
                if appendix:                     #each appendix restarts at 1 (\setcounter{equation}{0}): A1, B1, ...
                    converter.equation = 0
            else:
                subsection += 1
                level = 2
                number = (letters[section - 1] if appendix else str(section)) + '.' + str(subsection)
            card_id = ('app' if appendix else 's') + number.replace('.', '-')
        stop = heads[k + 1].start() if (k + 1 < len(heads)) else len(body)
        tex = body[start:stop]
        converter.leftovers = set()
        converter.card_number, converter.subsub = number, 0
        text_html = converter.inline(tex, assets)
        plain = re.sub(r'<[^>]+>', ' ', text_html)
        words = len(re.findall(r'[A-Za-z]{2,}', re.sub(r'\$[^$]*\$|\\begin\{.*?\\end\{[a-z*]+\}', ' ', plain, flags=re.S)))
        if (level == 1):
            section_title = ((number + ' ') if number else '') + re.sub(r'<[^>]+>', '', converter.inline(title, assets))
        if (text_html.strip() == ''):
            continue                         #a section that opens straight into its first subsection
        out += [{'id': card_id, 'section': section_title, 'order': len(out), 'number': number,
            'title': converter.inline(title, assets), 'level': level, 'appendix': appendix, 'html': text_html,
            'words': words, 'labels': re.findall(r'\\label\{([^}]*)\}', tex),
            'lines': [first_line + body[:start].count('\n'), first_line + body[:stop].count('\n')],
            'notebooks': sorted(set(re.findall(r'\\url\{(https://github\.com/[^}]*\.ipynb)\}', tex))),
            'tex': ' '.join(tex.split()), 'leftovers': sorted(converter.leftovers)}]
    for missing in sorted(converter.missing_locators):
        print('citations.json has no page for', missing, '(the ledger shows it unlinked)')
    if converter.mismatches:
        raise ValueError('numbers that differ from the PDF (label, here, PDF): ' + str(converter.mismatches))
    return out


if (__name__ == '__main__'):
    for card in cards():
        print(card['id'].ljust(9), card['number'].ljust(5), str(card['words']).rjust(5), 'words  ',
            card['title'][:50], '  lines', card['lines'], ' leftover:', ' '.join(card['leftovers']))
