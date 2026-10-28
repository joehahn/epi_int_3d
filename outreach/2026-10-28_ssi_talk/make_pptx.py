#make_pptx.py
#
#builds Using_Claude_Code_for_Research.pptx in the author's own deck style (../slide_style.md):
#4:3, light yellow, Calibri, a header line ending in a colon, nested round bullets, plots beside the text.
#run: uv run --with python-pptx python make_pptx.py

from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.oxml.ns import qn
from lxml import etree

FONT, MONO = 'Calibri', 'Courier New'
BLACK, GREY = RGBColor(0, 0, 0), RGBColor(0x59, 0x59, 0x59)
NAVY, GREEN = RGBColor(0x1F, 0x38, 0x64), RGBColor(0x1E, 0x56, 0x31)
YELLOW, WHITE = RGBColor(0xFF, 0xFE, 0xE5), RGBColor(0xFF, 0xFF, 0xFF)

def textbox(slide, x, y, w, h, framed=False):
    """framed: a white box with a thin grey border, like a screenshot."""
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    for side in ('left', 'right', 'top', 'bottom'):
        setattr(tf, 'margin_' + side, Inches(0.15 if framed else 0.05))
    if framed:
        tb.fill.solid()
        tb.fill.fore_color.rgb = WHITE
        tb.line.color.rgb = GREY
        tb.line.width = Pt(0.75)
    return tf

def runs(p, text, size, font=FONT, color=BLACK):
    """text: a string, or a list of (string, style) or (string, style, url) with style in {'', 'b', 'i', 'bi'}."""
    for t, st, *url in ([(text, '')] if isinstance(text, str) else text):
        for k, part in enumerate(t.split('\n')):      #\n = a line break within the paragraph
            if k > 0:
                p.add_line_break()
            r = p.add_run()
            r.text = part
            r.font.name, r.font.size, r.font.color.rgb = font, Pt(size), color
            r.font.bold, r.font.italic = 'b' in st, 'i' in st
            if url:
                r.hyperlink.address = url[0]

def bullet(p, level):
    """a real PowerPoint round bullet, indented by level (0, 1, 2)."""
    pPr = p._p.get_or_add_pPr()
    pPr.set('marL', str(int(Inches(0.3 + 0.4*level))))
    pPr.set('indent', str(int(-Inches(0.25))))
    for tag in ('a:buFont', 'a:buChar'):
        el = pPr.find(qn(tag))
        if el is not None:
            pPr.remove(el)
    etree.SubElement(pPr, qn('a:buFont')).set('typeface', 'Arial')
    etree.SubElement(pPr, qn('a:buChar')).set('char', '•')

def bullets(slide, x, y, w, h, items, size=22):
    """items: (level, text) pairs; level-0 items get extra space above."""
    tf = textbox(slide, x, y, w, h)
    for k, (level, text) in enumerate(items):
        p = tf.paragraphs[0] if k == 0 else tf.add_paragraph()
        bullet(p, level)
        p.space_before = Pt(0 if k == 0 else (16 if level == 0 else 6))
        runs(p, text, size - 4*level)

def header(slide, text, size=28, x=0.5):
    tf = textbox(slide, x, 0.35, 9.5 - x, 0.6)
    runs(tf.paragraphs[0], [(text, 'b')], size)

def plain(slide, x, y, w, h, lines, size=11, font=MONO, color=BLACK, framed=False):
    tf = textbox(slide, x, y, w, h, framed)
    for k, ln in enumerate(lines):
        p = tf.paragraphs[0] if k == 0 else tf.add_paragraph()
        runs(p, ln, size, font=font, color=color)

def picture(slide, path, x, y, w=None, h=None, framed=False):
    pic = slide.shapes.add_picture(path, Inches(x), Inches(y), Inches(w) if w else None, Inches(h) if h else None)
    if framed:
        pic.line.color.rgb = GREY
        pic.line.width = Pt(0.75)

def notes(slide, text):
    slide.notes_slide.notes_text_frame.text = text

prs = Presentation()
prs.slide_width, prs.slide_height = Inches(10), Inches(7.5)
blank = prs.slide_layouts[6]

def new_slide():
    s = prs.slides.add_slide(blank)
    s.background.fill.solid()
    s.background.fill.fore_color.rgb = YELLOW
    return s

#1. what Claude Code is
s = new_slide()
tf = textbox(s, 0.5, 0.2, 9.0, 0.6)
tf.paragraphs[0].alignment = 2                   #center
runs(tf.paragraphs[0], [('Using Anthropic’s ', ''), ('Claude Code', 'b'), (' for Research', '')], 32)
plain(s, 0.3, 0.85, 4.7, 0.4, ['Joe Hahn, Space Science Institute'], size=20, font=FONT, color=GREEN)
tf = textbox(s, 5.0, 0.85, 4.7, 0.4)
tf.paragraphs[0].alignment = 3                   #right
runs(tf.paragraphs[0], 'SSI AI Working Group, Oct 28, 2026', 20, color=NAVY)
bullets(s, 0.3, 1.35, 5.3, 6.0, [
    (0, [('Claude Code (CC) is an AI ', ''), ('agent', 'i'), ('\nthat works inside my project folder', '')]),
    (1, 'it reads my files, writes & runs code,\nand I approve each step'),
    (1, 'not a chat window'),
    (0, 'context is everything:'),
    (1, 'CC reads my code, notes, papers & data\nbefore it answers,\nso its outputs align with the project’s goals'),
    (1, 'it also manages all new code development,\nthe paper & this deck,\nand will eventually draft my next funding proposal'),
    (0, 'I’m the advisor, CC is my grad student:'),
    (1, 'CC is very good at math, code & data analysis'),
    (1, 'but needs direction: derive it this way not that way\nand test, test, test'),
    (1, 'like a student it can drift,\nand rarely says “I don’t know”')], size=19)
plain(s, 5.6, 1.45, 4.1, 0.4, [[('a CC exchange from this project, condensed:', 'i')]], size=14, font=FONT,
      color=GREY)
plain(s, 5.6, 1.85, 4.1, 2.8, [
    '> should "Δ² = d_r² + d_z² − c(d_r s_r + d_z s_z)" be written delta^2=dr^2+dz^2-(cN)^2 instead?',
    '● Yes: they are equal, and the new form is just Pythagoras, |d0|² = Δ² + |d1|².',
    '● edits Eqn. 21 in the paper',
    '● writes sympy code that re-derives Δ² from the geometry, runs it:',
    '  PASS G16: Δ² = d_r² + d_z² − (Nc)²',
    '● recompiles the paper'], size=12, framed=True)
plain(s, 5.6, 4.85, 4.1, 0.3, [[('paper Eqn. (21), before:', 'i')]], size=14, font=FONT, color=GREY)
picture(s, 'assets/eq21_before.png', 5.65, 5.15, w=4.0)
plain(s, 5.6, 5.8, 4.1, 0.3, [[('after:', 'i')]], size=14, font=FONT, color=GREY)
picture(s, 'assets/eq21_after.png', 5.65, 6.1, w=4.0)
notes(s, 'Claude Code is an AI agent that works inside my project folder: it reads my files, papers and data, '
      'writes and runs code, and I approve or redirect each step. It is not a chat window you copy and paste '
      'from; it takes actions, and I approve them.\n\nContext is everything: before it answers, it reads my '
      'code, notes, papers and data, so its work fits my project instead of being a generic guess. The better '
      'my files, the better its work.\n\nI am the research advisor and CC is my grad student: very good at math, '
      'code and data analysis, but it needs direction (derive it this way, not that way; why this and not '
      'that?) and it has to be told to test, test, test. Unlike a student it never tires, but it forgets '
      'between sessions unless it is written in the project files, and it rarely says I don’t know: it sounds '
      'just as confident when wrong.')

#2. literature
s = new_slide()
#upper left: header & first bullet
header(s, 'Literature search & citations:', size=26)
bullets(s, 0.4, 0.95, 4.55, 1.8, [
    (0, 'CC reads the papers in my project folder & finds what I need'),
    (1, 'eg the wire-gravity geometry in Fig. 3 of Hahn & Spitale (2013), which CC redrew & extended to 3D')],
    size=20)
picture(s, 'assets/fig1_transp.png', 5.7, 0.2, h=2.85)       #upper right: the new paper's redrawn figure
bullets(s, 0.4, 2.95, 9.2, 1.7, [
    (0, 'but AI can invent citations'),
    (1, 'CC confirms each citation’s authors, year, journal, volume & page in NASA ADS'),
    (0, 'what CC knows from training vs what it reads in my folder:'),
    (1, 'from training: well-known results & widely cited papers,\n'
        'eg the viscous overstability papers it cited when we chased a numerical instability,\n'
        'none of which were in my folder; all four checked out in ADS'),
    (1, 'what’s in my project folder: my papers, code & notation, code output, notes,\n'
        'all the specifics that CC will need to advance the research'),
    (1, 'training knowledge is broad but fuzzy on details,\n'
        'so I copy the key papers into my project folder, where CC reads them')], size=22)
notes(s, 'The four overstability papers: Schmit & Tscharnuter 1995, Schmidt et al. 2001, Latter & Ogilvie '
      '2008, Lehmann et al. 2017. Each bibliography entry carries a line like: % SciX check: done 2026-10-20, '
      'ADS bibcode 2025ApJ...986..220H.\n\nLanguage models can invent plausible-looking citations. So every new citation is checked against '
      'SciX / NASA ADS through the ADS API, and the entry records the bibcode and the date of the check.')

#3. physics
s = new_slide()
header(s, 'Physics & equations of motion:')
bullets(s, 0.4, 1.05, 5.25, 6.3, [
    (0, 'who does what:'),
    (1, 'me: I set up the problem to be solved'),
    (1, 'CC: then does the physical & mathematical reasoning'),
    (1, 'sympy (Python symbolic math): CC calls sympy, which re-derives each equation from its definitions & confirms the result'),
    (1, 'me: I check the equations by hand &\nsuggest alternative derivations'),
    (0, 'eg I asked CC: “I suspect d₁ ≪ Δ, try to derive a simple proof”'),
    (1, 'CC reasoned it through: |d₁| ≤ sΔ,\nwhere s ≪ 1 is the streamline’s slope'),
    (1, 'CC then put the proof in the paper &\nadded sympy checks'),
    (0, [('CC also reviewed my predecessor code\n', ''), ('epi_int_lite', 'i'),
         (': it found a few small bugs and two approximations worth improving upon', '')])], size=19)
plain(s, 5.75, 1.15, 4.0, 0.35, [[('Eqn. 7: solved via recall (a textbook result)', 'i')]], size=14, font=FONT,
      color=GREY)
picture(s, 'assets/eqn7.png', 5.8, 1.5, w=3.9, framed=True)
plain(s, 5.75, 3.4, 4.0, 0.6, [[('Eqns. 19–21: stepwise reasoning about the geometry (figure on previous page)',
                                 'i')]], size=14, font=FONT, color=GREY)
picture(s, 'assets/eqn19_21.png', 5.9, 4.0, h=3.25, framed=True)
notes(s, 'Claude derives; every equation gets a sympy check before it is used in code, tracked in a table; I '
      'hand-verify them later. The sympy check rebuilds the result from its definitions as symbols (eg solves for '
      'the near point, forms Delta^2) and asserts that it minus the paper formula simplifies to exactly zero. It '
      'confirms the algebra, not the physical assumptions; that is my job. One of the 9 issues it found in epi_int_lite changes the torque between '
      'streamlines by a factor of about m for the m >= 2 spiral waves I am after.')

#4. simulations
s = new_slide()
header(s, 'CC manages numerical simulations, tested against theory:', size=26)
bullets(s, 0.4, 1.2, 4.4, 6.0, [
    (0, 'CC ported my old Python 2 code to Python 3.14, then made it ~35× faster'),
    (1, [('using ', ''), ('Numba', 'i'), (', which compiles Python like Fortran or C & parallelizes it across all CPU cores; identical results', '')]),
    (0, 'every test simulation is compared to theory'),
    (1, 'viscous spreading, libration rates, nodal regression'),
    (0, 'regression tests flag any change in the answers, to 1 part in 10¹¹')], size=20)
bullets(s, 0.4, 5.9, 9.2, 1.0, [
    (0, 'trust in the derivations & code is earned by comparing simulations to\nwell-understood theory')], size=20)
picture(s, 'assets/libration.png', 4.9, 1.6, w=4.9)
plain(s, 5.1, 4.65, 4.6, 0.8, ['simulated ringlets (dots) vs Borderies, Goldreich & Tremaine theory (dashed)'],
      size=18, font=FONT)
notes(s, '447 ms to 12 ms per timestep at 1000 streamlines. In the plot the simulations follow theory except at '
      'the lowest masses, where the libration amplitude is large. Regression tests compare every final state '
      'with stored answers to 1 part in 10^11.')

#5. writing
s = new_slide()
picture(s, 'assets/paper_p3.png', -4.7, 0, h=11.7)  # zoomed, runs off the left & bottom edges; behind the text
header(s, 'Writing the science:', x=4.05)
bullets(s, 4.05, 1.1, 5.75, 5.0, [
    (0, 'CC manages the LaTeX/AASTeX, bibliography & the code that makes the figures'),
    (0, 'CC drafted a style guide from my last two papers, refined it with my edits,\n& writes in my voice'),
    (0, 'I edit by pasting text from the PDF, eg “shorten …”'),
    (1, 'before: …so the particles that bracket the perturbed particle’s longitude θ are easily identified: '
        'they are that particle and its neighbor on either side.'),
    (1, 'after: …so it and its neighbor on either side bracket θ.'),
    (0, 'my judgment, my voice; its typing')], size=20)
bullets(s, 4.05, 5.95, 5.75, 1.3, [
    (0, 'advice: start the research project by writing a project plan, verify everything along the way, '
        'and save every prompt & change')], size=20)
notes(s, 'It writes the plotting scripts, compiles the PDF and fixes the warnings. I edit the way I would with a '
      'co-author; the example is a real edit from this week, now on page 3 of the draft.')

#6. about me, added after the talk for the public (LinkedIn) version
s = new_slide()
header(s, 'About me:')
bullets(s, 0.4, 1.1, 9.2, 5.0, [
    (0, [('Joe Hahn', 'b'), (', planetary scientist at the ', ''),
         ('Space Science Institute', '', 'https://www.spacescience.org/bio.php?emp=JHAHN')]),
    (1, 'research interests: planetary rings, planetary dynamics & the zodiacal cloud;\n'
        'planet formation & circumstellar disks; computational dynamics'),
    (0, [('also an AI consultant & builder at ', ''), ('JMH DataSciences', 'b')]),
    (1, '12+ years building production AI & machine learning,\n'
        'eight at Oracle, with the last three in its AI Center of Excellence'),
    (1, 'I now run my own independent AI consulting shop:'),
    (2, 'I build AI tools that do work a machine should do'),
    (2, '& provide senior hands to AI teams that have more work than people'),
    (0, 'this project’s code, draft paper & every prompt I gave CC:'),
    (1, 'github.com/joehahn/epi_int_3d')], size=22)
tf = textbox(s, 0.4, 6.35, 9.2, 0.6)
tf.paragraphs[0].alignment = 2                   #center
runs(tf.paragraphs[0], [('jmh-datasciences.com', 'b'), ('     joe.hahn@jmh-datasciences.com', '')], 24, color=GREEN)
notes(s, 'Added after the talk, for the version posted publicly.')

prs.save('Using_Claude_Code_for_Research.pptx')
print('wrote Using_Claude_Code_for_Research.pptx')
