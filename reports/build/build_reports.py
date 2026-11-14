"""Build both reports (PDF via pdflatex, editable DOCX via pandoc) from their LaTeX sources.

    python reports/build/build_reports.py

Requires: a TeX Live installation with pdflatex (IEEEtran for the paper),
pandoc >= 3, and python-docx. Run `python run_all.py` first: the reports read
every number from reports/build/generated/values.tex and every table from
reports/build/generated/*.tex.

The DOCX files are produced from the same LaTeX source. Equation, table and
figure numbers in the DOCX are taken from the LaTeX .aux file so they match
the PDF exactly.
"""
from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt, Inches, RGBColor

BUILD = Path(__file__).resolve().parent
ROOT = BUILD.parents[1]
AUTHOR = "Nithin D H"
INSTITUTION = "Department of Chemical Engineering, Indian Institute of Technology Kanpur"
KEYWORDS = ("jacketed stirred tank; RK4; Monte Carlo; chance-constrained optimisation; Sobol sensitivity; "
            "adaptive operation; simulation study")

# Title and subject must match the \hypersetup block in each LaTeX source (tests/test_release.py checks both).
REPORTS = {
    "detailed_report": dict(src=ROOT / "reports" / "detailed_report" / "source" / "detailed_report.tex",
                            out=ROOT / "reports" / "detailed_report", columns=1,
                            header="Detailed Technical Report \u2014 Jacketed Stirred-Tank Cooling", font="Cambria", size=11,
                            title="Dynamic Modelling, Reliability-Based Optimisation and Adaptive Operation of a "
                                  "Jacketed Stirred-Tank Cooling Process",
                            subject=f"Detailed technical report, ESC113M extended project (V3), {INSTITUTION}"),
    "scientific_paper": dict(src=ROOT / "reports" / "scientific_paper" / "source" / "scientific_paper.tex",
                             out=ROOT / "reports" / "scientific_paper", columns=2, header=None, font="Times New Roman", size=10,
                             title="Reliability-Based Optimisation and Adaptive Operation of a Jacketed Stirred-Tank "
                                   "Cooling Process",
                             subject=f"Scientific paper, ESC113M extended project (V3), {INSTITUTION}"),
}


def run(cmd, cwd):
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-3000:], r.stderr[-2000:])
        raise SystemExit(f"command failed: {' '.join(cmd)}")
    return r


def stage(name, cfg):
    work = BUILD / "_work" / name
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    shutil.copy(cfg["src"], work / cfg["src"].name)
    shutil.copytree(BUILD / "generated", work / "generated")
    shutil.copytree(ROOT / "figures", work / "figures")
    return work


def pdflatex(work: Path, tex: str):
    for _ in range(3):
        run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", tex], work)
    log = (work / tex.replace(".tex", ".log")).read_text(errors="ignore")
    undefined = re.findall(r"(Reference|Citation) `[^']+' on page \d+ undefined", log)
    if undefined:
        raise SystemExit(f"undefined references: {undefined[:5]}")
    return work / tex.replace(".tex", ".pdf")


# ---------------------------------------------------------------- LaTeX -> pandoc-friendly LaTeX
def read_labels(aux: Path):
    labels = {}
    for line in aux.read_text(errors="ignore").splitlines():
        m = re.match(r"\\newlabel\{([^}]*)\}\{\{([^}]*)\}", line)
        if m:
            labels[m.group(1)] = m.group(2)
        m = re.match(r"\\bibcite\{([^}]*)\}\{([^}]*)\}", line)
        if m:
            labels["cite:" + m.group(1)] = m.group(2)
    return labels


def braced(s: str, i: int):
    """Return (content, end) for the brace group starting at s[i] == '{'."""
    depth = 0
    for j in range(i, len(s)):
        if s[j] == "{":
            depth += 1
        elif s[j] == "}":
            depth -= 1
            if depth == 0:
                return s[i + 1:j], j + 1
    raise ValueError("unbalanced braces")


def inline_inputs(tex: str, base: Path):
    return re.sub(r"\\input\{([^}]+)\}", lambda m: (base / (m.group(1) if m.group(1).endswith(".tex") else m.group(1) + ".tex")).read_text(), tex)


def for_pandoc(tex: str, labels: dict, kind: str) -> str:
    # keep the macro file, inline generated tables
    # equations: append the LaTeX number, remove the environment numbering
    def eq_repl(m):
        body = m.group(1)
        lab = re.search(r"\\label\{([^}]*)\}", body)
        num = labels.get(lab.group(1), "") if lab else ""
        body = re.sub(r"\\label\{[^}]*\}", "", body).strip()
        return f"\\[ {body} \\qquad\\qquad ({num}) \\]"
    def align_repl(m):
        rows = [r.strip() for r in re.split(r"\\\\(?![a-zA-Z])", m.group(1)) if r.strip()]
        out = []
        for r in rows:
            lab = re.search(r"\\label\{([^}]*)\}", r)
            num = labels.get(lab.group(1), "") if lab else ""
            r = re.sub(r"\\label\{[^}]*\}", "", r).replace("&", "").strip()
            out.append(f"\\[ {r} \\qquad\\qquad ({num}) \\]" if num else f"\\[ {r} \\]")
        return "\n".join(out)
    tex = re.sub(r"\\begin\{align\}(.*?)\\end\{align\}", align_repl, tex, flags=re.S)
    # aligned blocks inside one equation: put the rows on one line
    tex = re.sub(r"\\begin\{aligned\}(.*?)\\end\{aligned\}",
                 lambda m: " \\quad ".join(r.replace("&", "").strip() for r in re.split(r"\\\\(?![a-zA-Z])", m.group(1)) if r.strip()),
                 tex, flags=re.S)
    tex = re.sub(r"\\begin\{equation\}(.*?)\\end\{equation\}", eq_repl, tex, flags=re.S)
    tex = re.sub(r"\\eqref\{([^}]*)\}", lambda m: f"({labels.get(m.group(1), '?')})", tex)
    tex = re.sub(r"~?\\cite\{([^}]*)\}", lambda m: " [" + ", ".join(labels.get("cite:" + k.strip(), "?") for k in m.group(1).split(",")) + "]", tex)
    def bib(m):
        items = re.split(r"\\bibitem\{([^}]*)\}", m.group(1))[1:]
        body = "".join(f"\n\n[{labels.get('cite:' + items[i], '?')}] {items[i + 1].strip()}" for i in range(0, len(items), 2))
        return "\\section*{References}" + body + "\n\n"
    tex = re.sub(r"\\begin\{thebibliography\}\{[^}]*\}(.*?)\\end\{thebibliography\}", bib, tex, flags=re.S)
    tex = re.sub(r"\\ref\{([^}]*)\}", lambda m: labels.get(m.group(1), "?"), tex)
    # captions: prefix 'Table n:' / 'Figure n:' from the .aux numbers
    out, i = [], 0
    for m in re.finditer(r"\\caption\{", tex):
        if m.start() < i:
            continue
        content, end = braced(tex, m.end() - 1)
        lab = re.match(r"\s*\\label\{([^}]*)\}", tex[end:])
        num = labels.get(lab.group(1), "") if lab else ""
        env = tex.rfind("\\begin{table", 0, m.start()) > tex.rfind("\\begin{figure", 0, m.start())
        prefix = ("Table" if env else "Figure") + (f" {num}: " if num else ": ")
        out.append(tex[i:m.start()] + "\\caption{\\textbf{" + prefix + "}" + content + "}")
        i = end
    out.append(tex[i:])
    tex = "".join(out)
    # appendix headings (pandoc would continue the numbering)
    if "\\appendix" in tex:
        head, tail = tex.split("\\appendix", 1)
        letters = iter("ABCDEFG")
        tail = re.sub(r"\\section\{([^}]*)\}", lambda m: f"\\section*{{Appendix {next(letters)}. {m.group(1)}}}", tail)
        tex = head + tail
    tex = re.sub(r"\\noindent\\rule\{\\textwidth\}\{0\.4pt\}", "", tex)
    tex = re.sub(r"\\begin\{(enumerate|itemize)\}\[[^\]]*\]", lambda m: "\\begin{" + m.group(1) + "}", tex)
    if kind == "paper":
        tex = paper_frontmatter_for_pandoc(tex)
    else:
        tex = report_frontmatter_for_pandoc(tex)
    return tex


def _centre_block_to_paragraphs(body: str) -> str:
    body = re.sub(r"\\vspace\*?\{[^}]*\}", "", body)
    body = body.replace("\\begin{center}", "").replace("\\end{center}", "")
    body = re.sub(r"\\\\(\[[^\]]*\])?", "\n\n", body)
    return body


def report_frontmatter_for_pandoc(tex: str) -> str:
    """Cover page and page-1 title block -> marked paragraphs; postprocess() centres them,
    inserts the page break after the cover and the rule under the title block."""
    tex = re.sub(r"\\begin\{titlepage\}(.*?)\\end\{titlepage\}",
                 lambda m: "\n\nCOVERBEGIN\n\n" + _centre_block_to_paragraphs(m.group(1)) + "\n\nCOVEREND\n\n", tex, count=1, flags=re.S)
    tex = re.sub(r"\\setcounter\{page\}\{1\}\s*\\begin\{center\}(.*?)\\end\{center\}",
                 lambda m: "\n\nTITLEBEGIN\n\n" + _centre_block_to_paragraphs(m.group(1)) + "\n\nTITLEEND\n\n", tex, count=1, flags=re.S)
    return tex


def paper_frontmatter_for_pandoc(tex: str) -> str:
    """IEEEtran front matter -> plain LaTeX that pandoc understands."""
    author_block, _ = braced(tex, tex.index("\\author{") + len("\\author"))
    names = re.findall(r"\\IEEEauthorblockN\{(.*?)\}", author_block)
    affs = re.findall(r"\\IEEEauthorblockA\{(.*?)\}\s*\}?\s*$", author_block, flags=re.S)
    aff = re.sub(r"\\\\", ", ", affs[0]) if affs else ""
    start = tex.index("\\author{")
    _, end = braced(tex, start + len("\\author"))
    tex = tex[:start] + "\\author{" + ", ".join(names) + "}" + tex[end:]
    tex = re.sub(r"\\maketitle", lambda m: "\\maketitle\n\\begin{center}" + aff + "\\end{center}\n", tex)
    tex = re.sub(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", lambda m: "\\textbf{\\textit{Abstract}---}" + m.group(1).strip() + "\n\n", tex, flags=re.S)
    tex = re.sub(r"\\begin\{IEEEkeywords\}(.*?)\\end\{IEEEkeywords\}", lambda m: "\\textbf{\\textit{Keywords}---}" + m.group(1).strip() + "\n\nBODYSTARTMARKER\n\n", tex, flags=re.S)
    tex = tex.replace("\\begin{figure*}", "\\begin{figure}").replace("\\end{figure*}", "\\end{figure}").replace("\\balance", "")
    tex = tex.replace("\\begin{table*}", "\\begin{table}").replace("\\end{table*}", "\\end{table}")
    return tex


# ---------------------------------------------------------------- reference document and post-processing
def set_cell_border(cell, **edges):
    tcPr = cell._tc.get_or_add_tcPr()
    borders = tcPr.find(qn("w:tcBorders"))
    if borders is None:
        borders = OxmlElement("w:tcBorders"); tcPr.append(borders)
    for edge, (val, sz) in edges.items():
        el = borders.find(qn(f"w:{edge}"))
        if el is None:
            el = OxmlElement(f"w:{edge}"); borders.append(el)
        el.set(qn("w:val"), val); el.set(qn("w:sz"), str(sz)); el.set(qn("w:space"), "0"); el.set(qn("w:color"), "000000")


def add_page_field(paragraph):
    run = paragraph.add_run()
    for tag, text in (("begin", None), (None, "PAGE"), ("end", None)):
        if tag:
            el = OxmlElement("w:fldChar"); el.set(qn("w:fldCharType"), tag); run._r.append(el)
        else:
            el = OxmlElement("w:instrText"); el.set(qn("xml:space"), "preserve"); el.text = text; run._r.append(el)


def style(doc, name):
    """Case-insensitive style lookup (pandoc's reference.docx uses lower-case names for some styles)."""
    for s in doc.styles:
        if s.name and s.name.lower() == name.lower():
            return s
    return None


def reference_docx(path: Path, cfg):
    data = subprocess.run(["pandoc", "--print-default-data-file", "reference.docx"], capture_output=True).stdout
    path.write_bytes(data)
    doc = Document(str(path))
    font, size = cfg["font"], cfg["size"]
    for st in doc.styles:
        try:
            if st.type == 1:  # paragraph styles
                st.font.name = font
                rpr = st.element.get_or_add_rPr()
                rfonts = rpr.find(qn("w:rFonts"))
                if rfonts is None:
                    rfonts = OxmlElement("w:rFonts"); rpr.append(rfonts)
                for a in ("w:ascii", "w:hAnsi", "w:cs", "w:eastAsia"):
                    rfonts.set(qn(a), font)
                for a in ("w:asciiTheme", "w:hAnsiTheme", "w:eastAsiaTheme", "w:cstheme"):
                    rfonts.attrib.pop(qn(a), None)
                st.font.color.rgb = RGBColor(0, 0, 0)
        except Exception:
            pass
    for name in ("Normal", "Body Text", "First Paragraph", "Compact"):
        s = style(doc, name)
        if s is not None:
            s.font.size = Pt(size)
            s.paragraph_format.space_after = Pt(4); s.paragraph_format.space_before = Pt(0)
            if cfg["columns"] == 2:
                s.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    for name, sz in (("Heading 1", size + 3 if cfg["columns"] == 1 else size), ("Heading 2", size + 1 if cfg["columns"] == 1 else size),
                     ("Heading 3", size)):
        s = style(doc, name)
        if s is None:
            continue
        s.font.size = Pt(sz); s.font.bold = True; s.font.italic = False
        s.paragraph_format.space_before = Pt(10); s.paragraph_format.space_after = Pt(4)
    for name in ("Caption", "Table Caption", "Image Caption"):
        s = style(doc, name)
        if s is not None:
            s.font.size = Pt(size - 1.5); s.font.italic = False
    if style(doc, "Title") is not None:
        style(doc, "Title").font.size = Pt(18); style(doc, "Title").font.bold = True
    sec = doc.sections[0]
    sec.page_width, sec.page_height = Inches(8.5), Inches(11)
    m = Inches(1.0) if cfg["columns"] == 1 else Inches(0.75)
    sec.left_margin = sec.right_margin = sec.top_margin = sec.bottom_margin = m
    if cfg["header"]:
        hp = sec.header.paragraphs[0]
        hp.text = ""
        r = hp.add_run(cfg["header"]); r.font.size = Pt(9); r.font.name = font
        hp.add_run("\t")
        add_page_field(hp)
        pPr = hp._p.get_or_add_pPr()
        tabs = OxmlElement("w:tabs"); tab = OxmlElement("w:tab"); tab.set(qn("w:val"), "right"); tab.set(qn("w:pos"), "9360")
        tabs.append(tab); pPr.append(tabs)
        bdr = OxmlElement("w:pBdr"); b = OxmlElement("w:bottom")
        for k, v in (("w:val", "single"), ("w:sz", "4"), ("w:space", "1"), ("w:color", "000000")):
            b.set(qn(k), v)
        bdr.append(b); pPr.append(bdr)
    doc.save(str(path))


def _report_frontmatter(doc, cfg):
    paras = doc.paragraphs
    def block(begin, end):
        i = next((k for k, q in enumerate(paras) if q.text.strip() == begin), None)
        j = next((k for k, q in enumerate(paras) if q.text.strip() == end), None)
        return (i, j) if i is not None and j is not None else (None, None)
    for begin, end, top in (("COVERBEGIN", "COVEREND", True), ("TITLEBEGIN", "TITLEEND", False)):
        i, j = block(begin, end)
        if i is None:
            continue
        inner = [q for q in paras[i + 1:j] if q.text.strip()]
        for k, q in enumerate(inner):
            q.alignment = WD_ALIGN_PARAGRAPH.CENTER
            q.paragraph_format.first_line_indent = Pt(0)
            for r in q.runs:
                if k == 0:
                    r.font.size = Pt(18); r.font.bold = True
        if top:
            # the cover is its own section: no header, page numbers restart at 1 afterwards (as in the CHE453 formatting template; see detailed_report.tex)
            import copy
            inner[0].paragraph_format.space_before = Pt(220)
            final = doc.element.body.find(qn("w:sectPr"))
            cover = copy.deepcopy(final)
            for ref in cover.findall(qn("w:headerReference")) + cover.findall(qn("w:footerReference")):
                cover.remove(ref)
            typ = cover.find(qn("w:type"))
            if typ is None:
                typ = OxmlElement("w:type"); cover.insert(0, typ)
            typ.set(qn("w:val"), "nextPage")
            inner[-1]._p.get_or_add_pPr().append(cover)
            pg = final.find(qn("w:pgNumType"))
            if pg is None:
                pg = OxmlElement("w:pgNumType"); final.append(pg)
            pg.set(qn("w:start"), "1")
        else:
            pPr = inner[-1]._p.get_or_add_pPr()
            bdr = OxmlElement("w:pBdr"); b = OxmlElement("w:bottom")
            for kk, v in (("w:val", "single"), ("w:sz", "6"), ("w:space", "6"), ("w:color", "000000")):
                b.set(qn(kk), v)
            bdr.append(b); pPr.append(bdr)
        for q in (paras[i], paras[j]):
            q._p.getparent().remove(q._p)
        paras = doc.paragraphs


def postprocess(docx_path: Path, cfg):
    doc = Document(str(docx_path))
    # empty OMML separator characters make some renderers draw |..| instead of (..); Word ignores them
    for el in list(doc.element.body.iter(qn("m:sepChr"))):
        if el.get(qn("m:val")) == "":
            el.getparent().remove(el)
    # booktabs-style tables: top rule, rule under the header row, bottom rule, no vertical lines
    for t in doc.tables:
        rows = t.rows
        for ri, row in enumerate(rows):
            for cell in row.cells:
                edges = dict(left=("nil", 0), right=("nil", 0), top=("nil", 0), bottom=("nil", 0))
                if ri == 0:
                    edges["top"] = ("single", 12); edges["bottom"] = ("single", 6)
                if ri == len(rows) - 1:
                    edges["bottom"] = ("single", 12)
                set_cell_border(cell, **edges)
                for p in cell.paragraphs:
                    p.paragraph_format.space_after = Pt(0)
                    for r in p.runs:
                        r.font.size = Pt(cfg["size"] - 2 if cfg["columns"] == 1 else 7.5)
    for t in doc.tables:
        tblPr = t._tbl.tblPr
        w = tblPr.find(qn("w:tblW"))
        if w is None:
            w = OxmlElement("w:tblW"); tblPr.append(w)
        w.set(qn("w:type"), "pct"); w.set(qn("w:w"), "5000")
        lay = tblPr.find(qn("w:tblLayout"))
        if lay is None:
            lay = OxmlElement("w:tblLayout"); tblPr.append(lay)
        lay.set(qn("w:type"), "autofit")
        grid = t._tbl.tblGrid
        if grid is not None:
            for gc in grid.findall(qn("w:gridCol")):
                gc.attrib.pop(qn("w:w"), None)
        for row in t.rows:
            for cell in row.cells:
                tcW = cell._tc.get_or_add_tcPr().find(qn("w:tcW"))
                if tcW is not None:
                    tcW.set(qn("w:type"), "auto"); tcW.set(qn("w:w"), "0")
    # pictures: keep inside the text width / column width
    max_w = Inches(6.5) if cfg["columns"] == 1 else Inches(3.4)
    for shape in doc.inline_shapes:
        if shape.width > max_w:
            ratio = max_w / shape.width
            shape.width = int(shape.width * ratio); shape.height = int(shape.height * ratio)
    if cfg["columns"] == 1:
        _report_frontmatter(doc, cfg)
    if cfg["columns"] == 2:
        body = doc.element.body
        marker = next((p for p in doc.paragraphs if "BODYSTARTMARKER" in p.text), None)
        if marker is not None:
            # section break after the abstract/keywords: single column above, two columns below
            pPr = marker._p.get_or_add_pPr()
            sect = OxmlElement("w:sectPr")
            typ = OxmlElement("w:type"); typ.set(qn("w:val"), "continuous"); sect.append(typ)
            pg = OxmlElement("w:pgSz"); pg.set(qn("w:w"), "12240"); pg.set(qn("w:h"), "15840"); sect.append(pg)
            mar = OxmlElement("w:pgMar")
            for k in ("top", "bottom", "left", "right"):
                mar.set(qn(f"w:{k}"), "1080")
            for k in ("header", "footer", "gutter"):
                mar.set(qn(f"w:{k}"), "720" if k != "gutter" else "0")
            sect.append(mar)
            cols = OxmlElement("w:cols"); cols.set(qn("w:num"), "1"); sect.append(cols)
            pPr.append(sect)
            for r in marker.runs:
                r.text = ""
        final = body.find(qn("w:sectPr"))
        typ = final.find(qn("w:type"))
        if typ is None:
            typ = OxmlElement("w:type"); final.insert(0, typ)
        typ.set(qn("w:val"), "continuous")
        cols = final.find(qn("w:cols"))
        if cols is None:
            cols = OxmlElement("w:cols"); final.append(cols)
        cols.set(qn("w:num"), "2"); cols.set(qn("w:space"), "360")
    cp = doc.core_properties
    cp.title, cp.author, cp.subject, cp.keywords = cfg["title"], AUTHOR, cfg["subject"], KEYWORDS
    cp.last_modified_by, cp.comments, cp.category = AUTHOR, INSTITUTION, "Simulation study"
    doc.save(str(docx_path))


def build(name, cfg):
    work = stage(name, cfg)
    tex = cfg["src"].name
    pdf = pdflatex(work, tex)
    cfg["out"].mkdir(parents=True, exist_ok=True)
    shutil.copy(pdf, cfg["out"] / f"{name}.pdf")
    labels = read_labels(work / tex.replace(".tex", ".aux"))
    src = inline_inputs((work / tex).read_text(), work)
    pandoc_tex = for_pandoc(src, labels, "paper" if cfg["columns"] == 2 else "report")
    (work / "pandoc_input.tex").write_text(pandoc_tex)
    ref = work / "reference.docx"
    reference_docx(ref, cfg)
    out = cfg["out"] / f"{name}.docx"
    args = ["pandoc", "pandoc_input.tex", "-f", "latex", "-o", str(out), "--reference-doc", str(ref), "--resource-path", str(work)]
    if cfg["columns"] == 1:
        args.append("--number-sections")
    run(args, work)
    postprocess(out, cfg)
    print(f"built {name}: {cfg['out'] / (name + '.pdf')}, {out}")


if __name__ == "__main__":
    which = [a for a in sys.argv[1:] if not a.startswith("--")] or list(REPORTS)
    for n in which:
        build(n, REPORTS[n])
    if "--keep" not in sys.argv:
        shutil.rmtree(BUILD / "_work", ignore_errors=True)
