"""Release checks: licensing files, untouched original materials, report metadata, README structure."""
import hashlib
import importlib.util
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
AUTHOR = "Nithin D H"
INSTITUTION = "Department of Chemical Engineering, Indian Institute of Technology Kanpur"
TITLES = {
    "detailed_report": "Dynamic Modelling, Reliability-Based Optimisation and Adaptive Operation of a "
                       "Jacketed Stirred-Tank Cooling Process",
    "scientific_paper": "Reliability-Based Optimisation and Adaptive Operation of a Jacketed Stirred-Tank Cooling Process",
}
# MD5 of the original ESC113M group files as received (V1). They must stay byte-identical.
ORIGINAL_MD5 = {
    "ESC113M_Project_Code.m": "63f71ba1e5b831951bdc3052d4420740",
    "Original_ESC113M_Project_Report.pdf": "32c397dfd5126a76bf746aab2f746d05",
}


def test_licence_files_are_complete():
    mit = (ROOT / "LICENSE").read_text(encoding="utf-8")
    assert mit.startswith("MIT License") and f"Copyright (c) 2026 {AUTHOR}" in mit
    assert "THE SOFTWARE IS PROVIDED \"AS IS\"" in mit
    cc = (ROOT / "LICENSE-CC-BY-4.0.txt").read_text(encoding="utf-8")
    assert "Attribution 4.0 International" in cc and "Section 8" in cc
    lic = (ROOT / "LICENSING.md").read_text(encoding="utf-8")
    assert "docs/original_esc113m/" in lic and "ESC113M Group 03" in lic and "Not covered by either licence" in lic


def test_original_group_materials_are_unchanged():
    for name, md5 in ORIGINAL_MD5.items():
        assert hashlib.md5((ROOT / "docs" / "original_esc113m" / name).read_bytes()).hexdigest() == md5, name


def test_report_metadata_is_consistent_across_sources():
    spec = importlib.util.spec_from_file_location("build_reports", ROOT / "reports" / "build" / "build_reports.py")
    br = importlib.util.module_from_spec(spec)
    pytest.importorskip("docx")
    spec.loader.exec_module(br)
    assert br.AUTHOR == AUTHOR and br.INSTITUTION == INSTITUTION
    for name, title in TITLES.items():
        assert br.REPORTS[name]["title"] == title
        tex = br.REPORTS[name]["src"].read_text(encoding="utf-8")
        hs = re.search(r"\\hypersetup\{(.*?)\}\}\s*\n", tex, re.S).group(1)
        assert f"pdftitle={{{title}}}" in hs and f"pdfauthor={{{AUTHOR}}}" in hs and INSTITUTION in hs


@pytest.mark.parametrize("name", list(TITLES))
def test_built_report_files_carry_metadata(name):
    docx = pytest.importorskip("docx")
    pypdf = pytest.importorskip("pypdf")
    meta = pypdf.PdfReader(ROOT / "reports" / name / f"{name}.pdf").metadata
    assert meta.title == TITLES[name] and meta.author == AUTHOR and meta.get("/Institution") == INSTITUTION
    cp = docx.Document(str(ROOT / "reports" / name / f"{name}.docx")).core_properties
    assert cp.title == TITLES[name] and cp.author == AUTHOR and cp.comments == INSTITUTION


def test_readme_has_the_release_structure():
    text = (ROOT / "README.md").read_text(encoding="utf-8")
    headings = re.findall(r"^## (.+)$", text, re.M)
    expected = ["Project overview", "Engineering problem", "Model extension", "Main findings", "Repository structure",
                "Installation", "Reproducing the results", "Testing", "Assumptions and limitations",
                "References and attribution", "Licence"]
    assert headings == expected
    for phrase in ("simulation", "assumed", "not been validated experimentally", "auxiliary"):
        assert phrase in text


def test_no_machine_specific_paths_in_tracked_text():
    bad = re.compile(r"(/home/claude|/mnt/user-data|/tmp/claude|C:\\Users\\)")
    suffixes = {".py", ".m", ".md", ".tex", ".csv", ".json", ".yml", ".txt", ".cff", ".bib"}
    hits = [str(p.relative_to(ROOT)) for p in ROOT.rglob("*")
            if p.is_file() and p.suffix in suffixes and "_work" not in p.parts and ".venv" not in p.parts
            and p.name != "test_release.py" and bad.search(p.read_text(encoding="utf-8", errors="ignore"))]
    assert hits == []
