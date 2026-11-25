"""Build the paper PDF from paper/manuscript.md + results/results.json.

Every {{ expression }} in the Markdown is evaluated against RESULTS (R) and the helpers below, so no number is
typed by hand. <<NAME>> placeholders are filled from paper/sections/NAME.md when that file exists.
Then: pandoc (+ Nature CSL citations) -> HTML -> headless Edge/Chrome -> one PDF with the supplement and code.
Usage: python paper/build_paper.py"""
import glob, html, json, os, pathlib, re, subprocess, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
PAPER, BUILD = "paper", "paper/build"
os.makedirs(BUILD, exist_ok=True)
R = json.load(open("results/results.json"))


# ---- number formatting helpers used inside {{ }} -------------------------------------------------
def f0(x): return f"{x:.0f}"
def f1(x): return f"{x:.1f}"
def f2(x): return f"{x:.2f}"
def pc(x): return f"{100 * x:.0f}%"
def pct(k, n): return f"{int(k)} of {int(n)} ({100 * k / n:.0f}%)"

def P(p):
    """Nature-style P value: 'P = 1.1 × 10^−9^' (pandoc superscript), or 'P = 0.51' when not tiny."""
    if p >= 0.001:
        return f"*P* = {p:.2g}" if p < 0.1 else f"*P* = {p:.2f}"
    m, e = f"{p:.1e}".split("e")
    return f"*P* = {m} × 10^−{int(e[1:])}^"

def allele(hgvs):
    """The a2 allele record whose name contains this HGVS c. string."""
    hits = [a for a in R["a2_alleles"] if f" {hgvs} " in a["allele"] + " "]
    assert len(hits) == 1, f"{hgvs}: {len(hits)} matches"
    return hits[0]

def allele_n(hgvs): return allele(hgvs)["n_strains"]
def allele_gain(hgvs): return allele(hgvs)["parsimony_changes"]

def struct(gene, hgvs):
    """Structure-context row for one allele (single-allele rows only)."""
    rows = [r for r in R["core_structure"] if r["gene"] == gene and r["allele"].startswith(hgvs + " ")]
    assert rows, f"no structure row for {gene} {hgvs}"
    return rows[0]

def struct_keep(hgvs): return struct("FDC1" if hgvs != "c.305G>A" else "PAD1", hgvs)["fraction_retained"]
def struct_lost(gene, hgvs): return struct(gene, hgvs)["n_binding_residues_lost"]
def struct_tot(gene, hgvs): return struct(gene, hgvs)["n_binding_residues_total"]

def sig(label):
    """'k of n' strains with this allele signature that were measured POF-."""
    d = R["a3_allele_pheno"][label]
    return f"{d['POF-']} of {d['n']}"

def late_removed():
    """Largest fraction of Fdc1 removed by a late stop that is not called LOF."""
    return pc(1 - min(r["fraction_retained"] for r in R["core_structure"] if not r["LOF_called"]))

def ale_clade():
    """Best-supported bootstrap clade that is all POF- and carries Q154*."""
    c = [d for d in R["core_boot_clades"] if set(d["POF"]) == {"POF-"} and any("Q154" in a for a in d["LOF_alleles"])]
    return max(c, key=lambda d: d["n_strains"])

def sake_clade():
    """Best-supported bootstrap clade that is all POF- and mostly saké."""
    c = [d for d in R["core_boot_clades"] if set(d["POF"]) == {"POF-"} and d["groups"].get("Sake", 0) > d["n_strains"] / 2]
    return max(c, key=lambda d: d["support"])


def md_table(df):
    """Pandas DataFrame -> pandoc pipe table."""
    cols = list(df.columns)
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    out += ["| " + " | ".join(str(v) for v in row) + " |" for row in df.itertuples(index=False)]
    return "\n".join(out)

def origins_table():
    import pandas as pd
    d = pd.DataFrame(R["a2_alleles"])
    d["called LOF"] = d.called_LOF.map({True: "yes", False: "no (late stop)"})
    d["allele"] = d.allele.str.replace("*", r"\*", regex=False)
    return md_table(d[["allele", "n_strains", "n_haplotypes", "parsimony_changes", "called LOF"]]
                    .rename(columns={"n_strains": "strains", "n_haplotypes": "haplotypes", "parsimony_changes": "min. gains"}))

def sens_table():
    import pandas as pd
    d = pd.DataFrame(R["a3_sensitivity_rules"])
    for c in ["sensitivity", "specificity", "accuracy"]:
        d[c] = d[c].map(pc)
    return md_table(d[["rule", "n", "TP", "FN", "TN", "FP", "sensitivity", "specificity", "accuracy"]])

def cmh_table():
    import pandas as pd
    d = pd.read_csv("data/results/cmh_strata_by_clade.csv")
    d["used in CMH"] = d.informative.map({True: "yes", False: "no (one side only)"})
    return md_table(d.drop(columns="informative").rename(columns={"dom_POFneg": "dom. POF−", "dom_POFpos": "dom. POF+",
                                                                   "ref_POFneg": "ref. POF−", "ref_POFpos": "ref. POF+"}))


def render(text):
    """Fill <<NAME>> from paper/sections/NAME.md, then evaluate every {{ expr }}."""
    def fill(m):
        path = f"{PAPER}/sections/{m.group(1)}.md"
        return open(path, encoding="utf-8").read().strip() if os.path.exists(path) else f"[[PENDING: {m.group(1)}]]"
    for _ in range(2):                       # sections may themselves contain {{ }} tokens
        text = re.sub(r"<<([A-Z0-9_]+)>>", fill, text)
    env = {k: v for k, v in globals().items() if callable(v) or k == "R"}
    return re.sub(r"\{\{(.+?)\}\}", lambda m: str(eval(m.group(1), env)), text)


def literal_numbers(src):
    """Digits in the Markdown source outside {{ }}, citations and code — for the hand-check report."""
    body = re.sub(r"\{\{.+?\}\}|\[@[^\]]+\]|`[^`]*`|\(figures[^)]*\)|\{[^}]*\}", "", src)
    hits = []
    for i, line in enumerate(body.splitlines(), 1):
        for m in re.finditer(r"(?<![A-Za-z.\-*_/])\d[\d.,]*%?", line):
            hits.append((i, m.group(), line[max(0, m.start() - 40): m.end() + 30].strip()))
    return hits


def code_appendix():
    """HTML listing of every custom code file (Supplementary Material)."""
    files = sorted(glob.glob("src/yeastflavor/*.py")) + ["run_all.py", "paper/build_paper.py"] + sorted(glob.glob("tests/*.py"))
    parts = ["<section class='codeapp'><h2>Supplementary Code</h2><p>All custom code, as run to produce this paper "
             "(also in the GitHub repository). Steps are notebook-style cells run in order by <code>run_all.py</code>.</p>"]
    for f in files:
        parts.append(f"<h3>{html.escape(f)}</h3><pre>{html.escape(open(f, encoding='utf-8').read())}</pre>")
    return "\n".join(parts) + "</section>"


def find_browser():
    for b in [r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
              r"C:\Program Files\Google\Chrome\Application\chrome.exe", "chromium", "google-chrome"]:
        if os.path.exists(b) or b in ("chromium", "google-chrome"):
            return b


def main():
    src = open(f"{PAPER}/manuscript.md", encoding="utf-8").read()
    sup = open(f"{PAPER}/supplementary.md", encoding="utf-8").read() if os.path.exists(f"{PAPER}/supplementary.md") else ""
    md = render(src) + "\n\n" + render(sup)
    open(f"{BUILD}/manuscript_rendered.md", "w", encoding="utf-8").write(md)

    lits = literal_numbers(src + "\n" + sup)
    with open(f"{BUILD}/literal_numbers.txt", "w", encoding="utf-8") as fh:
        fh.write("Numbers typed directly in the Markdown (not from results.json). Each must be a literature value,\n"
                 "a method parameter, an accession, a section/figure number or a residue position.\n\n")
        for i, n, ctx in lits:
            fh.write(f"line {i}: {n:>8}   ...{ctx}...\n")
    pending = sorted(set(re.findall(r"\[\[(?:PENDING|TEAM|VERIFY)[^\]]*\]\]", md)))

    # pandoc: markdown -> standalone HTML with numbered Nature-style citations, resources embedded
    subprocess.run(["pandoc", f"{BUILD}/manuscript_rendered.md", "-o", f"{BUILD}/paper.html", "-s", "--citeproc",
                    "--bibliography", f"{PAPER}/references.bib", "--csl", f"{PAPER}/nature.csl",
                    "--css", f"{PAPER}/style.css", "--embed-resources", "--resource-path", os.pathsep.join([".", PAPER, BUILD]),
                    "--metadata", "link-citations=true"], check=True)
    page = open(f"{BUILD}/paper.html", encoding="utf-8").read()
    page = page.replace("</body>", code_appendix() + "</body>")
    open(f"{BUILD}/paper.html", "w", encoding="utf-8").write(page)

    out = os.path.abspath(f"{PAPER}/BCH339N_Sperling_brewing_yeast_flavor.pdf")
    if os.path.exists(out):
        os.remove(out)
    # Edge/Chrome headless printing; the launcher can return before the printing process writes the file
    r = subprocess.run([find_browser(), "--headless", "--disable-gpu", "--no-pdf-header-footer", f"--print-to-pdf={out}",
                        pathlib.Path(f"{BUILD}/paper.html").resolve().as_uri()],
                       capture_output=True, text=True, timeout=300)
    for _ in range(120):                     # Edge's launcher returns before the printing child finishes
        if os.path.exists(out) and os.path.getsize(out) > 0 and time.time() - os.path.getmtime(out) > 2:
            break
        time.sleep(1)
    assert os.path.exists(out), f"browser did not write the PDF: {r.stderr[-800:]}"
    n_pages = len(re.findall(rb"/Type\s*/Page[^s]", open(out, "rb").read()))
    print(f"wrote {out} ({n_pages} pages incl. supplement and code)")
    print(f"literal numbers to hand-check: {len(lits)} (paper/build/literal_numbers.txt)")
    print("open placeholders:", *pending, sep="\n  ") if pending else print("no open placeholders")


if __name__ == "__main__":
    main()
