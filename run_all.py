"""One command: data -> analysis -> figures -> results.json -> executed notebook.

Runs every src/yeastflavor/NN_*.py step in order in one shared namespace (the steps were
written as notebook cells, '# %%' markers), collects outputs into notebooks/*.ipynb, and
dumps every number the paper cites (the RESULTS dict the steps fill) to results/results.json.
Usage:  python run_all.py            (cache-first; set NCBI_EMAIL for live fetches)
        python run_all.py --steps 07 08   (only steps whose file name starts with these)
shortcut: in-process exec instead of a Jupyter kernel; re-run the notebook in Jupyter for a kernel-native copy."""
import ast, base64, contextlib, glob, io, json, os, sys, traceback
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.abspath(__file__))
os.chdir(ROOT)


def split_cells(path):
    cells, kind, buf = [], None, []
    def flush():
        if kind and "".join(buf).strip():
            src = buf if kind == "code" else [l[2:] if l.startswith("# ") else l.lstrip("#") for l in buf]
            cells.append((kind, "".join(src).strip("\n")))
    for line in open(path, encoding="utf-8"):
        if line.startswith("# %%"):
            flush(); buf = []
            kind = "markdown" if "[markdown]" in line else "code"
        elif kind:
            buf.append(line)
    flush()
    return cells


def run(src, ns):
    outs, stdout = [], io.StringIO()
    tree = ast.parse(src)
    last = tree.body.pop() if tree.body and isinstance(tree.body[-1], ast.Expr) else None
    try:
        with contextlib.redirect_stdout(stdout):
            exec(compile(tree, "<cell>", "exec"), ns)
            val = eval(compile(ast.Expression(last.value), "<cell>", "eval"), ns) if last else None
    except Exception:
        val = None
        outs.append({"output_type": "error", "ename": "Error", "evalue": "",
                     "traceback": traceback.format_exc().splitlines()})
        print(traceback.format_exc(), file=sys.stderr)
    if stdout.getvalue():
        outs.insert(0, {"output_type": "stream", "name": "stdout", "text": stdout.getvalue()})
    for num in plt.get_fignums():
        b = io.BytesIO(); plt.figure(num).savefig(b, format="png", dpi=110, bbox_inches="tight")
        outs.append({"output_type": "display_data", "metadata": {},
                     "data": {"image/png": base64.b64encode(b.getvalue()).decode(), "text/plain": "<Figure>"}})
    plt.close("all")
    if val is not None and not isinstance(val, matplotlib.artist.Artist):
        data = {"text/plain": repr(val)}
        if hasattr(val, "to_html"):
            data["text/html"] = val.to_html(max_rows=40)
        outs.append({"output_type": "execute_result", "metadata": {}, "data": data, "execution_count": None})
    return outs


def jsonable(o):
    # numpy scalars / arrays -> plain python so json.dump works
    if hasattr(o, "item") and getattr(o, "ndim", 1) == 0:
        return o.item()
    if hasattr(o, "tolist"):
        return o.tolist()
    return str(o)


def main(steps=None, out="notebooks/BCH339N_Brewing_Yeast_Flavor.ipynb"):
    plt.show = lambda *a, **k: None  # figures are collected after each cell
    os.makedirs("results", exist_ok=True); os.makedirs("notebooks", exist_ok=True)
    ns = {"__name__": "__main__", "RESULTS": {}}
    paths = sorted(glob.glob("src/yeastflavor/[0-9][0-9]_*.py"))
    if steps:  # partial run (dev only): step 00 + the named steps; reads earlier steps' outputs from data/results/
        paths = [p for p in paths if os.path.basename(p)[:2] in set(steps) | {"00"}]
    nb_cells, n, errors = [], 0, 0
    for path in paths:
        for kind, src in split_cells(path):
            if kind == "markdown":
                nb_cells.append({"cell_type": "markdown", "metadata": {}, "source": src})
                continue
            n += 1
            outs = run(src, ns)
            bad = any(o["output_type"] == "error" for o in outs)
            errors += bad
            for o in outs:
                if o["output_type"] == "execute_result":
                    o["execution_count"] = n
            nb_cells.append({"cell_type": "code", "metadata": {}, "execution_count": n,
                             "source": src, "outputs": outs})
            print(f"[{n}] {path}: {'ERROR' if bad else 'ok'}", flush=True)
    # partial runs never touch results.json, so parallel step development can't clobber it
    res_path = "results/results.json" if not steps else f"results/partial_{'_'.join(steps)}.json"
    json.dump(ns["RESULTS"], open(res_path, "w"), indent=1, sort_keys=True, default=jsonable)
    if not steps:
        nb = {"cells": nb_cells, "nbformat": 4, "nbformat_minor": 5,
              "metadata": {"kernelspec": {"name": "python3", "display_name": "Python 3", "language": "python"},
                           "language_info": {"name": "python"}}}
        json.dump(nb, open(out, "w"), indent=1)
        print(f"wrote {out}: {len(nb_cells)} cells")
    print(f"{res_path}: {len(ns['RESULTS'])} keys; {errors} cells with errors")
    return errors


if __name__ == "__main__":
    args = sys.argv[1:]
    sys.exit(main(args[args.index("--steps") + 1:] if "--steps" in args else None))
