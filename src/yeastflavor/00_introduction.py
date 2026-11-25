# %%
# Hops, Hooch, and Hidden Genes
import os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, matplotlib.pyplot as plt, seaborn as sns
import Bio
from Bio import Entrez
Entrez.email = os.environ.get("NCBI_EMAIL", "anonymous@example.org")   # always tell NCBI who you are (lecture 09); set NCBI_EMAIL
os.makedirs("figures", exist_ok=True)
print("Biopython", Bio.__version__, "| pandas", pd.__version__, "| numpy", np.__version__)

# %%
# Cache-first fetch for the published supplementary tables (data/raw/ is not in git). I look the file up in
# data/MANIFEST.csv, download it only if it is missing, and check its SHA-256 so a changed upstream file can't
# silently change our results.
import hashlib, io, time, zipfile, requests

def ensure_external(path):
    """Make sure a data/raw/external file exists (download from its MANIFEST URL if not) and matches its checksum."""
    man = pd.read_csv("data/MANIFEST.csv")
    row = man[man.local_path == path]
    assert len(row) == 1, f"{path} is not listed in data/MANIFEST.csv"
    url, sha = row.accession_or_url.iloc[0], row.sha256.iloc[0]
    if not os.path.exists(path):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        # Europe PMC sometimes answers 500/503, so I retry a few times with a growing wait before giving up
        for attempt in range(5):
            r = requests.get(url.split(" ")[0], timeout=300)
            if r.ok:
                break
            time.sleep(15 * (attempt + 1))
        r.raise_for_status()
        if "europepmc" in url:            # Europe PMC serves all supplementary files of a paper as one zip
            z = zipfile.ZipFile(io.BytesIO(r.content))
            data = z.read(next(n for n in z.namelist() if n.endswith(os.path.basename(path))))
        else:
            data = r.content
        open(path, "wb").write(data)
        print(f"downloaded {path}")
    digest = hashlib.sha256(open(path, "rb").read()).hexdigest()
    assert digest == sha, f"{path}: checksum differs from MANIFEST ({digest[:12]} vs {sha[:12]})"
    return path
