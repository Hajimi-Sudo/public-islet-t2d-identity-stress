#!/usr/bin/env python3
"""Download public GEO matrices, GPL annot files, and GTEx v8 median TPM into cache/.

Does not download FASTQ, RAW.tar, NHANES, MIMIC, or GWAS Catalog dumps.
If NCBI/GTEx fail under a proxy, retry once with proxies unset (host policy).
"""
from __future__ import annotations

import os
import ssl
import sys
import urllib.request
from pathlib import Path

from paths import ANNOT, CACHE_ROOT, GEO, GTEX_DIR

FILES = [
    {
        "path": GEO / "GSE50244" / "GSE50244_series_matrix.txt.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50244/matrix/GSE50244_series_matrix.txt.gz",
        "bytes": 5992,
    },
    {
        "path": GEO / "GSE50244" / "GSE50244_Genes_counts_TMM_NormLength_atLeastMAF5_expressed.txt.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50244/suppl/GSE50244_Genes_counts_TMM_NormLength_atLeastMAF5_expressed.txt.gz",
        "bytes": 9188836,
    },
    {
        "path": GEO / "GSE38642" / "GSE38642_series_matrix.txt.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE38nnn/GSE38642/matrix/GSE38642_series_matrix.txt.gz",
        "bytes": 1872293,
    },
    {
        "path": GEO / "GSE76894" / "GSE76894_series_matrix.txt.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76894/matrix/GSE76894_series_matrix.txt.gz",
        "bytes": 16514968,
    },
    {
        "path": GEO / "GSE20966" / "GSE20966_series_matrix.txt.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE20nnn/GSE20966/matrix/GSE20966_series_matrix.txt.gz",
        "bytes": 3247995,
    },
    {
        "path": ANNOT / "GPL570.annot.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL5nnn/GPL570/annot/GPL570.annot.gz",
        "bytes": 8471521,
    },
    {
        "path": ANNOT / "GPL6244.annot.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6244/annot/GPL6244.annot.gz",
        "bytes": 7193437,
    },
    {
        "path": ANNOT / "GPL1352.annot.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL1nnn/GPL1352/annot/GPL1352.annot.gz",
        "bytes": 8765182,
    },
    {
        "path": GTEX_DIR / "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",
        "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",
        "bytes": 6952331,
    },
]

PROXY_KEYS = [
    "HTTP_PROXY",
    "HTTPS_PROXY",
    "ALL_PROXY",
    "http_proxy",
    "https_proxy",
    "all_proxy",
    "DOWNLOAD_PROXY",
]


def log(msg: str) -> None:
    print(msg, flush=True)


def present(item: dict) -> bool:
    path: Path = item["path"]
    return path.exists() and path.stat().st_size == item["bytes"]


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".partial")
    ctx = ssl.create_default_context()
    req = urllib.request.Request(url, headers={"User-Agent": "public-islet-t2d-identity-stress/1.0"})
    with urllib.request.urlopen(req, context=ctx, timeout=180) as resp, tmp.open("wb") as out:
        while True:
            chunk = resp.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
    tmp.replace(dest)


def download_one(item: dict) -> None:
    path: Path = item["path"]
    if present(item):
        log(f"OK  {path.name} ({item['bytes']} bytes)")
        return
    log(f"GET {item['url']}")
    try:
        fetch(item["url"], path)
    except Exception as exc:
        log(f"first attempt failed ({exc}); retry without proxy")
        old = {k: os.environ.pop(k, None) for k in PROXY_KEYS}
        try:
            fetch(item["url"], path)
        finally:
            for k, v in old.items():
                if v is not None:
                    os.environ[k] = v
    if not present(item):
        size = path.stat().st_size if path.exists() else 0
        raise RuntimeError(f"{path} size {size} != expected {item['bytes']}")
    log(f"OK  {path.name}")


def main() -> None:
    CACHE_ROOT.mkdir(parents=True, exist_ok=True)
    log(f"cache: {CACHE_ROOT}")
    failed = []
    for item in FILES:
        try:
            download_one(item)
        except Exception as exc:
            failed.append((item["path"].name, str(exc)))
            log(f"FAIL {item['path'].name}: {exc}")
    if failed:
        log("Download incomplete:")
        for name, err in failed:
            log(f"  {name}: {err}")
        sys.exit(1)
    log("All public matrices present.")


if __name__ == "__main__":
    main()
