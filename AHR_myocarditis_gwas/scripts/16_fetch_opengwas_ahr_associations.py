#!/usr/bin/env python3
"""Fetch selected variant associations from OpenGWAS.

The OpenGWAS JWT is read from OPENGWAS_JWT. If the environment variable is not
set, the script reads one line from stdin. The token is never written to disk.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

import pandas as pd


DEFAULT_ID = "eqtl-a-ENSG00000106546"
DEFAULT_VARIANTS = ["rs4843270", "rs61825638", "rs3184504", "rs6540080", "rs10216901"]


def project_root() -> Path:
    return Path(__file__).resolve().parents[1]


def get_token() -> str:
    token = os.environ.get("OPENGWAS_JWT", "").strip()
    if token:
        return token
    token = sys.stdin.readline().strip()
    if not token:
        raise SystemExit("OPENGWAS_JWT is not set and no token was provided on stdin.")
    return token


def post_json(url: str, payload: dict[str, Any], token: str) -> Any:
    data = json.dumps(payload).encode("utf-8")
    request = urllib.request.Request(
        url,
        data=data,
        headers={
            "Authorization": f"Bearer {token}",
            "Content-Type": "application/json",
            "Accept": "application/json",
        },
        method="POST",
    )
    try:
        with urllib.request.urlopen(request, timeout=180) as response:
            text = response.read().decode("utf-8")
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")
        raise SystemExit(f"OpenGWAS HTTP {exc.code}: {body}") from exc
    return json.loads(text)


def normalise_response(data: Any) -> pd.DataFrame:
    if isinstance(data, dict):
        for key in ("data", "results", "associations"):
            if key in data:
                data = data[key]
                break
    if data in (None, [], {}):
        return pd.DataFrame()
    return pd.json_normalize(data)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--id", default=DEFAULT_ID, help="OpenGWAS dataset id.")
    parser.add_argument(
        "--variants",
        nargs="+",
        default=DEFAULT_VARIANTS,
        help="rsIDs or chr:pos variants to query.",
    )
    parser.add_argument(
        "--out",
        default=str(
            project_root()
            / "data/raw/exposure/opengwas/eqtl-a-ENSG00000106546.associations.tsv.gz"
        ),
        help="Output gzip TSV path.",
    )
    parser.add_argument("--proxies", type=int, default=0, choices=[0, 1])
    args = parser.parse_args()

    token = get_token()
    payload = {
        "variant": args.variants,
        "id": [args.id],
        "proxies": args.proxies,
        "r2": 0.8,
        "align_alleles": 1,
        "palindromes": 1,
        "maf_threshold": 0.3,
    }

    data = post_json("https://api.opengwas.io/api/associations", payload, token)
    df = normalise_response(data)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out, "wt", encoding="utf-8", newline="") as handle:
        df.to_csv(handle, sep="\t", index=False)

    missing = sorted(set(args.variants) - set(df.get("rsid", pd.Series(dtype=str)).astype(str)))
    status = {
        "out": str(out),
        "dataset": args.id,
        "requested_variants": args.variants,
        "n_rows": int(len(df)),
        "columns": list(df.columns),
        "missing_by_rsid_column": missing,
    }
    print(json.dumps(status, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
