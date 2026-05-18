#!/usr/bin/env python3
"""Download OpenGWAS dataset files returned by /gwasinfo/files."""

from __future__ import annotations

import argparse
import json
import os
import sys
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any


DEFAULT_ID = "eqtl-a-ENSG00000106546"


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
    request = urllib.request.Request(
        url,
        data=json.dumps(payload).encode("utf-8"),
        headers={
            "Authorization": f"Bearer {token}",
            "Content-Type": "application/json",
            "Accept": "application/json",
        },
        method="POST",
    )
    try:
        with urllib.request.urlopen(request, timeout=180) as response:
            return json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")
        raise SystemExit(f"OpenGWAS HTTP {exc.code}: {body}") from exc


def download(url: str, outdir: Path) -> dict[str, Any]:
    filename = url.split("/")[-1].split("?")[0]
    outfile = outdir / filename
    with urllib.request.urlopen(url, timeout=300) as response:
        content = response.read()
    outfile.write_bytes(content)
    return {"file": str(outfile), "bytes": len(content)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--id", default=DEFAULT_ID)
    parser.add_argument(
        "--outdir",
        default=str(project_root() / "data/raw/exposure/opengwas"),
    )
    args = parser.parse_args()

    token = get_token()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    data = post_json("https://api.opengwas.io/api/gwasinfo/files", {"id": [args.id]}, token)
    urls = data.get(args.id, [])
    if not urls:
        raise SystemExit(f"No files returned for {args.id}")

    downloaded = [download(url, outdir) for url in urls]
    print(json.dumps({"dataset": args.id, "downloaded": downloaded}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
