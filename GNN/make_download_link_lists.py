#!/usr/bin/env python3
import os
import sys
import time
from collections import deque
from html.parser import HTMLParser
from urllib.parse import urljoin, urldefrag, urlparse
from urllib.request import Request, urlopen

OUT_DIR = "/home/data1/wenhuai/GNN/download_links"

DATASETS = {
    "GSE167029": [
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167029/suppl/",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167029/matrix/",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167029/soft/",
    ],
    "GSE166489": [
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166489/suppl/",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166489/matrix/",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166489/soft/",
    ],
    "GSE183716": [
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183716/suppl/",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183716/matrix/",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183716/soft/",
    ],
    "GSE180045": [
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180045/suppl/",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180045/matrix/",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180045/soft/",
    ],
    "CNP0005824": [
        "https://ftp.cngb.org/pub/CNSA/data5/CNP0005824/",
    ],
    "STT0000127": [
        "https://ftp.cngb.org/pub/stomics/STT0000127/",
    ],
}

SKIP_EXT = (
    ".html", ".htm", ".css", ".js", ".png", ".gif", ".svg",
    ".ico", ".jpg", ".jpeg", ".webp"
)

SKIP_PARTS = (
    "/assets/",
    "/static/",
    "javascript:",
    "mailto:",
)

class LinkParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.hrefs = []

    def handle_starttag(self, tag, attrs):
        if tag.lower() != "a":
            return
        for k, v in attrs:
            if k.lower() == "href" and v:
                self.hrefs.append(v)

def fetch_html(url):
    req = Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urlopen(req, timeout=60) as r:
        data = r.read()
    return data.decode("utf-8", errors="ignore")

def is_dir_url(url):
    return urlparse(url).path.endswith("/")

def should_skip(url):
    low = url.lower()
    if any(x in low for x in SKIP_PARTS):
        return True
    path = urlparse(url).path.lower()
    if path.endswith(SKIP_EXT):
        return True
    return False

def crawl_base(base_url):
    base_url = base_url.rstrip("/") + "/"
    q = deque([base_url])
    seen_pages = set()
    files = set()

    while q:
        page = q.popleft()
        if page in seen_pages:
            continue
        seen_pages.add(page)

        try:
            html = fetch_html(page)
        except Exception as e:
            print(f"[WARN] cannot open: {page} | {e}", file=sys.stderr)
            continue

        parser = LinkParser()
        parser.feed(html)

        for href in parser.hrefs:
            if not href:
                continue

            u = urljoin(page, href)
            u, _frag = urldefrag(u)

            if not u.startswith(base_url):
                continue
            if should_skip(u):
                continue
            if u == page:
                continue

            if is_dir_url(u):
                if u not in seen_pages:
                    q.append(u)
            else:
                files.add(u)

        time.sleep(0.05)

    return sorted(files)

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    all_links = []

    for name, bases in DATASETS.items():
        print(f"\n========== {name} ==========")
        dataset_links = []

        for base in bases:
            print(f"[CRAWL] {base}")
            links = crawl_base(base)
            print(f"[FOUND] {len(links)} files")
            dataset_links.extend(links)

        dataset_links = sorted(set(dataset_links))
        all_links.extend(dataset_links)

        out_path = os.path.join(OUT_DIR, f"{name}_links.txt")
        with open(out_path, "w", encoding="utf-8") as f:
            for u in dataset_links:
                f.write(u + "\n")

        print(f"[WRITE] {out_path} | {len(dataset_links)} links")

    all_links = sorted(set(all_links))
    all_path = os.path.join(OUT_DIR, "ALL_links.txt")
    with open(all_path, "w", encoding="utf-8") as f:
        for u in all_links:
            f.write(u + "\n")

    print("\n========== SUMMARY ==========")
    print(f"[WRITE] {all_path} | {len(all_links)} links")
    print("[DONE] Link lists generated.")

if __name__ == "__main__":
    main()
