#!/usr/bin/env python
import argparse
import gc
import os
from pathlib import Path

import anndata as ad
import scanpy as sc
from scipy import sparse
from tqdm import tqdm

from project_paths import project_path


def force_csr(adata):
    """
    anndata.experimental.concat_on_disk 对 csc 支持不好。
    fallback memory concat 也统一转 CSR，避免后续写盘/合并不稳定。
    """
    if sparse.issparse(adata.X) and not sparse.isspmatrix_csr(adata.X):
        adata.X = adata.X.tocsr()

    for k in list(adata.layers.keys()):
        if sparse.issparse(adata.layers[k]) and not sparse.isspmatrix_csr(adata.layers[k]):
            adata.layers[k] = adata.layers[k].tocsr()

    return adata


def safe_unlink(path: Path):
    if path.exists():
        path.unlink()


def atomic_replace(tmp_file: Path, final_file: Path):
    final_file.parent.mkdir(parents=True, exist_ok=True)
    os.replace(str(tmp_file), str(final_file))


def try_concat_on_disk(files, out_file: Path, dataset: str):
    """
    只尝试 on-disk concat。
    这里故意写到 .on_disk_tmp.h5ad，不直接写正式 out_file。
    即使失败，也不会污染正式输出。
    """
    tmp_file = out_file.with_name(out_file.stem + ".on_disk_tmp.h5ad")
    safe_unlink(tmp_file)

    from anndata.experimental import concat_on_disk

    print(f"[try on-disk concat] {dataset} -> {tmp_file}")

    concat_on_disk(
        [str(x) for x in files],
        str(tmp_file),
        join="outer",
        label="sample_file",
        keys=[x.stem for x in files],
        index_unique=None,
        max_loaded_elems=50_000_000,
    )

    atomic_replace(tmp_file, out_file)
    print(f"[done on-disk] {out_file}")


def concat_memory_csr(files, out_file: Path, dataset: str):
    """
    稳定 fallback：
    1. 逐个 read_h5ad
    2. X/layers 转 CSR
    3. anndata.concat
    4. 写 memory_tmp
    5. 原子替换正式文件
    """
    tmp_file = out_file.with_name(out_file.stem + ".memory_tmp.h5ad")
    safe_unlink(tmp_file)

    print(f"[memory concat CSR] {dataset}: {len(files)} files")

    adatas = []
    for f in tqdm(files, desc=f"read {dataset}"):
        a = sc.read_h5ad(f)
        a = force_csr(a)

        # 保证 obs_names 唯一
        a.obs_names_make_unique()
        a.var_names_make_unique()

        adatas.append(a)

    print(f"[anndata.concat] {dataset}")
    merged = ad.concat(
        adatas,
        join="outer",
        label="sample_file",
        keys=[x.stem for x in files],
        index_unique=None,
        fill_value=0,
        merge="same",
        uns_merge="unique",
    )

    merged = force_csr(merged)
    merged.obs_names_make_unique()
    merged.var_names_make_unique()

    print(f"[write tmp] {tmp_file}")
    merged.write_h5ad(tmp_file, compression="gzip")

    atomic_replace(tmp_file, out_file)
    print(f"[done memory] {out_file}")

    del merged
    del adatas
    gc.collect()


def concat_one_dataset(dataset_dir: Path, out_file: Path, overwrite: bool = False, mode: str = "auto"):
    files = sorted(dataset_dir.glob("*.h5ad"))

    if not files:
        print(f"[skip] no h5ad in {dataset_dir}")
        return

    out_file.parent.mkdir(parents=True, exist_ok=True)

    if out_file.exists() and out_file.stat().st_size > 1024 and not overwrite:
        print(f"[skip] exists: {out_file}")
        return

    print("=" * 80)
    print(f"[dataset] {dataset_dir.name}")
    print(f"[n_files] {len(files)}")
    print(f"[out] {out_file}")

    # 不直接删除旧正式文件，避免失败后没有旧结果。
    # 成功写 tmp 后再 os.replace 替换。

    if mode in {"auto", "on_disk"}:
        try:
            try_concat_on_disk(files, out_file, dataset_dir.name)
            return
        except Exception as e:
            print(f"[warn] concat_on_disk failed for {dataset_dir.name}: {repr(e)}")
            tmp = out_file.with_name(out_file.stem + ".on_disk_tmp.h5ad")
            safe_unlink(tmp)

            if mode == "on_disk":
                raise

            print("[fallback] memory concat with CSR conversion")

    concat_memory_csr(files, out_file, dataset_dir.name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_root", default=str(project_path("data", "h5ad_per_sample")))
    ap.add_argument("--out_root", default=str(project_path("data", "h5ad_merged")))
    ap.add_argument("--datasets", nargs="*", default=[
        "GSE183716",
        "CNP0005824",
        "GSE166489",
        "GSE167029",
    ])
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument(
        "--mode",
        choices=["auto", "memory", "on_disk"],
        default="auto",
        help="auto: try concat_on_disk then fallback memory; memory: force memory CSR concat; on_disk: only on-disk concat",
    )
    args = ap.parse_args()

    for ds in args.datasets:
        concat_one_dataset(
            Path(args.in_root) / ds,
            Path(args.out_root) / f"{ds}.raw_merged.h5ad",
            overwrite=args.overwrite,
            mode=args.mode,
        )


if __name__ == "__main__":
    main()
