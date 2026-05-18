import re
from pathlib import Path
from typing import Dict, Optional

import pandas as pd


def safe_name(x: str) -> str:
    x = str(x)
    x = re.sub(r"[^\w.\-]+", "_", x)
    x = re.sub(r"_+", "_", x)
    return x.strip("_")


def extract_gsm(x: str) -> Optional[str]:
    m = re.search(r"(GSM\d+)", str(x))
    return m.group(1) if m else None


def clean_sample_id_from_path(path: str) -> str:
    p = Path(path)
    stem = p.name

    for suffix in [
        "_filtered_feature_bc_matrix",
        "-filtered_feature_bc_matrix",
        "_raw_feature_bc_matrix",
        "-raw_feature_bc_matrix",
        "_feature_bc_matrix",
        "-feature_bc_matrix",
    ]:
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]

    for suf in [
        ".h5ad", ".h5", ".hdf5", ".rds", ".RDS",
        ".mtx", ".mtx.gz", ".txt", ".txt.gz",
        ".tar", ".tar.gz", ".tgz", ".zip",
    ]:
        stem = stem.replace(suf, "")

    return safe_name(stem)


def parse_geo_soft_file(path: Path, dataset: str):
    rows = []
    if not path.exists():
        return rows

    current = {}
    with path.open("rt", errors="ignore") as f:
        for line in f:
            line = line.rstrip("\n")

            if line.startswith("^SAMPLE"):
                if current.get("gsm"):
                    rows.append(current)
                current = {
                    "dataset": dataset,
                    "gsm": "",
                    "title": "",
                    "source_name": "",
                    "characteristics": [],
                }

            elif line.startswith("!Sample_geo_accession"):
                current["gsm"] = line.split("=", 1)[-1].strip()

            elif line.startswith("!Sample_title"):
                current["title"] = line.split("=", 1)[-1].strip()

            elif line.startswith("!Sample_source_name_ch1"):
                current["source_name"] = line.split("=", 1)[-1].strip()

            elif line.startswith("!Sample_characteristics_ch1"):
                current.setdefault("characteristics", []).append(
                    line.split("=", 1)[-1].strip()
                )

    if current.get("gsm"):
        rows.append(current)

    for r in rows:
        r["characteristics"] = " | ".join(r.get("characteristics", []))

    return rows


def default_data_root() -> Path:
    return Path(__file__).resolve().parents[2] / "data"


def build_geo_sample_table(data_root: Optional[str] = None) -> pd.DataFrame:
    data_root = Path(data_root) if data_root is not None else default_data_root()
    all_rows = []

    for ds_root in sorted(data_root.iterdir()):
        if not ds_root.is_dir():
            continue

        dataset = ds_root.name
        for soft in ds_root.rglob("*family.soft"):
            all_rows.extend(parse_geo_soft_file(soft, dataset))

    if not all_rows:
        return pd.DataFrame(
            columns=["dataset", "gsm", "title", "source_name", "characteristics"]
        )

    return pd.DataFrame(all_rows).drop_duplicates()


def load_geo_sample_table(path: str) -> pd.DataFrame:
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return pd.DataFrame(
            columns=["dataset", "gsm", "title", "source_name", "characteristics"]
        )
    return pd.read_csv(p, sep="\t").fillna("")


def sample_context_from_geo(
    dataset: str,
    path: str,
    geo_df: Optional[pd.DataFrame] = None,
) -> Dict[str, str]:
    gsm = extract_gsm(path)
    out = {
        "gsm": gsm or "",
        "geo_title": "",
        "geo_source_name": "",
        "geo_characteristics": "",
    }

    if geo_df is None or geo_df.empty or not gsm:
        return out

    hit = geo_df[(geo_df["dataset"] == dataset) & (geo_df["gsm"] == gsm)]
    if len(hit) == 0:
        return out

    r = hit.iloc[0]
    out["geo_title"] = str(r.get("title", ""))
    out["geo_source_name"] = str(r.get("source_name", ""))
    out["geo_characteristics"] = str(r.get("characteristics", ""))

    return out


def infer_gse166489_patient_and_stage(sample_id: str, context: str):
    """
    GSE166489:
    P3.1 / P3.2 -> patient_id P3
    P4.1 / P4.2 -> patient_id P4
    .1 通常 acute，.2 通常 recovery
    """
    patient_id = sample_id
    stage = "unknown"

    m = re.search(r"(P\d+)\.(\d+)", context)
    if m:
        patient_id = m.group(1)
        suffix = m.group(2)

        if suffix == "1":
            stage = "acute"
        elif suffix == "2":
            stage = "recovery"

    m_hd = re.search(r"(C\.HD\d+|A\.HD\d+)", context)
    if m_hd:
        patient_id = m_hd.group(1)
        stage = "control"

    if re.search(r"recover|convalescent", context, re.I):
        stage = "recovery"

    return safe_name(patient_id), stage


def infer_group_severity_stage(dataset: str, sample_id: str, path: str, geo_context: Dict[str, str]):
    s = " ".join([
        dataset,
        sample_id,
        path,
        geo_context.get("geo_title", ""),
        geo_context.get("geo_source_name", ""),
        geo_context.get("geo_characteristics", ""),
    ])

    group = "unknown"
    severity = "unknown"
    stage = "unknown"

    if dataset == "CNP0005824":
        if re.search(r"\bControl\b", s, re.I):
            group, severity, stage = "control", "control", "control"
        elif re.search(r"\bCVB3", s, re.I):
            group, severity, stage = "CVB3_myocarditis", "disease", "acute"
        elif re.search(r"\bIVIG", s, re.I):
            group, severity, stage = "IVIG_treated", "treated", "post_treatment"

    elif dataset == "GSE183716":
        if re.search(r"Sample1", s, re.I):
            group, severity, stage = "MIS-C", "mis_c", "acute"
        elif re.search(r"Sample2", s, re.I):
            group, severity, stage = "MIS-C", "mis_c", "acute"
        elif re.search(r"Sample3", s, re.I):
            group, severity, stage = "MIS-C", "mis_c", "recovery"
        elif re.search(r"Sample4", s, re.I):
            group, severity, stage = "MIS-C_or_KD", "unknown", "unknown"

    elif dataset == "GSE166489":
        if re.search(r"C\.HD|pediatric healthy|child healthy|healthy child", s, re.I):
            group, severity, stage = "pediatric_healthy", "control", "control"
        elif re.search(r"A\.HD|adult healthy|healthy adult", s, re.I):
            group, severity, stage = "adult_healthy", "adult_control", "control"
        elif re.search(r"MIS[-_ ]?C|P\d+\.\d+", s, re.I):
            group, severity = "MIS-C", "unknown"
            _, inferred_stage = infer_gse166489_patient_and_stage(sample_id, s)
            stage = inferred_stage

        if re.search(r"severe", s, re.I):
            severity = "severe"
        elif re.search(r"moderate", s, re.I):
            severity = "moderate"

        if re.search(r"recover|convalescent", s, re.I):
            stage = "recovery"
            if severity == "unknown":
                severity = "recovered"

    elif dataset == "GSE167029":
        if re.search(r"\bC\d+|control|healthy|CTL", s, re.I):
            group, severity, stage = "control", "control", "control"

        if re.search(r"MIS[-_ ]?C[_ -]?MYO|MISC[_ -]?MYO|myocarditis|MYO", s, re.I):
            group, severity, stage = "MIS-C_MYO", "severe_myo", "acute"
        elif re.search(r"MIS[-_ ]?C|MISC", s, re.I):
            group, severity, stage = "MIS-C", "mis_c", "acute"
        elif re.search(r"acute infection|acute[-_ ]?inf|infection", s, re.I):
            group, severity, stage = "acute_infection", "non_mis_c", "acute"
        elif re.search(r"\bKD\b|Kawasaki", s, re.I):
            group, severity, stage = "KD", "disease_control", "acute"

        if group == "unknown":
            if re.search(r"GSM\d+_C\d+", s, re.I):
                group, severity, stage = "control", "control", "control"
            elif re.search(r"GSM\d+_P\d+", s, re.I):
                group, severity, stage = "patient", "unknown", "acute"

    elif dataset == "GSE180045":
        group, severity, stage = "ICI_myocarditis_extension", "extension", "unknown"
        if re.search(r"PBMC|pbmc", s):
            stage = "pbmc"
        if re.search(r"MCE|MCL", s):
            group = "ICI_myocarditis"
        if re.search(r"\bA\d+|\bB\d+", s):
            group = "control_or_reference"

    return group, severity, stage


def infer_meta_from_path(
    dataset: str,
    path: str,
    geo_df: Optional[pd.DataFrame] = None,
) -> Dict[str, str]:
    sample_id = clean_sample_id_from_path(path)

    gsm = extract_gsm(path)
    if gsm and not sample_id.startswith(gsm):
        sample_id = safe_name(f"{gsm}_{sample_id}")

    geo_context = sample_context_from_geo(dataset, path, geo_df)
    group, severity, stage = infer_group_severity_stage(
        dataset=dataset,
        sample_id=sample_id,
        path=path,
        geo_context=geo_context,
    )

    patient_id = sample_id
    context = " ".join([
        sample_id,
        path,
        geo_context.get("geo_title", ""),
        geo_context.get("geo_source_name", ""),
        geo_context.get("geo_characteristics", ""),
    ])

    if dataset == "GSE166489":
        patient_id, inferred_stage = infer_gse166489_patient_and_stage(sample_id, context)
        if stage == "unknown":
            stage = inferred_stage

    elif dataset == "GSE167029":
        m = re.search(r"_([CP]\d+(?:_bis)?)", sample_id)
        if m:
            patient_id = safe_name(m.group(1))

    elif dataset == "CNP0005824":
        patient_id = sample_id

    elif dataset == "GSE183716":
        if "Sample1" in sample_id:
            patient_id = "TRICOV19"
        elif "Sample2" in sample_id:
            patient_id = "TRICOV20"
        elif "Sample3" in sample_id:
            patient_id = "TRICOV20"
        elif "Sample4" in sample_id:
            patient_id = "TRICOV22"

    return {
        "dataset": dataset,
        "source_file": path,
        "gsm": geo_context.get("gsm", "") or gsm or "",
        "sample_id": safe_name(sample_id),
        "patient_id": safe_name(patient_id),
        "group": group,
        "severity": severity,
        "stage": stage,
        "batch": dataset,
        "geo_title": geo_context.get("geo_title", ""),
        "geo_source_name": geo_context.get("geo_source_name", ""),
        "geo_characteristics": geo_context.get("geo_characteristics", ""),
    }
