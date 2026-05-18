#!/usr/bin/env python3
from pathlib import Path
import pandas as pd

rows = []

def add(dataset, gsm, sample_id, patient_id, group_v0, severity, stage, include_v0=True, note=""):
    rows.append(dict(
        dataset=dataset,
        gsm=gsm,
        sample_id=sample_id,
        patient_id=patient_id,
        group_v0=group_v0,
        severity=severity,
        stage=stage,
        include_v0=int(include_v0),
        note=note,
    ))

# GSE167029: children hospitalized with suspicion of SARS-CoV-2 infection;
# main V0 contrast: CTL vs MIS-C vs MIS-C_MYO.
for gsm, sid in [
    ("GSM5091035", "C24"), ("GSM5091036", "C26"), ("GSM5091037", "C27"),
    ("GSM5091038", "C27_bis"), ("GSM5091039", "C28"), ("GSM5091040", "C29"),
    ("GSM5091041", "C30"), ("GSM5091042", "C30_bis"), ("GSM5091043", "C31"),
    ("GSM5091044", "C32"), ("GSM5091045", "C39"),
]:
    add("GSE167029", gsm, sid, sid.replace("_bis", ""), "CTL", "control", "control")

add("GSE167029", "GSM5091046", "P3", "P3", "AcuteInf_CoV2neg", "acute_infection", "acute")
for gsm, sid in [("GSM5091047", "P5"), ("GSM5091048", "P8"), ("GSM5091049", "P12"), ("GSM5091050", "P13")]:
    add("GSE167029", gsm, sid, sid, "AcuteInf_CoV2pos", "acute_infection", "acute")
for gsm, sid in [("GSM5091051", "P17"), ("GSM5091052", "P20")]:
    add("GSE167029", gsm, sid, sid, "MIS-C", "mild_or_nonmyocarditis", "acute")
for gsm, sid in [("GSM5091053", "P23"), ("GSM5091054", "P31"), ("GSM5091055", "P33"), ("GSM5091056", "P37"), ("GSM5091057", "P39"), ("GSM5091058", "P43")]:
    add("GSE167029", gsm, sid, sid, "MIS-C_MYO", "severe_myocarditis", "acute")
for gsm, sid in [("GSM5091059", "P47"), ("GSM5091060", "P53"), ("GSM5091061", "P55")]:
    add("GSE167029", gsm, sid, sid, "KD", "kawasaki_reference", "acute", True, "Reference inflammatory pediatric disease")

# GSE166489: use GEX samples first; ADT/TCR/BCR are ignored in V0.0.
for gsm, sid, sev in [
    ("GSM5073055", "P1.1", "severe"), ("GSM5073056", "P2.1", "severe"),
    ("GSM5073057", "P3.1", "severe"), ("GSM5073058", "P4.1", "moderate"),
    ("GSM5073059", "P5.1", "moderate"), ("GSM5073060", "P6.1", "severe"),
    ("GSM5073061", "P7.1", "severe"),
]:
    patient = sid.split(".")[0]
    add("GSE166489", gsm, sid, patient, f"MIS-C_{sev}", sev, "acute")
for gsm, sid, sev in [("GSM5073062", "P3.2", "severe"), ("GSM5073063", "P4.2", "moderate")]:
    patient = sid.split(".")[0]
    add("GSE166489", gsm, sid, patient, f"MIS-C_{sev}_recovered", sev, "recovered")
for n, gsm in enumerate(["GSM5073064", "GSM5073065", "GSM5073066", "GSM5073067", "GSM5073068", "GSM5073069"], start=1):
    add("GSE166489", gsm, f"C.HD{n}", f"C.HD{n}", "pediatric_healthy", "control", "control")
for n, gsm in enumerate(["GSM5073070", "GSM5073071", "GSM5073072"], start=1):
    add("GSE166489", gsm, f"A.HD{n}", f"A.HD{n}", "adult_healthy", "adult_control", "control", False, "Excluded by default because V0 is pediatric-focused")

# GSE183716: minimal public sample labels from GEO title. Confirm clinical mapping before final manuscript.
add("GSE183716", "GSM5569718", "Sample1_TRICOV19_A", "TRICOV19", "MIS-C_acute", "acute", "acute")
add("GSE183716", "GSM5569719", "Sample2_TRICOV20_A", "TRICOV20", "MIS-C_acute", "acute", "acute")
add("GSE183716", "GSM5569720", "Sample3_TRICOV20_F", "TRICOV20", "MIS-C_followup", "recovery", "followup", True, "Follow-up/recovery sample by sample code F; verify with sample sheet before publication")
add("GSE183716", "GSM5569721", "Sample4_TRICOV22_K", "TRICOV22", "KD_reference", "kawasaki_reference", "acute", True, "Likely Kawasaki reference by code K; verify with sample sheet before publication")

out = Path(__file__).resolve().parents[1] / "metadata" / "metadata_public_v0.tsv"
out.parent.mkdir(parents=True, exist_ok=True)
pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
print(out)
