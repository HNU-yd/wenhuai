"""
Kynurenine-AHR V0 gene sets.

V0-beta 设计原则：
1. 模块拆开看，不只依赖 composite score。
2. human/mouse 用同一套 canonical symbols，通过大小写不敏感匹配适配。
3. AHR_regulon_proxy_score 只是 proxy，不等于正式 SCENIC regulon。
"""

DATASET_SPECIES = {
    "CNP0005824": "mouse",
    "GSE166489": "human",
    "GSE167029": "human",
    "GSE183716": "human",
    "GSE180045": "human",
    "STT0000127": "mouse",
}

GENE_SETS_HUMAN = {
    "KYN_metabolism_score": [
        "IDO1", "IDO2", "TDO2", "AFMID", "KMO", "KYNU", "HAAO", "QPRT", "AADAT"
    ],
    "AHR_response_score": [
        "AHR", "AHRR", "ARNT", "CYP1A1", "CYP1B1",
        "NQO1", "TIPARP", "PTGS2", "ALDH1A1"
    ],
    "AHR_regulon_proxy_score": [
        "AHR", "AHRR", "ARNT", "CYP1A1", "CYP1B1",
        "NQO1", "TIPARP", "PTGS2", "ALDH1A1",
        "IL1B", "NFKBIA", "STAT1", "IRF1"
    ],
    "myeloid_inflammation_score": [
        "IL1B", "TNF", "NFKBIA", "IRF1", "STAT1",
        "S100A8", "S100A9", "CCL2", "CCL3", "CCL4", "CXCL8"
    ],
    "chemotaxis_score": [
        "CCR1", "CCR2", "CCR5", "CXCR4",
        "CXCL8", "CCL2", "CCL3", "CCL4", "CCL5",
        "ITGAM", "ITGB2"
    ],
    "antigen_presentation_score": [
        "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1",
        "HLA-DQA1", "CD74", "CIITA"
    ],
}

# 这里保留 human canonical 名称，实际匹配时会大小写不敏感；
# mouse 数据里常见 Ido1 / Ahr / S100a8 会被 present_genes() 正确匹配到。
GENE_SETS_MOUSE = GENE_SETS_HUMAN

COMPOSITE_SCORE_CONFIGS = {
    # 原 composite：KYN + AHR response + AHR regulon proxy + myeloid inflammation
    "Kyn_AHR_myeloid_score": [
        "KYN_metabolism_score",
        "AHR_response_score",
        "AHR_regulon_proxy_score",
        "myeloid_inflammation_score",
    ],

    # 新增：不含 broad inflammation 的版本，避免广义炎症主导全部结果。
    "Kyn_AHR_axis_score": [
        "KYN_metabolism_score",
        "AHR_response_score",
        "AHR_regulon_proxy_score",
    ],

    # 新增：myeloid effector 状态，不混 KYN/AHR。
    "myeloid_effector_score": [
        "myeloid_inflammation_score",
        "chemotaxis_score",
        "antigen_presentation_score",
    ],
}

CELLTYPE_MARKERS = {
    "T_cell": ["CD3D", "CD3E", "TRAC", "IL7R", "CCR7"],
    "CD8_T": ["CD8A", "CD8B", "GZMK", "GZMH"],
    "NK": ["NKG7", "GNLY", "KLRD1", "PRF1"],
    "B_cell": ["MS4A1", "CD79A", "CD79B", "BANK1"],
    "Plasmablast": ["MZB1", "XBP1", "JCHAIN", "IGHG1"],
    "Monocyte": ["LYZ", "LST1", "S100A8", "S100A9", "FCN1"],
    "CD16_Monocyte": ["FCGR3A", "MS4A7", "LST1", "CTSS"],
    "DC": ["FCER1A", "CLEC10A", "CST3", "CD1C"],
    "pDC": ["LILRA4", "GZMB", "IRF7", "TCF4"],
    "Platelet": ["PPBP", "PF4", "NRGN"],
}

MYELOID_CELLTYPES = {
    "Monocyte",
    "CD16_Monocyte",
    "DC",
    "pDC",
}

MODULE_SCORE_COLUMNS = [
    "KYN_metabolism_score",
    "AHR_response_score",
    "AHR_regulon_proxy_score",
    "myeloid_inflammation_score",
    "chemotaxis_score",
    "antigen_presentation_score",
]

COMPOSITE_SCORE_COLUMNS = list(COMPOSITE_SCORE_CONFIGS.keys())

ALL_SCORE_COLUMNS = MODULE_SCORE_COLUMNS + COMPOSITE_SCORE_COLUMNS


def get_gene_sets_for_dataset(dataset: str):
    species = DATASET_SPECIES.get(dataset, "human")
    if species == "mouse":
        return GENE_SETS_MOUSE
    return GENE_SETS_HUMAN
