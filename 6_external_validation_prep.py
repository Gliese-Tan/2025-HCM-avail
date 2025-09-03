import scanpy as sc
import pandas as pd
import numpy as np
from scipy.stats import zscore

adata = sc.read_h5ad("2/2hcm_dcm.h5ad")
meta = pd.read_csv("2/clinicalmeta.csv")

meta.columns = meta.columns.str.strip().str.lower().str.replace(" ", "_")
meta = meta.rename(columns={"individual": "patient_id"})
meta["patient_id"] = meta["patient_id"].astype(str).str.strip()
adata.obs["donor_id"] = adata.obs["donor_id"].astype(str)
adata.obs["patient_id"] = adata.obs["donor_id"].str.replace("P", "", regex=False)
adata.obs = adata.obs.merge(meta, on="patient_id", how="left")
if "disease_y" in adata.obs.columns:
    adata.obs["disease"] = adata.obs["disease_y"]
    adata.obs = adata.obs.drop(columns=[c for c in ["disease_x", "disease_y"] if c in adata.obs.columns])

def get_gene_expr(a, gene):
    names = pd.Index(a.var_names)
    m = names.str.upper() == gene.upper()
    if not m.any():
        raise ValueError(f"{gene} not found")
    g = names[m][0]
    vals = a[:, g].X
    if hasattr(vals, "toarray"):
        vals = vals.toarray()
    return np.asarray(vals).ravel()

cm_adata = adata[adata.obs["cell_type_leiden0.6"].astype(str).str.startswith("Cardiomyocyte")].copy()
cm_adata.obs["HBEGF"] = get_gene_expr(cm_adata, "HBEGF")
hb_expr = (
    cm_adata.obs[["donor_id", "HBEGF"]]
    .groupby("donor_id").mean()
    .rename(columns={"HBEGF": "HBEGF_expr"})
)

fb_adata = adata[adata.obs["cell_type_leiden0.6"].astype(str).str.startswith("Fibroblast")].copy()
fb_adata.obs["THBS1"] = get_gene_expr(fb_adata, "THBS1")
thbs_expr = (
    fb_adata.obs[["donor_id", "THBS1"]]
    .groupby("donor_id").mean()
    .rename(columns={"THBS1": "THBS1_expr"})
)

clinical_vars = [
    "disease",
    "left_ventricular_ejection_fraction",
    "heart_weight",
    "left_ventricularmass",
    "left_ventricular_end-diastolic_dimension",
    "posterior_wall_thickness",
]
clinical_df = (
    adata.obs[["donor_id"] + clinical_vars]
    .drop_duplicates("donor_id")
    .set_index("donor_id")
)

merged = pd.concat([hb_expr, thbs_expr, clinical_df], axis=1).dropna()

records = []
for gene in ["HBEGF_expr", "THBS1_expr"]:
    for var in clinical_vars[1:]:
        for donor_id, row in merged.iterrows():
            records.append({
                "Gene": gene.replace("_expr", ""),
                "Clinical_Variable": var,
                "Expression": row[gene],
                "Clinical_Value": row[var],
                "Disease": row["disease"]
            })
plot_df = pd.DataFrame(records)
plot_df.to_csv("20250507/scatter_plot_data.csv", index=False)

hb_only_vars = [
    "left_ventricular_ejection_fraction",
    "heart_weight",
    "left_ventricularmass",
    "left_ventricular_end-diastolic_dimension",
    "posterior_wall_thickness",
]
cm_clinical = (
    cm_adata.obs[["donor_id"] + hb_only_vars + ["disease"]]
    .drop_duplicates("donor_id")
    .set_index("donor_id")
)
hb_plot = pd.concat([hb_expr, cm_clinical], axis=1).dropna()
zscored = hb_plot.drop(columns="disease").apply(zscore)
zscored = zscored.sort_values("HBEGF_expr", ascending=True)
disease_ordered = hb_plot.loc[zscored.index, "disease"]
zscored.to_csv("20250507/HBEGF_clustermap_zscore_matrix.csv")
disease_ordered.to_csv("20250507/HBEGF_clustermap_disease_labels.csv")

median_expr = merged["HBEGF_expr"].median()
grp = pd.Categorical(["Low" if x < median_expr else "High" for x in merged["HBEGF_expr"]], categories=["Low", "High"], ordered=True)
boxplot_df = pd.DataFrame({
    "donor_id": merged.index,
    "HBEGF_group": grp,
    "disease": merged["disease"],
    "left_ventricular_ejection_fraction": merged["left_ventricular_ejection_fraction"],
})
boxplot_df = boxplot_df.rename(columns={"donor_id": "Patient", "disease": "Disease", "left_ventricular_ejection_fraction": "LVEF"})
boxplot_df.to_csv("20250507/boxplot_LVEF_by_HBEGFgroup.csv", index=False)
