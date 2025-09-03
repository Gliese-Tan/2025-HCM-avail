import scvelo as scv
import scanpy as sc
import pandas as pd

samples = [
    "velocyto/H1_AL.loom","velocyto/H1_VA.loom","velocyto/H1_VS.loom",
    "velocyto/H2_AL.loom","velocyto/H2_VA.loom","velocyto/H2_VS.loom",
    "velocyto/H3_AL.loom","velocyto/H3_VA.loom","velocyto/H3_VS.loom",
    "velocyto/H4_AL.loom","velocyto/H4_VA.loom","velocyto/H4_VS.loom",
    "velocyto/H5_AL.loom","velocyto/H5_VA.loom","velocyto/H5_VS.loom",
    "velocyto/H6_AL.loom","velocyto/H6_VA.loom","velocyto/H6_VS.loom",
    "velocyto/H7_AL.loom","velocyto/H7_VA.loom","velocyto/H7_VS.loom",
    "velocyto/H8_AL.loom","velocyto/H8_VA.loom","velocyto/H8_VS.loom"
]

adata_list = [scv.read(p, cache=True) for p in samples]
for ad in adata_list:
    ad.var_names_make_unique()

adata = sc.concat(
    adata_list,
    label="sample",
    keys=[
        "H1_AL","H1_VA","H1_VS","H2_AL","H2_VA","H2_VS",
        "H3_AL","H3_VA","H3_VS","H4_AL","H4_VA","H4_VS",
        "H5_AL","H5_VA","H5_VS","H6_AL","H6_VA","H6_VS",
        "H7_AL","H7_VA","H7_VS","H8_AL","H8_VA","H8_VS"
    ]
)
adata = adata[:, ~adata.var_names.duplicated()]

fib_df = pd.read_csv("20241129/fibroblast_umap.csv")
fib_cells = fib_df["CellID"].tolist()
fibroblast_adata = adata[adata.obs_names.isin(fib_cells)].copy()

scv.pp.filter_genes(fibroblast_adata, min_cells=20)
scv.pp.filter_and_normalize(fibroblast_adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(fibroblast_adata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(fibroblast_adata)
scv.tl.velocity_graph(fibroblast_adata)

umap_df = fib_df.set_index("CellID").loc[fibroblast_adata.obs_names, ["UMAP1","UMAP2"]]
fibroblast_adata.obsm["X_umap"] = umap_df.values

anno = pd.read_csv("20241129/fibroblast_annotations.csv").set_index("fibroblast_cells")["cell_subtypes"]
fibroblast_adata.obs["cell_subtypes"] = fibroblast_adata.obs_names.map(anno)
fibroblast_adata = fibroblast_adata[~fibroblast_adata.obs["cell_subtypes"].isna()].copy()

fibroblast_adata.write("scvelo/1fibroblast_adata_with_velocity.h5ad")
