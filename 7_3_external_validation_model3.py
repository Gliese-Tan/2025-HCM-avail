import scanpy as sc
import pandas as pd
import numpy as np
import joblib

adata = sc.read_h5ad('2/2hcm_dcm.h5ad')

target_types = [
    'Cardiomyocyte_I', 'Cardiomyocyte_II', 'Cardiomyocyte_III',
    'Fibroblast_I', 'Fibroblast_II', 'Activated_fibroblast',
    'Macrophage', 'Proliferating_macrophage'
]
adata_filtered = adata[adata.obs['cell_type_leiden0.6'].isin(target_types)].copy()

sc.pp.normalize_total(adata_filtered, target_sum=1e4)
sc.pp.log1p(adata_filtered)
sc.pp.highly_variable_genes(adata_filtered, n_top_genes=2000, subset=True)
sc.pp.scale(adata_filtered)
sc.tl.pca(adata_filtered, n_comps=30)

model = joblib.load('va_nn_model3.pkl')
X_external = adata_filtered.obsm['X_pca'][:, :30]
classes = list(model.classes_)
va_idx = classes.index('VA') if 'VA' in classes else None
va_prob = model.predict_proba(X_external)[:, va_idx] if va_idx is not None else np.full(X_external.shape[0], np.nan)

adata_filtered.obs['donor_id'] = adata_filtered.obs['donor_id'].astype(str)
adata_filtered.obs['VA_score'] = va_prob

donor_scores = (
    adata_filtered.obs.groupby('donor_id')['VA_score']
    .mean()
    .reset_index()
    .rename(columns={'VA_score': 'avg_VA_score'})
)
donor_scores['donor_id'] = donor_scores['donor_id'].str.replace('^P', '', regex=True)

clinical = pd.read_csv('2/clinicalmeta.csv')
clinical.columns = clinical.columns.str.lower().str.strip().str.replace(' ', '_')
clinical['donor_id'] = clinical['individual'].astype(str).str.strip()

merged = donor_scores.merge(
    clinical[['donor_id', 'left_ventricular_ejection_fraction']],
    on='donor_id', how='inner'
).rename(columns={'left_ventricular_ejection_fraction': 'lvef'})

merged.dropna().to_csv('VA_score_LVEF_model3_table.csv', index=False)
