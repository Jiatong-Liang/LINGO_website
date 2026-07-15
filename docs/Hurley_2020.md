# Hurley

The [*Hurley 2020* dataset](installation.md) profiles cells from zebrafish hearts at multiple time points.

To see more details regarding the experiment and features of the dataset, please refer to: _Hurley, K., Ding, J., Villacorta-Martin, C., Herriges, M. J., Jacob, A., Vedaie, M., ... & Kotton, D. N. (2020). Reconstructed single-cell fate trajectories define lineage plasticity windows during differentiation of human PSC-derived distal lung progenitors. Cell stem cell, 26(4), 593-608._

The example will walk through the **clonal plasticity prediction,** which involves predicting the plasticity of a clone from its observed cells.

The manuscript discusses many variations of the LINGO model, but we recommend users to the *KG Trainable Embedding* instead of the *KG Scratch* and *KG Frozen*.

## Data Preprocessing

The first thing to do is to verify whether your dataset has been merged with or modified by other sources. Since our analysis focuses exclusively on Hu 2022, we subset the data to retain only samples from that study. In addition, we applied the following quality control to subset the data: nFeature_RNA >= 200, nCount_RNA >= 500, percent.mito <= 20. In the end we retained 178,634 out of 197,258 cells. 

```python
data_path = '/home/User/LINGO/data/Hurley_2020_NatureGenetics.h5ad'

print(f"  Reading: {data_path}")
with warnings.catch_warnings():
    warnings.filterwarnings(
        "ignore",
        message="Variable names are not unique.*",
        category=UserWarning,
    )
    adata = sc.read_h5ad(data_path)
adata.var_names_make_unique()

if "dataset" in adata.obs.columns:
    hu_mask = adata.obs["dataset"].astype(str).str.contains("Hu_2022", na=False)
    if hu_mask.sum() > 0:
        adata = adata[hu_mask].copy()
        print(f"  Filtered by dataset=Hu_2022: {adata.n_obs:,} cells")

print(f"  Loaded {adata.n_obs:,} cells x {adata.n_vars:,} genes")

qc_summary = {
    "enabled": False,
    "min_features": None,
    "min_counts": None,
    "max_percent_mito": None,
    "n_cells_before_qc": int(adata.n_obs),
    "n_cells_after_qc": int(adata.n_obs),
}
qc_masks = []
min_features = 200
min_counts = 500
max_percent_mito = 20

if "nFeature_RNA" in adata.obs.columns:
    feat_mask = np.asarray(adata.obs["nFeature_RNA"], dtype=float) >= float(min_features)
    qc_masks.append(feat_mask)
if "nCount_RNA" in adata.obs.columns:
    count_mask = np.asarray(adata.obs["nCount_RNA"], dtype=float) >= float(min_counts)
    qc_masks.append(count_mask)
if "percent.mito" in adata.obs.columns:
    mito_mask = np.asarray(adata.obs["percent.mito"], dtype=float) <= float(max_percent_mito)
    qc_masks.append(mito_mask)

if len(qc_masks) > 0:
    keep_mask = np.logical_and.reduce(qc_masks)
    n_before_qc = int(adata.n_obs)
    n_after_qc = int(keep_mask.sum())
    adata = adata[keep_mask].copy()
    qc_summary = {
        "enabled": True,
        "min_features": int(min_features),
        "min_counts": int(min_counts),
        "max_percent_mito": float(max_percent_mito),
        "n_cells_before_qc": n_before_qc,
        "n_cells_after_qc": n_after_qc,
    }
    print("  Applied Hu QC filter:")
    print(
        "    "
        f"nFeature_RNA >= {min_features}, "
        f"nCount_RNA >= {min_counts}, "
        f"percent.mito <= {max_percent_mito}"
    )
    print(f"    Retained {n_after_qc:,} / {n_before_qc:,} cells")
else:
    print("  QC columns not found, skipped cell-level Hu QC filter")
adata.uns["lineage_qc_filter_summary"] = qc_summary
```

[!Note] If your `adata` object does not contain the `dataset` column or if you are certain you are not working with a merged dataset, you may skip this step.

The second thing is to identify the column name that contains the barcodes, timepoints, and cell fate.

```python
from utilities_clonal_plasticity import _parse_timepoint_value

# The lineage task must use real barcode/clone fields and forgery is prohibited.
barcode_candidates = [
    'barcodes', 'barcode', 'clone', 'clone_id', 'CloneID', 'cell_barcode', 'lentiviral_barcode', 'lentiBC'
]

barcode_col = None
obs_cols = list(adata.obs.columns)
for c in barcode_candidates:
    if c in obs_cols:
        barcode_col = c
        break

time_candidates = [
    "timepoint", "time_point", "time point", "time", "timestamp",
    "time_info", "tp", "day", "dpt", "stage", "age", "Aging_Time"
]

time_col = None
for t in time_candidates:
    if t in obs_cols:
        time_col = t
        break

fate_candidates = [
    "celltypes", "celltype", "subclustered.celltypes",
]

fate_col = None
for f in fate_candidates:
    if f in obs_cols:
        fate_col = f
        break

print(f"  Using barcode column: {barcode_col}")
print(f"  Using time column: {time_col}")
print(f"  Using fate column: {fate_col}")

vals = adata.obs[time_col].astype(str).values
parsed = [_parse_timepoint_value(v) for v in vals]
all_timepoints = sorted({x for x in parsed if x is not None})
print(f"All available timepoints: {all_timepoints}")

# To inspect the barcodes, run the following: adata.obs[barcode_col]
# To inspect the times, run the following: adata.obs[time_col]
# To inspect the fates, run the following: adata.obs[fate_col]
```

[!Note] We provided a list of common column names for the barcodes in `barcode_candidates` and time in `time_candidates` to search for. If those common names are not found, the user should manually go through the dataset to identify the correct barcode/clone and time column name. 

The third thing to do is to classify major cell fate groups. Please see `prepare_major_fate_column` function to understand how it's done for our endothetial cell example.

```python
from utilities_fate_prediction import prepare_major_fate_column
adata_raw, fate_col = prepare_major_fate_column(adata, raw_fate_col)
# This function will create a new fate column
```

The fourth thing to do is to filter the data to remove any entry where there is not at least `min_cells_per_tp` many cells at any unique timepoint.

```python
from utilities_fate_prediction import filter_hu_barcodes_for_cross_timepoint
min_cells_per_tp = 2
filtered_adata = filter_hu_barcodes_for_cross_timepoint(
    adata_raw,
    barcode_col=barcode_col,
    time_col=time_col,
    min_cells_per_tp=min_cells_per_tp,
)
```

[!Note] We provided a list of common tokens for missing barcodes in `INVALID_CLONE_TOKENS`. For different datasets, the user needs to check for possible differently defined tokens. 

The fifth thing to do is identify the timepoints of interest. 

```python
# You can choose to use all possible timepoint pairings e.g.
# day 3 vs day 7, day 3 vs day 30, day 7 vs 30
selected_pairs = [
    (all_timepoints[i], all_timepoints[j])
    for i in range(len(all_timepoints) - 1)
    for j in range(i + 1, len(all_timepoints))
]

# For this particular example we will just pick day 3 vs day 30
selected_pairs = [(3.0, 30.0)]

## Handle stage_filter
stage_filter = ""
if stage_filter:
    selected_pairs = [
        (t_early, t_late)
        for (t_early, t_late) in selected_pairs
        if make_stage_tag(t_early, t_late) == stage_filter
    ]
    if not stage_pairs:
        raise ValueError(f"LINEAGE_STAGE_TAG_FILTER={stage_filter} matched no stage pairs")
print(f"Planned stage runs: {selected_pairs}")
```

(Optional) The user may additionally filter for different criterias.

## Running LINGO for Fate Predictions

After preprocessing the data, we are ready to obtain fate predictions.

```python
# All of these folder paths were set inside the LINGO folder
embedding_path = '/home/User/LINGO/finetune_codes/embeddings/saved_results_Danio_rerio/embeddings/kg_gene_embedding_matrix.pt'
gene_to_idx_path = '/home/User/LINGO/finetune_codes/embeddings/saved_results_Danio_rerio/embeddings/gene_to_idx.pkl'
output_path = '/home/User/LINGO/hu_output_folder'
selected_pairs = [(3.0, 30.0)] ######!
model = "finetuned_model"
```

The manuscript discusses many variations of the LINGO model, but we recommend users to the *KG Trainable Embedding* instead of the *KG Scratch* and *KG Frozen*.

Run the following:
```python
from utilities_fate_prediction import fate_prediction
fate_prediction(output_path, embedding_path, gene_to_idx_path, filtered_adata, barcode_col, time_col, fate_col, selected_pairs, model=model)

```
The inputs required are the `output_path` where you wish to store the outputs, `embedding_path` where the pretrained embeddings are stored, `gene_to_idx_path` is the path to the gene alignments, `filtered_adata` is the preprocessed data, `barcode_col` is the column name that contains the barcodes, `time_col` is the column name that contains the timepoints, `fate_col` contains the major cell fate categories, `selected_pairs` are the two time points you wish to use (e.g. day 3 vs day 30), and `model` is the 3 possible models you can select from.  

## Interpreting the outputs
To be completed.


## Reproducing manuscript figures

From a user's perspective the `fate_prediction` function is all you need to know, but if you wish to reproduce to reproduce the figure and downstream benchmark results reported in the manuscript, run the Jupyter notebook code in the folder `LINGO/reproduce_codes/Hu_2022.ipynb`.

This code in the nodebook:
- activates Danio rerio
- loads pretrained zebrafish embeddings
- builds cross-timepoint lineage-origin tasks
- reports macro AUROC and macro AUPRC style outputs consistent with the manuscript

The output will be stored in `LINGO/reproduce_outputs/hu_output`