# Weinreb

The [*Weinreb 2020* dataset](installation.md) tracks mouse hematopoiesis, which is the process in which blood stem cells differentiate into mature blood cells. This tutorial walks through the results of the Weinreb 2020 dataset discussed in the manuscript. 

To see more details regarding the experiment and features of the dataset, please refer to: _Ghosn, E., Yoshimoto, M., Nakauchi, H., Weissman, I. L., & Herzenberg, L. A. (2019). Hematopoietic stem cell-independent hematopoiesis and the origins of innate-like B lymphocytes. Development, 146(15), dev170571._

The example will walk through the **Cross-timepoint linkage prediction task,** which involves predicting whether two cells observed at *different* time points belong to the same clone.

## Data Preprocessing

The first thing to do is to verify whether your dataset has been merged with or modified by other sources. Since our analysis focuses exclusively on Weinreb 2020, we subset the data to retain only samples from that study.

```python
import scanpy as sc

adata = sc.read_h5ad('/home/User/LINGO/data/Weinreb_2020_Science.h5ad') # This is just an example path

if 'dataset' in adata.obs.columns:
    weinreb_mask = adata.obs['dataset'].astype(str).str.contains('Weinreb_2020', na=False)
    if weinreb_mask.sum() > 0:
        adata = adata[weinreb_mask].copy()
        print(f"  Filtered by dataset=Weinreb_2020: {adata.n_obs:,} cells")
```
[!Note] If your `adata` object does not contain the `dataset` column or if you are certain you are not working with a merged dataset, you may skip this step.

The second thing is to identify the column name that contains the barcodes and timepoints.

```python
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

print(f"  Using barcode column: {barcode_col}")
print(f"  Using time column: {time_col}")
# To inspect the barcodes, run the following: adata.obs[barcode_col]
# To inspect the times, run the following: adata.obs[time_col]
```

[!Note] We provided a list of common column names for barcodes in `barcode_candidates` and timepoints in `time_candidates` to search for. If those common names are not found, the user should manually go through their dataset to identify the correct barcode/clone and timepoint column name. 

The third thing to do is to identify the timepoints of interest. To view all timepoints available please run the following:

```python
from utilities_linkage_prediction import (
    _parse_timepoint_value, 
)

vals = adata.obs[time_col].astype(str).values
parsed = [_parse_timepoint_value(v) for v in vals]
parsed_valid = [x for x in parsed if x is not None]
all_timepoints = sorted(set(parsed_valid))
print(f"All parsed timepoints: {all_timepoints}")

# Example:
selected_pairs = [(2.0, 4.0)]
```

The user needs to select which pair of time points to use. For example, if one wants to do *cross-timepoint* linkage prediction of day 2 vs day 4 see `selected_pairs` defined above. If one wants to do *within-timepoint* linkage prediction, then one can set `selected_pairs = [(2.0, 2.0)]` to indicate that only day 2 will be considered.

The fourth thing to do is to filter the data to remove any missing barcodes.

```python
def _valid_clone_mask(series):
    s = series.astype("string")
    return ~(s.isna() | s.str.lower().isin(INVALID_CLONE_TOKENS))

INVALID_CLONE_TOKENS = {"none", "nan", "na", "negative", "neg", "-1", ""}

valid_bc = _valid_clone_mask(bcodes).values
filtered_adata = adata[valid_bc, :].copy()

print(f"Original data size: {adata.shape}")
print(f"Filtered data size: {filtered_adata.shape}")
# In our example, there are no invalid tokens so we are all set
```

[!Note] We provided a list of common tokens for missing barcodes in `INVALID_CLONE_TOKENS`. For different datasets, the user needs to check for possible differently defined tokens. 

(Optional) The user may additionally filter for different criterias, for example in the [Yang Example](Yang_2022.md) we will filter for only tumor cells and ignore all other cells. 

## Model Hyperparameters
Here we provide model hyperparameters that the user can consider. We do not suggest changing any hyperparameter unless the user has familiarity with model.

```python
from config import Config

def build_experiment_hparams():
    """Centrally manage key hyperparameters in the main process to facilitate explicit parameter adjustment."""
    return {
        "shared_pairs": {
            "n_positive": 50000,
            "neg_to_pos_ratio": 1,
            "min_clone_size": 1,
            "train_ratio": 0.7,
            "val_ratio": 0.15,
        },
        "train_common": {
            "use_focal_loss": True,
            "n_positive": 50000,
            "neg_to_pos_ratio": 1,
            "min_clone_size": 3,
            "n_epochs": Config.FINETUNE_EPOCHS,
            "batch_size": 256,
            "mlp_lr": 1e-2,
            "emb_lr": 1e-2,
            "dropout": 0.0,
            "batch_size": 128,
            "reg_to_init_lambda": 0.0,
            "warmup_epochs": 3,
            "resample_interval": 0,
            "emb_freeze_epochs": 0,
            "model_select_metric": "auprc_smooth",
        },
        "finetuned_model": {
            "freeze_embeddings": False,
            "use_adapter": True,
            "downstream_mode": "embedding1",
        },
        "evaluation": {
            "neg_to_pos_ratio": 1,
            "min_clone_size": 2,
            "batch_size": 128,
        },
    }

hparams = build_experiment_hparams()
```

## Running LINGO for Linkage Predictions

After preprocessing the data, we are ready to obtain linkage predictions.

```python
# All of these folder paths were set inside the LINGO folder
embedding_path = '/home/User/LINGO/finetune_codes/embeddings/saved_results_Mus_musculus/embeddings/kg_gene_embedding_matrix.pt'
gene_to_idx_path = '/home/User/LINGO/finetune_codes/embeddings/saved_results_Mus_musculus/embeddings/gene_to_idx.pkl'
output_path = '/home/User/LINGO/weinreb_output_folder'
selected_pairs = [(2.0, 4.0)]
model = "finetuned_model"
```

The manuscript discusses many variations of the LINGO model, but we recommend users to the *KG Trainable Embedding* instead of the *KG Scratch* and *KG Frozen*.

Run the following:
```python
from utilities_linkage_prediction import linkage_prediction
linkage_prediction(output_path, embedding_path, gene_to_idx_path, filtered_adata,
               barcode_col, time_col, selected_pairs, hparams, model)

```
The inputs required are the `output_path` where you wish to store the outputs, `embedding_path` where the pretrained embeddings are stored, `gene_to_idx_path` is the path to the gene alignments, `filtered_adata` is the preprocessed data, `barcode_col` is the column name that contains the barcodes, `time_col` is the column name that contains the timepoints, `selected_pairs` are the two time points you wish to use (e.g. day 2 vs day 4). `hparams` are the hyperparameters, and `model` is the 3 possible models you can select from.  

[!Note] To do *Within-timepoint linkage predictions*, please see the [Yang example](Yang_2022.md). One can specify the same day twice, e.g. `selected_pairs = [(4.0, 4.0)]` to represent that you only care about within-timepoint predictions for day 4. 

## Interpreting the outputs
The `linkage_prediction` function begins by splitting a subset of the Weinreb 2020 dataset into training, validation, and test sets for model development. Once the model is trained, it randomly samples 100,000 unseen cell pairs and predicts cross-timepoint linkages, reporting both AUROC and AUPRC as performance metrics.

After the script completes, the `output_path` folder will contain the following files:

- `*_shared_pairs_cache.pkl` — caches the train/validation/test splits. This allows you to rerun the code with the same random seed without recomputing the splits from scratch.
- `*_ablation_study_results.json` — includes a finetuned_model field that stores the final AUROC and AUPRC for the 100,000 cell pairs, under the keys `shared_test_auc` and `shared_test_auprc`.
- finetuned_model_results/ — a subfolder containing:
    - evaluation_metrics — model performance during training,
    - finetuned_gene_embedding — the learned gene embeddings,
    - finetuned_model_checkpoint — the saved model weights,
    - pair_predictions_for_analysis.csv — a CSV with `label` (true linkage) and `pred_label_best_f1` (predicted linkage) columns.

For most users, the most critical outputs are:

- The final AUROC and AUPRC values, which indicate how well the model fits the data. 
- The pair_predictions_for_analysis.csv file, which provides the predicted cross-timepoint linkages needed for downstream analyses. The columns `label` and `pred_label_best_f1`, which are the true and predicted cross-timepoint linkage results, will be used for downstream analysis. 

## Downstream Analysis


## Reproducing manuscript figures

From a user's perspective the `linkage_prediction` function is all you need to know, but if you wish to reproduce the figure and downstream benchmark results reported in the manuscript, run the Jupyter notebook code in the folder `LINGO/reproduce_codes/Weinreb_2020.ipynb`.

This code in the nodebook:
- activates Mus musculus
- loads pretrained mouse embeddings
- reconstructs cross-timepoint lineage relationships
- writes stage-wise results to the mouse output directory

The output will be stored in `LINGO/reproduce_outputs/weinreb_output`