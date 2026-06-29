# Weinreb

The [*Weinreb 2020* dataset](installation.md) tracks mouse hematopoiesis, which is the process in which blood stem cells differentiate into mature blood cells. This tutorial walks through the results of the Weinreb 2020 dataset discussed in the manuscript. 

To see more details regarding the experiment and features of the dataset, please refer to: _Ghosn, E., Yoshimoto, M., Nakauchi, H., Weissman, I. L., & Herzenberg, L. A. (2019). Hematopoietic stem cell-independent hematopoiesis and the origins of innate-like B lymphocytes. Development, 146(15), dev170571._

The example will walk through the **Cross-timepoint linkage prediction task,** which involves predicting whether two cells observed at *different* time points belong to the same clone.

The manuscript discusses many variations of the LINGO model, but we recommend users to the *KG Trainable Embedding* instead of the *KG Scratch* and *KG Frozen*.

Run the following:
```python
run_experiment3(data_path='/home/jkliang/LINGO/data/Weinreb_2020_Science.h5ad', output_path='/home/jkliang/LINGO/weinreb_results', stage_pairs=[(2.0, 4.0)])
```
The inputs required are the `data_path` to the Weinreb 2020 dataset, the `output_path` is where you wish to store the outputs, and `stage_pairs` are the two time points you wish to use (e.g. day 2 vs day 4). If you want to run multiple different time points with one line of code, one can set `stage_pairs` to something like 

```python
stage_pairs = [(2.0, 4.0), (4.0, 6.0), (2.0, 6.0)]
```

to get comparisons of day 2 vs day 4, day 4 vs day 6, and day 2 vs day 6, etc.

## Interpreting the outputs
The `run_experiment3` code will take a portion of the Weinreb 2020 dataset and split in into train, validation, and test datasets to create the model. After creating the model, 100,000 cell pairings will be randomly taken and be used to predict cross-timepoint linkage. 

After running the code, within the `output_path` folder you should see the `experiment_3_finetuned_improved` folder which contains `day2_day4_pair_predictions_for_analysis.csv`. Within the csv file, the columns `label` and `pred_label_best_f1` are the true and predicted cross-timepoint linkage results.The code also reports the AUC and AUPRC.


## Reproducing manuscript figures

From a user's perspective the `run_experiment3` function is all you need to know, but if you wish to reproduce to reproduce the figure and downstream benchmark results reported in the manuscript, run the following code from the repository root (LINGO folder) in your terminal:

Set the raw data path:
```bash
export LINEAGE_WEINREB_PATH=/absolute/path/to/Weinreb_2020_Science.h5ad
```

Run:
```bash
python finetune_codes/10_lineage_reconstruction_Weinreb_2020.py
```

This script:
- activates Mus musculus
- loads pretrained mouse embeddings
- reconstructs cross-timepoint lineage relationships
- writes stage-wise results to the mouse output directory