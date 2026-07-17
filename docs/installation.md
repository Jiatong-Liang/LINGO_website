# Installation/Quickstart

## 🚀 Quickstart

To reproduce the downstream benchmark experiments described in the manuscript:

1. Create the Conda environment from `environment.yaml`
2. Download **LibKG** from Zenodo
3. Download the four downstream raw `.h5ad` benchmark files
4. Use the pretrained embeddings already provided in `LINGO/finetune_codes/embeddings`
5. Run the corresponding scripts in `LINGO/finetune_codes/`. Please see the examples. 

---
## Installation

Start by cloning the github:

```bash
git clone git@github.com:Young0222/LINGO.git
```

Please use virtual environment to manage dependencies for `LINGO` to avoid conflicts with other Python packages. 

Using `venv`:

```bash
python -m venv lingo-env
source lingo-env/bin/activate 
pip install LINGO_requirements.txt
```

Using `conda`:

```bash
conda env create -f environment.yaml
conda activate lingo-pretrain
```

For a sanity check:

```bash
python -c "import torch, scanpy, torch_geometric; print('Environment is ready')"
```

[NOTE!] The provided `environment.yaml` and `requirements.txt` is intended to reproduce the LINGO pipeline. If your machine uses a different CUDA version, you may need to adjust the PyTorch installation accordingly.

---

## 📦 Data Download

The released **LibKG** dataset is available on Zenodo, please download `LibKG.zip`:

[LibKG on Zenodo](https://zenodo.org/records/20434093)

This Zenodo release contains the **species-specific LibKG resources** used by `LINGO`, including graph node tables, edge tables, and metadata for:

- **Mus musculus**
- **Homo sapiens**
- **Danio rerio**

These task-specific raw `.h5ad` files are available at:

[https://zenodo.org/records/12176634](https://zenodo.org/records/12176634)

Please download the following four files:

- `Hu_2022_NatureGenetics.h5ad`
- `Hurley_2020_CellStemCell.h5ad`
- `Weinreb_2020_Science.h5ad`
- `Yang_2022_Cell.h5ad`

In addition, please download this data file from the following [https://zenodo.org/records/14002044](https://zenodo.org/records/14002044):
- MECOM_KO.seurat.rds

[!IMPORTANT] The LibKG Zenodo record contains the knowledge graph resources only. The downstream finetuning benchmark `.h5ad` files are provided separately at the [Zenodo link](https://zenodo.org/records/12176634). The [pretrained embedding](https://github.com/Young0222/LINGO/tree/main/finetune_codes/embeddings) files used by finetuning are already available in the GitHub repository.

After downloading, unpack the archive and place the files according to the layout below.

```text
LINGO/
├── environment.yaml
├── requirements.txt
├── merged_by_species/ (this folder is optional)
    ├── merged_Mus_musculus.h5ad
    ├── merged_Homo_sapiens.h5ad
    ├── merged_Danio_rerio.h5ad
├── reproduce_codes/
├── pretrain_codes/
├── finetune_codes/
├── data/
    ├── Weinreb_2020_Science.h5ad
    ├── Yang_2022_Cell.h5ad
    ├── Hu_2022_NatureGenetics.h5ad
    ├── Hurley_2020_CellStemCell.h5ad
    ├── MECOM_KO.seurat.rds
```

[Note!] The pretraining pipeline expects three merged `.h5ad` files under `merged_by_species/` but this folder is **optional**. You can continue without this folder if you do not want to train the model from scratch.  