# Installation

Start by cloning the github:

```bash
git clone git@github.com:Young0222/LINGO.git
```

We recommend using a virtual environment to manage dependencies for `LINGO`. This
helps avoid conflicts with other Python packages and ensures a clean installation.

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

> [!NOTE]
> The provided `environment.yaml` is intended to reproduce the LINGO pipeline. If your machine uses a different CUDA version, you may need to adjust the PyTorch installation accordingly.

