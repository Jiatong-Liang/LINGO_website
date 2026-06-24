# Welcome to LINGO documentation

<p align="center">
  <strong>A Knowledge Graph-Grounded Foundation Model for Single-Cell Lineage Inference</strong>
</p>

<p align="center">
  Multi-species pretraining • biological knowledge graph • lineage-aware transfer learning
</p>

---

## 🧬 Overview

LINGO is a knowledge graph-grounded foundation model for **single-cell lineage inference**.   It is designed to learn transferable gene and cell representations from large-scale single-cell lineage tracing data, while explicitly incorporating structured biological priors.

The method for LINGO is described in the following preprint: (insert citation)

LINGO is described as:

- a **two-stage framework** with graph-based pretraining and task-specific finetuning
- a model pretrained on **18 tissue types** across **mouse, human, and zebrafish**
- a framework evaluated on four held-out systems:
  - **mouse hematopoiesis** from Weinreb 2020
  - **mouse tumor evolution** from Yang 2022
  - **human cell line clone plasticity** from Hurley 2020
  - **zebrafish heart lineage origin prediction** from Hu 2022

LINGO supports three core downstream tasks:

- **cell linkage prediction**
- **lineage origin prediction**
- **clonal plasticity prediction**

To get started, please refer to [`Installation/Quickstart`](installation.md).

```{toctree}
:hidden:
:maxdepth: 4

Installation/Quickstart <installation>
```