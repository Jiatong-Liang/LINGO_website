=====
momi3
=====

momi (short for MOran Models for Inference) is a Python package for computing the expected sample frequency spectrum (SFS)—a key summary statistic in population genetics—and using it to infer demographic history.

This third version is a complete rewrite of the original momi and momi2 packages. It introduces several major improvements, including greater flexibility in model specification and enhanced performance and scalability.

Additionally, momi3 is now being built as a core component of the broader ``demestats`` package, providing a more comprehensive framework for demographic inference.

To navigate the ``momi3`` documentation, please refer to the ``Notation`` section first to understand momi3's representation of parameters within a demographic model.

Given any ``demes`` formatted demographic model:

- The ``Tutorial`` goes over all of the core functions of ``momi3``, including how to output and modify model constraints, compute the expected SFS, and compute the likelihood and its gradient.
- The ``Random Projection`` section teaches users how to use random projections to compute an approximation of the full expected SFS and discusses its benefits.
- The ``Optimization`` section demonstrates how to construct custom inference pipelines using ``scipy.minimize``, highlighting momi3's modular design where each component — from objective functions to model constraints — can be tailored to specific research requirements. This flexibility enables researchers to implement specialized optimization strategies.

``momi3`` is described in the following preprint:

Dilber, E., & Terhorst, J. (2024, March 29). Faster inference of complex demographic models from large allele frequency spectra [Preprint]. bioRxiv. https://doi.org/10.1101/2024.03.26.586844
