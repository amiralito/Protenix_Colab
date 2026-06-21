# Protenix-v2 on Colab: Accessible High-Accuracy Biomolecular Structure Prediction on Google Colab

Two Google Colab notebooks and a helper script for running [Protenix-v2](https://github.com/bytedance/Protenix)
(ByteDance, AlphaFold3-class) end-to-end on Colab, using the **ColabFold MSA server** so **no local
sequence databases** are needed.

Everything runs from forms, predictions and the MSA cache persist to
Google Drive, and the outputs are wired straight into downstream interface scoring (**ipSAE**),
PAE export for **ChimeraX**, and an interactive **MolView** model browser.

* Protenix-v2 is still not publicly accessible but has been integrated into these Colab notebooks. You can still use the previous versions available in the same notebooks.

## Contents

| File | What it is |
| --- | --- |
| `protenix_v2_single.ipynb` | One complex at a time, with full per-model inspection (3D view, pLDDT, interface contacts, PAE, ipSAE). |
| `protenix_v2_batch.ipynb`  | Many complexes in one pass, with aggregate tables (confidence, per-job seed means, ipSAE, interface contacts). |
| `make_protenix_batch.py`   | Generates Protenix batch-input JSON(s) from FASTA files for common screen designs. |

## Open in Colab

[![Single](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/amiralito/Protenix_Colab/blob/main/protenix_v2_single.ipynb) &nbsp;**single prediction**

[![Batch](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/amiralito/Protenix_Colab/blob/main/protenix_v2_batch.ipynb) &nbsp;**batch predictions**

## Requirements

- A Colab GPU runtime. Protenix is memory-light for small inputs (≈6 GB at 500 tokens, ≈18 GB at
  1000), so a **T4** handles monomers and small dimers; use an **A100 (40 GB)** for large complexes
  / resistosomes (≈67 GB at 2000 tokens).
- A Google Drive (optional but recommended — outputs and the MSA cache survive a disconnect).
- Nothing to install locally; `pip install protenix` and the MSA server are handled in the notebook.

The first run on a fresh runtime downloads the model weights automatically.

---

## Single notebook — `protenix_v2_single.ipynb`

Predict and inspect one complex.

1. **1–3b** Setup: GPU check, install, (Blackwell torch if needed), helpers, Drive mount.
2. **4a–4d** Pick one input form: monomer, multimer (≤4 distinct chains), homomer (N copies), or protein + ligand (+ ion).
3. **5** Run settings: model, seeds, samples per seed, dtype, MSA mode.
4. **5b** *(optional)* Enable dense-PAE capture — run **before** inference.
5. **6** MSA (ColabFold server) + inference.
6. **7** Rank predictions by confidence.
7. **8** **MolView** viewer — dropdowns to pick any model, colored by pLDDT, with its scores.
8. **9 / 9b** Per-residue pLDDT; inter-chain interface contacts with cross-seed recurrence (+ per-job CSV).
9. **9c / 9d** PAE export (ChimeraX/ipSAE `.npz` + heatmap) and ipSAE interface scores.
10. **10 / 10b** Zip & download / archive the run.

## Batch notebook — `protenix_v2_batch.ipynb`

Predict and analyze many complexes; unique sequences are MSA-searched once and reused across jobs.

1. **1–3b** Setup (same as above).
2. **5** Run settings.
3. **B1** Load a batch list (built with `make_protenix_batch.py`).
4. **B2-pre** *(optional)* Enable dense-PAE capture — run **before** B2.
5. **B2** MSA + inference over the whole batch.
6. **B3** Per-job summary (best sample, ranked by interface ipTM).
7. **B4** Master confidence table — every model in the batch (CSV).
8. **B5** PAE export — best per seed (ChimeraX/ipSAE `.npz` + heatmap).
9. **B6** ipSAE interface scores (all chain pairs, CSV).
10. **B7** Per-job means across seeds — pTM / ipTM / interface-ipTM / ipSAE (mean ± std).
11. **B8** MolView viewer across all jobs.
12. **B9** Interface contacts — one table per job.
13. **10 / 10b** Zip & download / archive.

---

## `make_protenix_batch.py`

A Protenix input file is a top-level **list** of job dicts; `protenix pred` / `protenix msa`
iterate over every entry. This script builds that list from FASTA files for the common screen
designs, then **B1** in the batch notebook loads it.

```bash
# one job per sequence (monomers)
python make_protenix_batch.py monomer   --fasta seqs.fasta -o batch_inputs

# homo-oligomer: each sequence as an N-mer (e.g. a hexameric resistosome cap)
python make_protenix_batch.py homomer   --fasta nlrs.fasta --copies 6 -o batch_inputs

# every effector × NLR pair (cartesian product of two FASTAs)
python make_protenix_batch.py all_pairs --fasta_a effectors.fasta --fasta_b nlrs.fasta -o batch_inputs

# explicit pairs from a TSV/CSV (columns: idA, idB; sequences pulled from the FASTAs)
python make_protenix_batch.py pairs     --pairs pairs.tsv --fasta_a effectors.fasta --fasta_b nlrs.fasta -o batch_inputs
```

**Modes**

| Mode | Builds | Needs |
| --- | --- | --- |
| `monomer`   | one job per sequence | `--fasta` |
| `homomer`   | each sequence as an `--copies` N-mer | `--fasta --copies` |
| `all_pairs` | cartesian product of two FASTAs (or all within-set pairs of one) | `--fasta_a [--fasta_b]` |
| `pairs`     | explicit `idA,idB` rows from a TSV/CSV | `--pairs --fasta_a [--fasta_b]` |

**Common options** (applied to every job)

- `--ligand CCD_ATP` — a CCD code, a SMILES string, or `FILE_/path.sdf`; `--ligand_copies N`
- `--ion MG` — a bare CCD ion code (e.g. `MG`, `ZN`, `NA`); `--ion_copies N`
- `--copies_a` / `--copies_b` — copy number per partner in `all_pairs` / `pairs`
- `--include_self` — in single-FASTA `all_pairs`, also pair each sequence with itself
- `--chunk 50` — split a large screen into runnable files (`batch_001.json`, `batch_002.json`, …)
- `-o` / `--out_dir`, `--name` — output directory and base filename

Each run also writes a `*_manifest.csv` (json file × job name × number of entities). Job names are
sanitized and de-duplicated automatically.

---

## Outputs

Per run (under `OUT_DIR/run_<date>_<name>/`):

- `<job>/seed_<n>/predictions/` — Protenix `.cif` models + `*_summary_confidence_*.json`
- `<job>/seed_<n>/pae/` — `pae_<model>.npz` (per-token, ChimeraX ≥1.10 + ipSAE), a ColabFold-style
  `.json`, and a heatmap `.png`
- `batch_confidence_all_models.csv`, `batch_confidence_best_per_job.csv` *(batch)*
- `batch_perjob_seed_means.csv` *(batch — pTM/ipTM/ipSAE mean ± std across seeds)*
- `ipsae_all_chain_pairs.csv`, `ipsae_best_per_job.csv`
- `interface_contacts_all_jobs.csv` (+ per-job `interface_contacts.csv`)

**Viewing PAE in ChimeraX** (≥1.10, ligand-safe):
```
open <model>.cif
open pae_<model>.npz structure #1
```
For older ChimeraX / protein-only complexes, open the `*_pae.json` with `format pae` instead.

---

## Notes & tips

- **Blackwell GPUs (sm_120).** Stock wheels lack sm_120 kernels and Protenix's fused LayerNorm
  is precompiled without it. The notebook detects sm_120, sets `LAYERNORM_TYPE=torch`, disables
  `--enable_fusion`, and points you to a one-time nightly cu128 torch install (cell 2c). An **A100**
  needs none of this and is the reliable fallback.
- **MSA.** Uses the public ColabFold server (`api.colabfold.com`); MSAs are cached per sequence and
  reused across jobs and reruns.
- **PAE is per-token.** For protein–protein complexes token = residue, so ChimeraX/ipSAE line up.
  Complexes with ligands have extra atom-tokens — the Boltz-style `.npz` route handles this in both
  ChimeraX (≥1.10) and ipSAE; the per-residue `.json` route does not.
- **ipSAE** runs in Boltz mode (`.npz` + `.cif`), default cutoffs 10/10 (adjustable in the cell).

## Built on

- **Protenix** — ByteDance · <https://github.com/bytedance/Protenix>
- **ColabFold MSA server** — Mirdita et al., *Nat. Methods* (2022)
- **ipSAE** — Dunbrack, *bioRxiv* (2025) · <https://github.com/DunbrackLab/IPSAE>
- **MolView** — Steven Yu · <https://github.com/54yyyu/molview>
- **ChimeraX** — Pettersen et al., *Protein Sci.* (2021)


## References

1. **Protenix-v2.** ByteDance AML AI4Science Team (Protenix Team). *Protenix-v2 Technical Report* (2026). https://github.com/bytedance/Protenix/blob/main/docs/PX2.pdf
2. **Protenix (original).** ByteDance AML AI4Science Team. *Protenix: Advancing Structure Prediction Through a Comprehensive AlphaFold3 Reproduction.* bioRxiv (2025). https://doi.org/10.1101/2025.01.08.631967 <!-- verify DOI before release -->
3. **Protenix-Mini.** Protenix Team. *Protenix-Mini: Efficient Structure Prediction with Compact Architecture and Few-Step Diffusion.* arXiv:2507.11839 (2025). https://arxiv.org/abs/2507.11839
4. **AlphaFold 3.** Abramson J, Adler J, Dunger J, *et al.* Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature* **630**, 493–500 (2024). https://doi.org/10.1038/s41586-024-07487-w
5. **ColabFold.** Mirdita M, Schütze K, Moriwaki Y, Heo L, Ovchinnikov S, Steinegger M. ColabFold: making protein folding accessible to all. *Nature Methods* **19**, 679–682 (2022). https://doi.org/10.1038/s41592-022-01488-1
6. **ipSAE.** Dunbrack RL Jr. *Increasing the accuracy of protein–protein interface confidence estimation with ipSAE.* bioRxiv (2025). https://doi.org/10.1101/2025.02.10.637595 · https://github.com/DunbrackLab/IPSAE <!-- verify title/DOI before release -->
7. **MolView.** Yu S. *MolView: a Mol\*-based molecular viewer for Jupyter/Colab.* Software. https://github.com/54yyyu/molview
8. **UCSF ChimeraX.** Pettersen EF, Goddard TD, Huang CC, *et al.* UCSF ChimeraX: Structure visualization for researchers, educators, and developers. *Protein Science* **30**, 70–82 (2021). https://doi.org/10.1002/pro.3943

## Citation

If you use these notebooks in your work, please cite this repository:

<!-- PLACEHOLDER — fill in after creating the Zenodo release -->
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

> Toghani, A. (2026). *Protenix-v2 on Colab: Accessible High-Accuracy Biomolecular Structure Prediction on Google Colab.* Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX

```bibtex
@software{toghani_protenix_colab_2026,
  author    = {Toghani, AmirAli},
  title     = {{Protenix-v2 on Colab: Accessible High-Accuracy Biomolecular Structure Prediction on Google Colab}},
  year      = {2026},
  publisher = {Zenodo},
  version   = {v1.0.0},
  doi       = {10.5281/zenodo.XXXXXXX},
  url       = {https://doi.org/10.5281/zenodo.XXXXXXX}
}
```

Please also cite **Protenix** [1, 2] (and [3] if you use the Mini model), which these notebooks wrap, and the tools you make use of: **ColabFold** [5] for MSAs, **ipSAE** [6] for interface scoring, **MolView** [7] for visualization, and **ChimeraX** [8] for PAE/structure inspection. Protenix is a reproduction of **AlphaFold 3** [4].

## License

MIT. The bundled tools (Protenix, ipSAE, MolView, ColabFold) keep their own licenses;
please cite them if you use this in published work.
