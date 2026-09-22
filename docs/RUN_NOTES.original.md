# VCC2026 context A/B/C — local-data adaptation

Server: `casia1` (`cbmi_casia@123.184.7.208`). Working directory: `/ssd2/liuqi/xiexiu`.

## Input provenance and differences

The prediction algorithm in `code/model.py` and `code/predict.py` is unchanged from
`/home/cbmi_casia/liuqi/vcc/model/x1/Virtual-Cell-Challenge-2026`.
This is a local-data adaptation, not an exact reproduction of the repository leaderboard score.

| Input | Preparation |
|---|---|
| Official controls, gene and target panels | `/ssd1/PubData/vcc2026-val-1`; controls linked into `data/` |
| H1 2025 training/validation/test | Original repository preparation, local raw files verified against repository SHA256 values |
| CD4 | Original publisher `GWCD4i.DE_stats.h5ad`; original repository preparation and quality filters; SHA256 verified |
| Promoter neighbours | Original GENCODE v47 annotation and original repository code; SHA256 verified |
| K562 | Local `Replogle_K562_gwps/01_preprocess/adata_processed.h5ad`, `layers/counts`; post-QC 1,795,341 cells and 8,248 gene columns (8,246 unique symbols). Barcode numeric suffix restores batch IDs. Duplicate symbols are summed. All retained controls are used. CPM denominator is the retained gene row sum. |
| HCT116 / HEK293T | Local `XAtlas_*/01_preprocess/adata_processed.h5ad`, `layers/counts`; post-QC cells. First 250 retained controls per sample in local row order. Original remote `cell_integer_id` ordering is unavailable. Gene projection and sample-matched control weighting follow the repository formulas. |

All three adapted sources use the union of official targets and H1 training targets,
as the original preparation does. Full available source panels are retained for
source centering. Missing genes are explicitly marked unmeasured rather than
treated as observed zero responses. Sources with fewer than 20 target cells are
excluded by the original model.

No additional raw datasets are downloaded by these preparation or prediction scripts.
The original repository copy retains its download function, but the active workflows
use local files only. `prepare_original_local.py` replaces the download function
with an existence, size and SHA256 checker.

The local source adapter reads coalesced count intervals around selected cells
to avoid decompressing unrelated perturbations. A direct full-slice comparison
verified that this returns identical count and gene-index values.

## Execution

Data preparation/inference runtime: `env/bin/python` (Python 3.10, NumPy 2.2.6,
SciPy 1.15.3, pandas 2.3.3, anndata 0.11.4, h5py 3.16.0, numba 0.67.0).
These are not the repository's pinned versions. Input and output validation is
therefore performed explicitly. VCC packaging uses an independent Python 3.11
environment, `vcc-env`, with `vcc-cli==0.2.0`.

Preparation scripts:

```bash
cd /ssd2/liuqi/xiexiu
env/bin/python prepare_original_local.py
env/bin/python run_local_sources.py
```

The pipeline waits for local source preparation, audits all inputs, runs the
original predictor, and runs official VCC prep directly on the original CSR H5AD:

```bash
OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 NUMBA_NUM_THREADS=4 \
  env/bin/python run_prediction_pipeline.py
```

Do not start another pipeline while the existing one is running. It refuses to
overwrite an unfinished prediction and records failures in its status JSON.

## Expected deliverables and verification

`output/prediction_ABC_local.h5ad` contains all three contexts. Each context has
300 targets with 400 predicted cells each: 360,000 cells total, 18,533 genes.
`output/prediction_ABC_local.vcc` is the submission package. The user subsequently
authorized official submission using the current Wang Yuhang team account.
`submit_official.py` waits for validated artifacts, verifies the approved account,
and submits under `Xiexiu-X1-LocalAtlas-ABC-20260915`. It records an attempt marker
to prevent accidental duplicate submissions after an ambiguous response.

Optional observation-reordering compaction was stopped because it was very slow
on this host. No compact H5AD is a deliverable. Official `pack.py` supports the
original CSR H5AD and performs the same official checks without this optional step.
Packaging uses a task-specific temporary directory under `/dev/shm` to avoid
unnecessary intermediate-file I/O on the network-mounted `/ssd2`.

Checks include:

- Input control labels, gene order, nonnegative integer counts and positive depths.
- Source probability normalization, finite effects, measured gene masks and target coverage.
- Predictor moment-fitting tolerances, exact per-cell depth and repaired pseudobulk counts.
- Official VCC prep checks.
- Final 900 context/target groups, exactly 400 cells per group, official gene order,
  and output SHA256 hashes.

Authoritative completion evidence is `logs/final_validation.json` with
`state: complete`, not the presence of a partial output file. Current stages are
recorded in `logs/local_sources_status.json` and `logs/prediction_pipeline_status.json`.
Source coverage is written to `data/source_coverage.csv` after input validation.

