# Gosai 2024 Synthetic MPRA Evaluator (GAME)

A [GAME](https://github.com/de-Boer-Lab/GAME) Evaluator that requests predictions for **~51k synthetic MPRA sequences from Gosai et al. 2024** (https://www.nature.com/articles/s41586-024-08070-z).
It scores how well a model predicts regulatory activity (expression) across three
cell types — **K562**, **HepG2**, and **SK-N-SH** — by correlating the model's
predictions against measured log2(RNA/DNA) activity.

---

## How It Works

The evaluator:

1. Loads the measured MPRA sequences from `41586_2024_8070_MOESM14_ESM_clean.txt` and the plasmid backbone from `upstream_downstream_backbone.txt`
2. Filters to the synthetic sequences only (`Simulated_Annealing`, `FastSeqProp`, `AdaLead`), which avoids train/test leakage for models trained on the human genome
3. Creates JSON to send to Predictor
4. Negotiates serialization format (JSON or MsgPack) with the Predictor
5. POSTs all sequences to the Predictor's `/predict` endpoint in a single request, across three seperate prediction tasks (K562, HepG2, SK-N-SH)
6. Computes Pearson *r* between predicted expression (log scale) and measured log₂(RNA/DNA) values for each cell type
7. Saves raw predictions and a final evaluation summary CSV

Additional details about the Evaluator's data and the sequence file used can be found in the evaluator_data folder.

---

### Run Evaluator using Apptainer container

Download the container and evaluator data from Zenodo: `<LINK>`

```bash
apptainer run --containall \
  -B /absolute/path/to/evaluator_data:/evaluator_data \
  -B /absolute/path/to/predictions:/predictions \
  gosai_evaluator.sif <predictor_ip> <predictor_port> /predictions
```

---

## Arguments

| Argument | Description |
|---|---|
| `predictor_ip` | IP address or hostname of the predictor server |
| `predictor_port` | Port the predictor is listening on |
| `output_dir` | Directory where prediction JSONs and metric CSVs are written |

---

## Request Structure

A single POST request is sent containing all ~51,000 synthetic sequences. The backbone upstream/downstream flanking sequences and the promoter coordinates (used as `prediction_ranges`) are read from `upstream_downstream_backbone.txt`.

Each request payload follows this structure:

```json
{
  "readout": "point",
  "prediction_tasks": [
    {
      "name": "gosai_synthetic_sequences_k562",
      "type": "expression",
      "cell_type": "K562 (erythroid precursors)",
      "scale": "log",
      "species": "homo_sapiens"
    },
    {
      "name": "gosai_synthetic_sequences_hepg2",
      "type": "expression",
      "cell_type": "HepG2 (hepatocytes)",
      "scale": "log",
      "species": "homo_sapiens"
    },
    {
      "name": "gosai_synthetic_sequences_sknsh",
      "type": "expression",
      "cell_type": "SK-N-SH (neuroblastoma)",
      "scale": "log",
      "species": "homo_sapiens"
    }
  ],
  "upstream_seq": "<plasmid upstream flanking sequence>",
  "downstream_seq": "<plasmid downstream flanking sequence>",
  "sequences": {
    "<seq_id>": "<200 nt sequence>",
    ...
  },
  "prediction_ranges": {
    "<seq_id>": [start, end],
    ...
  }
}
```

**Sequence structure:** Each sequence is 200nt long and corresponds to a synthetic MPRA design. The `prediction_ranges` coordinates correspond to the promoter region within the backbone, as defined in `upstream_downstream_backbone.txt`, and are shared across all sequences.

---

## Outputs

All outputs are written to `output_dir`:

| File | Description |
|---|---|
| `gosai_2024_synthetic_mpra_<timestamp>_predictions_from_<predictor_name>.json` | Raw predictions returned by the predictor |
| `evaluation_summary_gosai_2024_synthetic_mpra_<timestamp>.csv` | Pearson *r* per cell type, appended across runs |

**Metrics computed:**

The measured log₂(RNA/DNA) values for each sequence are stored in the `K562_l2fc`, `HepG2_l2fc`, and `SKNSH_l2fc` columns of the input sequence file in the `evaluator_data` directory.

- **Pearson *r*** — predicted vs. measured log₂(RNA/DNA) for each of the three cell types (K562, HepG2, SK-N-SH)

Cell-type specificity is *not* computed for this dataset.

---

## Directory Structure

```
Gosai_2024_synthetic/
├── evaluator_RestAPI.py                # Main entry point
├── config.py                           # Settings (name, paths, formats, retries)
├── data_loader.py                      # Loads sequence file + backbone, builds request payload
├── evaluator_content_handler.py        # Format negotiation, HTTP POST, deserialization
├── evaluator_metrics_calculator.py     # Pearson r and CSV output
├── gosai_synthetic_evaluator.def       # Apptainer container definition
└── evaluator_data/
    ├── 41586_2024_8070_MOESM14_ESM_clean.txt
    └── upstream_downstream_backbone.txt
```


## The dataset & metric

- **Sequences:** **50,748** synthetic designs from Gosai et al. 2024, each exactly
  **200 bp**, restricted to the `Simulated_Annealing`, `FastSeqProp`, and `AdaLead` design
  methods. These are model-generated rather than genomic, which avoids train/test leakage
  for models trained on the human genome. The breakdown by design method is:
  - `FastSeqProp` — 26,890
  - `Simulated_Annealing` — 11,935
  - `AdaLead` — 11,923

- **Measured values:** per-cell-type log2(RNA/DNA) activity (`K562_l2fc`, `HepG2_l2fc`,
  `SKNSH_l2fc`).
- **Metric:** **Pearson r** between predicted and measured values, computed independently
  per cell type. A task with zero variance in either vector is reported as `r = 0`.
- **Output:** a tab-separated `evaluation_summary_<evaluator_name>.csv` written to the
  mounted predictions directory, with one row per task per metric.

---

## Notes

- This evaluator covers the **synthetic** Gosai sequences only.
- `game_schema_version`: **1.0**
- Authors: Ishika Luthra and Satyam Priyadarshi