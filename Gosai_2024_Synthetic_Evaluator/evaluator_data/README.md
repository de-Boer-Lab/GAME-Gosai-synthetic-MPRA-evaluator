# Gosai Synthetic MPRA Evaluator Data

This folder holds the input data the Evaluator sends to a Predictor and the measured
values it scores against. Two files are used at runtime:

- **`41586_2024_8070_MOESM14_ESM_clean.txt`** — the sequences and their measured
  log2(RNA/DNA) activity.
- **`upstream_downstream_backbone.txt`** — the plasmid backbone (flanking sequences and the
  reporter coordinates) that provides biological context around each mpra probe sequence.

The sequence file was extracted from the raw Gosai et al.[*Machine-guided design of cell-type-targeting cis-regulatory elements*](https://www.nature.com/articles/s41586-024-08070-z). supplementary data.

## `41586_2024_8070_MOESM14_ESM_clean.txt`

The measured-activity table, filtered to drop probes with no measurements.

Starting from the paper's supplementary measurement table (Supplementary
file MOESM14, **77,157** rows), any row with no measured activity in *any* of the three cell
types — i.e. `K562_l2fc`, `HepG2_l2fc`, and `SKNSH_l2fc` all `NA` — was removed. Rows missing
only some cell types were kept. This leaves **76,193** rows.

**Columns** (tab-separated; the first column is the original row index):

| Column | Description |
|---|---|
| `ID` | Unique sequence identifier; used as the prediction key sent to the Predictor |
| `sequence` | The 200 bp variable region |
| `origin` | How the sequence was produced (see below) |
| `target_cell`, `round`, `track_ID` | Design metadata from the paper |
| `K562_l2fc`, `HepG2_l2fc`, `SKNSH_l2fc` | **Measured** log2(RNA/DNA) activity per cell type — the scoring ground truth |
| `MinGap` | Design metadata |
| `K562_lfcSE`, `HepG2_lfcSE`, `SKNSH_lfcSE` | Standard errors of the measured activity |
| `K562_prediction`, `HepG2_prediction`, `SKNSH_prediction` | The paper's own model predictions (not used by this Evaluator) |

**`origin` breakdown** (all 76,193 rows) from Supplementary Table 4 [41586_2024_8070_MOESM6_ESM.xlsx](https://www.nature.com/articles/s41586-024-08070-z):

| origin | count | used for scoring? |
|---|---|---|
| `FastSeqProp` | 26,890 | ✅ synthetic |
| `Simulated_Annealing` | 11,935 | ✅ synthetic |
| `AdaLead` | 11,923 | ✅ synthetic |
| `Malinois_natural` | 11,691 | ❌ genomic |
| `DHS_natural` | 11,616 | ❌ genomic |
| `control` | 2,138 | ❌ control |

At request time the Evaluator keeps **only the three synthetic generators**
(`FastSeqProp`, `Simulated_Annealing`, `AdaLead`) — **50,748** sequences. These are
model-designed rather than genomic, which avoids train/test leakage for Predictors trained
on the human genome. The natural and control sequences are not scored.

---

## `upstream_downstream_backbone.txt`

The plasmid context each 200 bp test sequence is embedded in, plus the coordinate window the
point prediction is read out over. It is a 3-row, tab-separated table:

| Row key | Contents |
|---|---|
| `upstream` | 2,079 bp upstream (5′) plasmid arm |
| `downstream_padded` | 1,397 bp downstream (3′) arm, including a 20 bp `N` placeholder for the barcode region |
| `gene_coordinates` | `[2419, 3171]` — start/end of the reporter (TSS → end of GFP) within the assembled construct |

Derived from the full plasmid vector sequence (3,676 bp). In that vector the
variable insert is marked by a long run of `N`s, so the sequence was split on that run into
the **upstream arm** (everything before the N-run, 2,079 bp) and the **downstream arm**
(everything after, 1,397 bp). The reporter window was located by finding the Puffin TSS motif
(`AATCCG`, at position 2419) and extending by the fixed reporter length of 752 bp (TSS through
the end of GFP), giving `gene_coordinates = [2419, 3171]`.

**Assembled construct.** Predicts for `upstream + sequence + downstream_padded` and reads its point prediction out over
`gene_coordinates`. The same window applies to every sequence.

---

## Notes

- Measured activity is reported on a **log** scale (log2 RNA/DNA); the Evaluator only scores
  a task whose predictions are also returned on the log scale.
- The `_clean` file keeps all sequence origins; the synthetic-only filtering happens in the
  Evaluator at request time, not in this file.
- The paper's Supplementary Table S4 (`41586_2024_8070_MOESM6_ESM.xlsx`) documents how the
  synthetic library was designed. It is kept for reference and is not read by the Evaluator.