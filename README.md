# 3DFin LiDAR Thinning Pipeline — TLS / MLS / Hand Comparison

An end-to-end R pipeline that takes [3DFin](https://github.com/3DFin/3DFin) tree-detection outputs from **Terrestrial Laser Scanning (TLS)** and **Mobile Laser Scanning (MLS)**, applies competition-index-based thinning decisions, aligns all data sources into a common geodetic reference, and produces comparison maps and agreement statistics against hand-recorded inventory data.

Developed as part of forest inventory research at Freie Universität Berlin.

---

## What the pipeline does

```
TLS CSV ──┐
           ├─► DBH imputation (RF) ─► Competition indices ─► Thinning decisions ─┐
MLS CSV ──┘    (missing values filled       (Hegyi, Braathe,      (worst 30%       │
               by random forest)             RK1–RK4)              local rank)      │
                                                                                    │
Hand CSV ───────────────────────────────────────────────────────────────────────►  │
(or Excel,                                                                          ▼
 polar or UTM)                                                         Similarity ICP alignment
                                                                       (scale + rotate + shift)
                                                                                    │
                                                                                    ▼
                                                                        Overlap map & bar plot
                                                                        Proximity-merged thin layer
                                                                        Cross-media agreement stats
```

### Step-by-step

| Section | What happens |
|---------|-------------|
| **01 DBH imputation / interpolation** | Trees with `DBH ≤ dbh_min_valid_m` (default 1 cm) are treated as missing. Gaps are filled by one of three strategies in priority order: (1) a cross-sensor TLS→MLS transfer model (`tls_transfer`), (2) a local Random Forest fit on measured trees in the same dataset (`rf_imputed`), or (3) the minimum observed DBH in the stand (`fallback_min`). Trees with valid original measurements keep `DBH_source = "measured"`. Final columns: `DBH_meas` (raw), `DBH_amended` (filled value), `DBH_source` (strategy used), `DBH` (analysis-ready; equals `DBH_meas` if measured, else `DBH_amended`). |
| **02 Competition indices** | `TreeCompR` computes six competition indices per tree. Each tree receives a *local neighbourhood percentile rank* within a configurable radius. Trees ranked in the worst 30% are flagged `"remove"`. |
| **03 UTM32 conversion** | Local scan-origin offsets (`x_local`, `y_local`) are converted to absolute UTM32N coordinates by adding a user-supplied origin (`E0`, `N0`). |
| **04 Hand data load** | Supports three input formats: (a) CSV with UTM columns, (b) CSV/Excel with local `x_m`/`y_m` offsets, (c) polar coordinates (`retning`/`avstand` in cm). |
| **05 Similarity ICP alignment** | TLS and MLS point clouds are aligned to the hand reference using an iterative closest point algorithm with similarity (scale + rotation + translation) degrees of freedom. Mutual nearest-neighbour pairs are used at each iteration; scale is clamped to a user-defined guardrail (`icp_scale_min`/`icp_scale_max`). |
| **06 Proximity merge** | Aligned trees from all three sources are matched within `proximity_max_match_m` metres and merged into a single layer recording which media detected / thinned each tree. |
| **07 Visualizations** | A comparison map (colour-coded by detection combination) and a layered bar chart (proportion per method) are written as 300 dpi PNGs. |
| **08 Diagnostics** | Per-medium index agreement (pairwise kappa, Jaccard, Fleiss kappa, bootstrap CIs), alignment quality flags, between-media confusion matrices (sensitivity / specificity / McNemar test) and distance-binned disagreement tables are all written to a timestamped `diagnostics_*/` subfolder. |

---

## Requirements

R ≥ 4.2.0 with the following packages:

```r
install.packages(c(
  "TreeCompR",   # competition indices
  "ranger",      # random forest DBH imputation
  "sf",          # spatial I/O (GeoPackage, CRS transforms)
  "RANN",        # fast nearest-neighbour search (ICP)
  "readr",
  "dplyr",
  "tidyr",
  "purrr",
  "stringr",
  "tibble",
  "rlang",
  "ggplot2",
  "readxl",
  "scales"       # axis formatting in bar plots
))

# Optional — adds north arrow and scale bar to the map:
install.packages("ggspatial")

# Optional — enables interactive mapview preview:
install.packages("mapview")
```

> **Note:** `TreeCompR` may need to be installed from its GitHub repository if not on CRAN.

---

## Input data

Place all input files in a subfolder of `input/` (e.g. `input/Stand 91/`).

### TLS / MLS CSV (3DFin output)

Standard 3DFin export with semicolon or comma delimiter. Required columns (case-insensitive, auto-detected):

| Column | Description |
|--------|-------------|
| `TH` | Tree height (metres) |
| `DBH` | Diameter at breast height (metres, 0 = missing) |
| `X` / `x` | Local easting offset from scan origin (metres) |
| `Y` / `y` | Local northing offset from scan origin (metres) |

The first (unnamed) column is interpreted as the tree ID.

### Hand / field measurement CSV or Excel

Three formats are accepted (auto-detected in `load_hand()`):

1. **Direct UTM columns** — any of `E_utm32`, `E`, `Easting` + `N_utm32`, `N`, `Northing`
2. **Local metric offsets** — `x_m` / `y_m` (combined with `E0`/`N0` to produce UTM)
3. **Polar coordinates** — `retning` (bearing, degrees) + `avstand` (distance, default cm). The raw CRS is set via `crs_hand_raw` (default EPSG:25833 / NTM5) and reprojected to `crs_out`.

An Excel file (`.xlsx` / `.xls`) with polar columns is also supported via `hand_sheet`.

---

## Configuration

All tunable parameters are in the **`#### 00 USER SETTINGS ####`** block at the top of the script. Key settings:

```r
# Paths
tls_in   <- file.path("input", "Stand 91", "TLS_Stand91.csv")
mls_in   <- file.path("input", "Stand 91", "MLS_Stand91.csv")
hand_in  <- file.path("input", "Stand 91", "Hand_Stand91.csv")
out_dir  <- file.path("output", "Stand 91")

# Scan-origin offset for local -> UTM32N
E0 <- 612780.911   # easting of local (0,0)
N0 <- 6618078.542  # northing of local (0,0)

# Thinning: worst 30% by local neighbourhood rank
local_q <- 0.70    # trees ranked >= 0.70 are flagged "remove"

# Alignment guardrails
icp_scale_min <- 0.85
icp_scale_max <- 1.15
match_dmax_m  <- 5.0   # max distance for final hand-to-sensor matching
```

---

## DBH imputation / interpolation

Trees detected by 3DFin may have `DBH = 0` or very small values when the algorithm cannot resolve the stem diameter. The pipeline imputes these in three strategies applied in the following priority order:

| Priority | Strategy | `DBH_source` value | When used |
|----------|----------|--------------------|-----------|
| 1 | **TLS→MLS transfer model** | `tls_transfer` | MLS only; enabled when `mls_use_tls_transfer_model = TRUE`. A Random Forest is trained on TLS trees (with valid DBH) and applied to MLS trees with missing DBH. Avoids propagating MLS-specific scan geometry bias. Predictors are controlled by `mls_transfer_predictors` (default: `TH`, `X`, `Y`). |
| 2 | **Local RF imputation** | `rf_imputed` | Random Forest fit on measured trees within the same dataset. Requires at least 5 training rows; predictors: `TH`, optionally `X`/`Y`. |
| 3 | **Minimum DBH fallback** | `fallback_min` | Assigned when fewer than 5 measured trees are available for RF training. Uses the smallest valid DBH observed in the stand. |

Trees that had a valid original DBH (`DBH > dbh_min_valid_m`, default `0.01` m) are never imputed and receive `DBH_source = "measured"`.

### Output columns produced

| Column | Description |
|--------|-------------|
| `DBH_meas` | Raw DBH from 3DFin (0 or near-0 for unresolved trees) |
| `DBH_amended` | Imputed/interpolated value (filled by one of the three strategies above) |
| `DBH_source` | Strategy that produced the final value: `"measured"`, `"rf_imputed"`, `"tls_transfer"`, or `"fallback_min"` |
| `DBH` | Analysis-ready diameter: equals `DBH_meas` for measured trees, `DBH_amended` otherwise |

All downstream competition index calculations and spatial outputs use `DBH`.

---

## Running the script

```r
# From the repository root (no setwd() required):
source("ver_Braathe_Merged_Script_3DFin_Int_Aft_Combo_Coord_to_UTM32_to_map - Github.R")
```

The script runs sequentially from top to bottom and prints all output file paths on completion.

---

## Outputs

All outputs are written to `out_dir` (default `output/Stand <id>/`).

| File | Description |
|------|-------------|
| `Stand_<id>_TLS_30p_augmented.csv` | TLS trees with imputed DBH, all CI values, local ranks and thinning flags |
| `Stand_<id>_MLS_30p_augmented.csv` | Same for MLS |
| `Stand_<id>_TLS_30p_UTM32.csv` | TLS trees in UTM32N coordinates |
| `Stand_<id>_MLS_30p_UTM32.csv` | Same for MLS |
| `Stand_<id>_HAND_30p_UTM32.csv` | Hand-measured trees in UTM32N |
| `Stand_<id>_TLS_30p_aligned.csv` | TLS trees after ICP alignment to hand reference |
| `Stand_<id>_MLS_30p_aligned.csv` | Same for MLS |
| `trees_aligned_density_Stand<id>_30p_<timestamp>.gpkg` | GeoPackage with layers `HAND_base`, `TLS_aligned`, `MLS_aligned`, `THIN_merged` |
| `Stand_<id>_merged_thin_by_proximity.gpkg / .csv` | Proximity-merged cross-media thinning layer |
| `Stand_<id>_overlap_points.csv` | Per-hand-point detection flags (MLS/TLS/Hand present) |
| `map_thinned_overlap_Stand<id>.png` | Comparison map coloured by detection combination |
| `bar_thinned_overlap_Stand<id>.png` | Layered bar chart of proportions per recording method |
| `CI_plots_local_rank_TLS/` | Local-rank CI maps per competition index (TLS) |
| `CI_plots_local_rank_MLS/` | Same for MLS |
| `diagnostics_<timestamp>/` | Full diagnostics suite (see §08 above) |
| `trees_aligned_..._<timestamp>_settings.txt` | Run manifest with all settings and output paths |

---

## Competition indices

Six indices from `TreeCompR` are computed. All use DBH in metres and height in metres.

| Index | Description |
|-------|-------------|
| `CI_Hegyi` | Distance-weighted size ratio to competitors (Hegyi 1974) |
| `CI_Braathe` | Crown-overlap based index (Braathe 1951) |
| `CI_RK1` – `CI_RK4` | Relative-competition variants included in `TreeCompR` |

Each index gets a `*_rank_local` column (0–1 percentile within the local neighbourhood) and a `thin_*` column (`"remove"` / `"keep"` / `NA`).

---

## Notes on the alignment

The ICP step solves for a **2D similarity transform** (uniform scale + rotation + translation). It uses mutual nearest-neighbour pairs and a trimmed least-squares fit at each iteration. If the estimated scale falls outside `[icp_scale_min, icp_scale_max]` the algorithm falls back to a rigid (scale = 1) solution and sets `fallback_rigid = TRUE` in the diagnostics.

After alignment, trees are clipped to the hand dataset's bounding box and matched to hand trees within `match_dmax_m` metres using a greedy unique-assignment algorithm.

---

## Project structure

```
.
├── ver_Braathe_Merged_Script_3DFin_Int_Aft_Combo_Coord_to_UTM32_to_map - Github.R
├── input/
│   └── Stand 91/
│       ├── TLS_Stand_91_pur.csv
│       ├── MLS_Friday_0045_sc_Katze_Stand91.csv
│       └── 22012026_Stand91_Bäume_mit_XY_Hand.csv
└── output/
    └── Stand 91/            ← created automatically
```

---

## License

MIT License — see [LICENSE](LICENSE).

---

## Author

Janek Wuigk — Freie Universität Berlin
