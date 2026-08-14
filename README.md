# enzymatic_reactions

Analysis scripts for enzymatic assays run on a BioTek Synergy 384-well plate
reader. Each assay type lives in its own top-level folder; today the
repository contains the alcohol dehydrogenase (ADH) DCPIP-linked kinetic
screen.

## Repository layout

```
enzymatic_reactions/
└── adh/                     # DCPIP-linked ADH activity assay
    ├── scripts/             # analysis pipelines
    │   ├── dcpip_standard_curve.py
    │   ├── kcat_plate.py
    │   └── deprecated/      # older single-file versions kept for reference
    ├── test/                # example plate-reader .xlsx exports (input)
    └── outputs/             # generated PDFs / PNGs / CSVs / stats
```

Every script writes into a sibling `outputs/` folder next to the input
`.xlsx`, so the same directory convention works whether you run scripts from
the repo or against a plate stored elsewhere.

## Requirements

- Python 3.10+
- `numpy`, `matplotlib`, `openpyxl`

```bash
pip install numpy matplotlib openpyxl
```

## adh/ — DCPIP-linked alcohol dehydrogenase assay

The assay couples alcohol oxidation to reduction of the redox dye DCPIP
(2,6-dichlorophenolindophenol), which is followed as loss of absorbance at
600 nm. Two scripts turn raw plate-reader exports into publication-ready
panels and QC statistics.

### 1. [adh/scripts/dcpip_standard_curve.py](adh/scripts/dcpip_standard_curve.py)

Fits DCPIP standard curves from a single 384-well plate containing three
independent pH blocks (columns 1–8 = pH 6, 9–16 = pH 7, 17–24 = pH 9). Each
column is a vertical serial dilution with row A at the starting concentration
and row P as the per-column 0 µM blank. A Beer–Lambert line is fit through
the origin, $A_{\text{corr}} = \varepsilon_{\text{app}} \cdot [\text{DCPIP}]$,
independently for each pH. The pH 7 slope is exposed as the primary
`eps_app` that the kcat pipeline consumes.

Serial dilution is specified as a `sample:buffer` ratio (parts transferred vs.
parts already in the destination well). The default is `40:60` — 40 µL passed
between wells into 60 µL of fresh buffer, giving a factor of 0.4 per step.
With the default 300 µM starting concentration, rows A→O of each column carry:

| Row | [DCPIP] (µM) | Row | [DCPIP] (µM) |
| :-: | -----------: | :-: | -----------: |
|  A  |      300     |  I  |      0.197   |
|  B  |      120     |  J  |      0.0786  |
|  C  |       48     |  K  |      0.0315  |
|  D  |       19.2   |  L  |      0.0126  |
|  E  |        7.68  |  M  |      0.00503 |
|  F  |        3.07  |  N  |      0.00201 |
|  G  |        1.23  |  O  |      0.000805|
|  H  |        0.492 |  P  |      0 (blank) |

Example:

```bash
python adh/scripts/dcpip_standard_curve.py adh/test/260731_dcpip_standard_curve.xlsx \
    --start-conc 300 --dilution 40:60
```

Outputs (into `adh/outputs/`):

- `<stem>_standard_curve.pdf` / `.png` — overlaid pH 6/7/9 fits, editable text
- text summary of `eps_app`, standard error, R² and per-column blanks per pH

### 2. [adh/scripts/kcat_plate.py](adh/scripts/kcat_plate.py)

Per-well `k_cat` from a 384-well kinetic export:

1. Parse per-well absorbance-vs-time traces.
2. Subtract the mean water-control trace (row P negative controls) to remove
   baseline drift.
3. Fit an initial-velocity slope `v0` over the linear region — the window
   grows from 4 points up to `--window-s` seconds (or half the run) and keeps
   the longest window that satisfies R² ≥ 0.9.
4. Convert to turnover: `v0_µM/s = −slope / eps_app`,
   `k_cat = v0 / [E]` in s⁻¹.
5. Auto-exclude wells whose ΔA does not clearly exceed the water-control
   distribution ("no enzyme added"; flagged red on the traces grid).
6. Optionally auto-detect a plate loaded rotated 180° from the water/TCEP
   control-row signature and remap well IDs (`--rotate auto|off|180`).

Plate layouts:

- `--layout test` — built-in QC plate: Ca + ethanol everywhere; row P cols
  1–12 = water (negative), cols 13–24 = TCEP (positive; reduces DCPIP
  directly).
- `--layout screen` — real screening plate: 24 metals × 15 alcohols, row P =
  no substrate. Edit the `METALS` / `ALCOHOLS` lists at the top of the script
  to match your protocol.

Quality control includes a `k_cat`-based Z' (positive wells vs. no-substrate)
and a raw-absorbance Z' (TCEP vs. water controls over the first 300 s).

Example runs:

```bash
# Single plate, take eps_app from a paired standard-curve file:
python adh/scripts/kcat_plate.py adh/test/260731_full_plate_test.xlsx \
    --standard-curve adh/test/260731_dcpip_standard_curve.xlsx \
    --layout test --label test_plate

# Triplicate average (well-wise mean/std of k_cat):
python adh/scripts/kcat_plate.py plate1.xlsx plate2.xlsx plate3.xlsx --average \
    --eps 0.006035 --label ADH1_pH7

# Cross-condition compare from a manifest {condition: [plate.xlsx, ...]}:
python adh/scripts/kcat_plate.py --compare manifest.json --eps 0.006035
```

Outputs (into `<xlsx_parent>/../outputs/`):

- `<label>_traces_grid.pdf` / `.png` — 16×24 grid of per-well traces with the
  raw curve, baseline-corrected curve, chosen fit window and `k_cat` printed
  in each panel; control wells shaded (water = blue, TCEP = green).
- `<label>_kcat_heatmap.pdf` / `.png` — Nature-panel-sized heat map with
  row/column condition labels; auto-excluded wells shown red.
- `<label>_wells.csv` — per-well table (`v0`, `k_cat`, R², fit window, n).
- `<label>_stats.txt` — plate-wide QC (means, CVs, both Z' factors, rotation
  status, auto-excluded wells).
- With `--average`: `<label>_avg_heatmaps.pdf` — side-by-side mean / std.
- With `--compare`: `compare_conditions.pdf` and a tidy
  `compare_conditions.csv`.

## Notes

- PDFs are written with editable text (`pdf.fonttype = 42`), so figures can
  be tweaked in Illustrator or Inkscape without re-rasterising.
- `adh/scripts/deprecated/` holds earlier one-shot versions of both scripts;
  they are not imported anywhere and are kept only for reference.