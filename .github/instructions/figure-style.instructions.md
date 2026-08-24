---
description: "Use when creating new Python scripts that produce figures or plots. Enforces Nature-journal figure style, matplotlib rcParams, color palette, file output format, and panel sizing conventions used throughout this project."
applyTo: "**/*.py"
---

# Figure Style Guidelines

All figures produced by scripts in this project must conform to Nature's figure
preparation guidelines and match the visual style established in `adh/scripts/kcat_plate.py`.

---

## 1. rcParams — always apply at the top of every plotting function

```python
NATURE_RC = {
    "svg.fonttype": "none",        # editable text in Inkscape/Illustrator
    "pdf.fonttype": 42,            # TrueType fonts embedded in PDF
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 6,
    "axes.labelsize": 6,
    "axes.titlesize": 6,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size": 2.0,
    "ytick.major.size": 2.0,
    "axes.spines.top": False,      # no top spine
    "axes.spines.right": False,    # no right spine
}
```

Call `plt.rcParams.update(NATURE_RC)` at the start of every plotting function,
before creating any figure or axes.

---

## 2. Figure sizes

| Use case | Size (inches) | Size (mm) | Notes |
|---|---|---|---|
| Standard single/double-column panel | `(7.08, 5.20)` | 180 × 132 | Double-column, half-page |
| Tall grid (e.g. 16×24 well plate) | `(7.08 * 1.5, 7.08 * 1.5 * 16/24)` | ~270 × 180 | Keep aspect ratio |
| Side-by-side panels (×2) | `(7.08 * 2, 5.20)` | 360 × 132 | Each panel same width |
| Side-by-side panels (×3) | `(7.08 * 3, 5.20)` | — | Shared colour scale |

Single-column figure maximum width is 88 mm (3.46 in); double-column is 180 mm
(7.08 in). Never exceed 180 mm wide or 170 mm tall for a main figure panel.

---

## 3. Colour palette

Use these project-standard colours consistently across all scripts:

| Role | Hex | Notes |
|---|---|---|
| Primary data line / fit | `#1f77b4` | matplotlib tab blue |
| Raw / secondary trace | `"0.55"` | mid-grey string |
| Fit window shading | `"0.8"` | light grey, alpha 0.5 |
| Z′ region shading | `#7a9cc6` | muted blue, alpha 0.35 |
| Water-control background | `#dce9f7` | light blue |
| TCEP-control background | `#d9ecd0` | light green |
| Excluded / bad-well background | `#ffd6d6` | light red |
| Control-row separator / divider | `"0.3"` | dark grey, lw 0.8, ls "--" |
| Control-row band fill | `#e8e8e8` | very light grey |
| Bad-data NaN cells in heatmap | `#ff5252` | red |

For heatmaps use `gray_r` as the default colormap with `cmap.set_bad("#ff5252")`.
Use `viridis` when colouring individual text annotations by value magnitude.

---

## 4. Line weights and marker sizes

| Element | Setting |
|---|---|
| Data traces | `lw=0.4` |
| Fit / model lines | `lw=0.7` |
| Fit scatter points | `ms=0.8`, `mec="none"` |
| Axis spine / frame | `linewidth=0.5` (set via rcParams) |
| Tick / colorbar outline | `linewidth=0.4` |
| Control-row divider | `lw=0.8`, `ls="--"` |
| Legend line examples | `lw=1.0` |

---

## 5. Annotations and text

- All in-panel text (k_cat values, labels): `fontsize=4`
- Figure-level annotations (Z′, rotation warnings): `fontsize=7`
- Suptitle / panel title: `fontsize=7`
- Colorbar label: set via `cb.set_label(...)` — inherits rcParams size (6 pt)
- Use LaTeX for physical quantities: `k$_{cat}$`, `Z'$_{raw}$`, `µM`, `s$^{-1}$`
- Control-row y-tick: italic, colour `"0.4"`, label text `"controls"`

---

## 6. Axes and layout

- Remove top and right spines (enforced by rcParams; do not re-add them).
- Use `fig.tight_layout()` for single-panel figures.
- Use explicit `fig.subplots_adjust(wspace=0.05, hspace=0.05, top=0.93, bottom=0.02, left=0.03, right=0.99)` for dense grid figures.
- Colourbar: `shrink=0.85`, `pad=0.02`; tick params `labelsize=6`; outline `linewidth=0.4`.
- Legends: `frameon=False`; place with `loc="upper center"` + `bbox_to_anchor` for grid figures.

---

## 7. File output

Always save **both** PDF and PNG from every plotting function:

```python
out_path.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(out_path, format="pdf")
fig.savefig(out_path.with_suffix(".png"), format="png", dpi=600)
plt.close(fig)
```

- PDF: vector, editable fonts (`pdf.fonttype=42`), for submission.
- PNG: 600 dpi raster, for quick review and presentations.
- Always call `plt.close(fig)` to free memory.
- Output directory: `<plate_dir>/../outputs/` (i.e. `_default_outdir` convention).

---

## 8. Nature figure checklist

Before finalising any figure:

- [ ] All text ≥ 5 pt at final print size (6 pt default satisfies this).
- [ ] No text rendered as outlines — `svg.fonttype="none"`, `pdf.fonttype=42`.
- [ ] Fonts: Helvetica or Arial only (DejaVu Sans as fallback).
- [ ] Figure width ≤ 180 mm (7.08 in) for double-column; ≤ 88 mm for single-column.
- [ ] No top or right spines.
- [ ] Colourbars labelled with units.
- [ ] NaN / missing data shown in a distinct colour (`#ff5252`), not white.
- [ ] Control wells / rows visually separated from assay data.
- [ ] Z′ value annotated on the figure when available.
- [ ] Both PDF and PNG saved at 600 dpi.
