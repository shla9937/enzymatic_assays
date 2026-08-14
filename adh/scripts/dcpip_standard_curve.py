"""Fit DCPIP absorbance standard curves from a BioTek 384-well plate export.

Layout assumed (per updated user protocol):
    * Vertical titrations down each column (16 wells per column).
    * Row A holds the starting concentration; rows B..O are serial dilutions;
      row P is 0 [DCPIP] (per-column blank).
    * Columns 1-8   -> pH 6      (8 replicates)
    * Columns 9-16  -> pH 7      (8 replicates)
    * Columns 17-24 -> pH 9      (8 replicates)
    * Non-numeric readings (e.g. "OVRFLW") are dropped automatically.

Beer-Lambert is fit as A_corrected = eps_app * [DCPIP], forced through origin,
independently for each pH.

Serial dilution:
    Each step transfers ``sample`` parts of the previous (concentrated) well
    into ``buffer`` parts of fresh diluent in the next well, so
    C_new = C_old * sample / (sample + buffer). Passed as ``--dilution sample:buffer``
    (default ``1:1``, i.e. equal parts → factor 1/2 = 2-fold dilution per step).

Usage:
    python dcpip_standard_curve.py path/to/standard_curve.xlsx
    python dcpip_standard_curve.py path/to/standard_curve.xlsx \\
        --start-conc 300 --dilution 6:1
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import openpyxl


DEFAULT_START_CONC_UM: float = 300.0
DEFAULT_DILUTION_RATIO: str = "1:1"  # sample:buffer -> factor 1/2 per step (2-fold)
DEFAULT_OUTPUT_SUBDIR: str = "outputs"

# (pH, first_col, last_col_inclusive); pH 7 is used as the "primary" fit that
# downstream code (kcat_plate.py) picks up via StandardCurveResult.eps_app.
PH_BLOCKS: list[tuple[float, int, int]] = [
    (6.0, 1, 8),
    (7.0, 9, 16),
    (9.0, 17, 24),
]
PRIMARY_PH: float = 7.0

PLATE_ROWS: list[str] = list("ABCDEFGHIJKLMNOP")  # 16 rows
BLANK_ROW: str = "P"                                # 0 [DCPIP]
N_ROWS: int = len(PLATE_ROWS)                       # 16

PH_COLORS: dict[float, str] = {
    6.0: "#1f77b4",  # blue
    7.0: "#2ca02c",  # green
    9.0: "#d62728",  # red
}


# =============================================================================
# Data classes
# =============================================================================

@dataclass
class PhFit:
    """Standard-curve fit for a single pH block."""

    ph: float
    concentrations_uM: np.ndarray       # per-well concentrations included
    absorbance_corrected: np.ndarray    # per-well (raw - column blank)
    replicate_labels: np.ndarray        # column number per point (str)
    eps_app: float                      # slope, AU / µM
    eps_app_stderr: float
    r_squared: float
    n_points: int
    blanks: dict[int, float]            # per-column row-P absorbance actually used


@dataclass
class StandardCurveResult:
    """Aggregate of per-pH fits; exposes primary-pH eps_app for downstream reuse."""

    per_ph: dict[float, PhFit] = field(default_factory=dict)
    primary_ph: float = PRIMARY_PH
    start_conc_uM: float = DEFAULT_START_CONC_UM
    dilution_ratio: str = DEFAULT_DILUTION_RATIO
    dilution_factor: float = 1.0 / 7.0
    source: Path | None = None

    def _primary(self) -> PhFit:
        if self.primary_ph in self.per_ph:
            return self.per_ph[self.primary_ph]
        # Fall back to any available fit so callers don't blow up.
        return next(iter(self.per_ph.values()))

    @property
    def eps_app(self) -> float:
        return self._primary().eps_app

    @property
    def r_squared(self) -> float:
        return self._primary().r_squared

    def absorbance_to_concentration(self, absorbance, blank: float = 0.0,
                                    ph: float | None = None):
        eps = self.per_ph[ph if ph is not None else self.primary_ph].eps_app
        return (np.asarray(absorbance) - blank) / eps


# =============================================================================
# CLI / main
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("xlsx", type=Path, help="Path to the plate-reader .xlsx export.")
    parser.add_argument("-o", "--output", type=Path, default=None,
                        help=f"Output PDF path (default: <xlsx_parent>/../"
                             f"{DEFAULT_OUTPUT_SUBDIR}/<stem>_standard_curve.pdf).")
    parser.add_argument("--start-conc", type=float, default=DEFAULT_START_CONC_UM,
                        help=f"Starting [DCPIP] in row A (µM). Default {DEFAULT_START_CONC_UM}.")
    parser.add_argument("--dilution", type=str, default=DEFAULT_DILUTION_RATIO,
                        help="Per-step dilution as 'sample:buffer' ratio. "
                             f"Default {DEFAULT_DILUTION_RATIO} (equal parts → factor 1/2, "
                             "i.e. 2-fold per step).")
    parser.add_argument("--sheet", type=str, default=None,
                        help="Sheet name (default: first sheet).")
    parser.add_argument("--title", type=str, default="",
                        help="Optional plot title.")
    args = parser.parse_args()

    result = fit_standard_curve(
        args.xlsx,
        start_conc_uM=args.start_conc,
        dilution_ratio=args.dilution,
        sheet=args.sheet,
    )

    out_pdf = args.output
    if out_pdf is None:
        out_dir = args.xlsx.resolve().parent.parent / DEFAULT_OUTPUT_SUBDIR
        out_pdf = out_dir / f"{args.xlsx.stem}_standard_curve.pdf"

    plot_standard_curve(result, out_pdf, title=args.title or args.xlsx.stem)
    _print_summary(result)
    print(f"\nWrote editable PDF: {out_pdf}")


# =============================================================================
# Dilution helper
# =============================================================================

def parse_dilution_ratio(s: str) -> float:
    """Parse a ``sample:buffer`` ratio into a per-step factor ``C_new / C_old``.

    Each step adds ``sample`` parts of the previous well into ``buffer`` parts
    of fresh diluent, so ``C_new = C_old * sample / (sample + buffer)``.
    Examples: ``'1:1'`` -> 1/2, ``'1:5'`` -> 1/6, ``'1:9'`` -> 1/10.
    A bare number in (0, 1) is taken as the factor directly.
    """
    txt = s.strip()
    if ":" in txt:
        sample_str, buffer_str = txt.split(":", 1)
        sample = float(sample_str)
        buf = float(buffer_str)
        if sample <= 0 or buf < 0:
            raise ValueError(f"Bad dilution ratio: {s!r}")
        return sample / (sample + buf)
    factor = float(txt)
    if not (0 < factor < 1):
        raise ValueError(f"Dilution factor must be in (0, 1): {s!r}")
    return factor


# =============================================================================
# Plate loader (BioTek Synergy export)
# =============================================================================

def _find_plate_grid(ws) -> tuple[int, int]:
    """Return (header_row, first_data_col) for the '1..24' column header row."""
    for row_idx, vals in enumerate(ws.iter_rows(values_only=True), start=1):
        for start in range(len(vals) - 1):
            if vals[start] == 1 and vals[start + 1] == 2:
                # Row labels live one column to the left of the "1" column.
                return row_idx, start
    raise ValueError("Could not find plate column-header row (1, 2, ...) in sheet.")


def _load_plate(xlsx_path: Path, sheet: str | None = None) -> dict[str, list]:
    """Load the plate as {row_label: [24 raw cell values]}."""
    wb = openpyxl.load_workbook(xlsx_path, data_only=True, read_only=True)
    ws = wb[sheet] if sheet else wb[wb.sheetnames[0]]
    header_row, first_data_col = _find_plate_grid(ws)

    plate: dict[str, list] = {}
    for row in ws.iter_rows(min_row=header_row + 1, values_only=True):
        label = row[first_data_col - 1]
        if not isinstance(label, str) or len(label) != 1 or not label.isalpha():
            continue
        plate[label.upper()] = list(row[first_data_col:first_data_col + 24])
    return plate


def _to_float_or_nan(x) -> float:
    if isinstance(x, (int, float)) and not isinstance(x, bool):
        return float(x)
    return np.nan


# =============================================================================
# Fit
# =============================================================================

def _fit_single_ph(
    plate: dict[str, list],
    ph: float,
    col_lo: int,
    col_hi: int,
    concs_by_row: dict[str, float],
) -> PhFit:
    """Fit A = eps_app * [DCPIP] through the origin for one pH block."""
    xs: list[float] = []
    ys: list[float] = []
    labels: list[str] = []
    blanks: dict[int, float] = {}

    for col in range(col_lo, col_hi + 1):
        # Per-column blank from row P.
        blank_val = _to_float_or_nan(plate[BLANK_ROW][col - 1])
        if np.isnan(blank_val):
            raise ValueError(f"Blank well {BLANK_ROW}{col} (pH {ph}) is not numeric.")
        blanks[col] = blank_val

        for row_label in PLATE_ROWS:
            if row_label == BLANK_ROW:
                continue
            raw = _to_float_or_nan(plate[row_label][col - 1])
            if np.isnan(raw):
                continue
            xs.append(concs_by_row[row_label])
            ys.append(raw - blank_val)
            labels.append(str(col))

    x = np.asarray(xs)
    y = np.asarray(ys)

    # Linear regression forced through origin: y = m*x.
    sxx = float(np.sum(x * x))
    sxy = float(np.sum(x * y))
    slope = sxy / sxx if sxx > 0 else np.nan
    residuals = y - slope * x
    dof = len(x) - 1
    resid_var = float(np.sum(residuals ** 2)) / dof if dof > 0 else np.nan
    slope_se = float(np.sqrt(resid_var / sxx)) if dof > 0 and sxx > 0 else np.nan
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    ss_res = float(np.sum(residuals ** 2))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    return PhFit(
        ph=ph,
        concentrations_uM=x,
        absorbance_corrected=y,
        replicate_labels=np.asarray(labels),
        eps_app=slope,
        eps_app_stderr=slope_se,
        r_squared=r_squared,
        n_points=len(x),
        blanks=blanks,
    )


def fit_standard_curve(
    xlsx_path: str | Path,
    start_conc_uM: float = DEFAULT_START_CONC_UM,
    dilution_ratio: str = DEFAULT_DILUTION_RATIO,
    ph_blocks: list[tuple[float, int, int]] = PH_BLOCKS,
    primary_ph: float = PRIMARY_PH,
    sheet: str | None = None,
) -> StandardCurveResult:
    """Load the plate, split into pH blocks, and fit each independently."""
    xlsx_path = Path(xlsx_path)
    plate = _load_plate(xlsx_path, sheet=sheet)

    dilution_factor = parse_dilution_ratio(dilution_ratio)
    # Row A = step 0 (start conc); Row O = step 14; Row P = 0 (blank).
    concs_by_row: dict[str, float] = {}
    for step, row_label in enumerate(PLATE_ROWS[:-1]):  # A..O
        concs_by_row[row_label] = start_conc_uM * (dilution_factor ** step)
    concs_by_row[BLANK_ROW] = 0.0

    per_ph: dict[float, PhFit] = {}
    for ph, col_lo, col_hi in ph_blocks:
        per_ph[ph] = _fit_single_ph(plate, ph, col_lo, col_hi, concs_by_row)

    return StandardCurveResult(
        per_ph=per_ph,
        primary_ph=primary_ph,
        start_conc_uM=start_conc_uM,
        dilution_ratio=dilution_ratio,
        dilution_factor=dilution_factor,
        source=xlsx_path,
    )


# =============================================================================
# Plot
# =============================================================================

def plot_standard_curve(result: StandardCurveResult, out_pdf: Path,
                        title: str = "") -> None:
    """Render an editable PDF with all pH curves overlaid + per-pH stat box."""
    plt.rcParams.update({
        "svg.fonttype": "none",   # editable text in SVG (Illustrator/Inkscape)
        "pdf.fonttype": 42,       # editable text in PDF (Illustrator/Inkscape)
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
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    # ~90 x 68 mm; fits one quadrant of a 2x2 Nature double-column figure.
    fig, ax = plt.subplots(figsize=(3.55, 2.68))

    x_max = 0.0
    stats_lines: list[str] = []
    for ph in sorted(result.per_ph):
        fit = result.per_ph[ph]
        color = PH_COLORS.get(ph, "black")
        if fit.concentrations_uM.size:
            x_max = max(x_max, float(np.max(fit.concentrations_uM)))
            ax.scatter(
                fit.concentrations_uM,
                fit.absorbance_corrected,
                s=8, facecolors="none", edgecolors=color, linewidths=0.5,
                label=f"pH {ph:g}", zorder=3,
            )
        stats_lines.append(
            f"pH {ph:g}: "
            f"$\\varepsilon$={fit.eps_app:.4g}"
            f"±{fit.eps_app_stderr:.2g} AU/µM, "
            f"R²={fit.r_squared:.4f}, n={fit.n_points}"
        )
    if x_max <= 0:
        x_max = result.start_conc_uM

    x_line = np.linspace(0.0, x_max * 1.05, 200)
    for ph in sorted(result.per_ph):
        fit = result.per_ph[ph]
        if not np.isfinite(fit.eps_app):
            continue
        ax.plot(x_line, fit.eps_app * x_line,
                color=PH_COLORS.get(ph, "black"), lw=0.6, zorder=2)

    stats_text = (
        "\n".join(stats_lines)
        + f"\nstart = {result.start_conc_uM:g} µM, "
        + f"dilution {result.dilution_ratio} (factor {result.dilution_factor:.4g})"
    )
    ax.text(
        0.98, 0.02, stats_text,
        transform=ax.transAxes, ha="right", va="bottom",
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                  edgecolor="lightgray", linewidth=0.4),
    )

    ax.set_xlabel("[DCPIP] (µM)")
    ax.set_ylabel("A$_{600}$ (blank-subtracted)")
    ax.set_title(title or "DCPIP standard curve")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.legend(loc="upper left", frameon=False)

    fig.tight_layout()
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf, format="pdf")
    fig.savefig(out_pdf.with_suffix(".png"), format="png", dpi=600)
    plt.close(fig)


# =============================================================================
# Text summary
# =============================================================================

def _print_summary(result: StandardCurveResult) -> None:
    print(f"start conc:          {result.start_conc_uM:g} µM (row A)")
    print(f"dilution:            {result.dilution_ratio}  "
          f"(per-step factor {result.dilution_factor:.6g})")
    print(f"primary pH for eps:  {result.primary_ph:g}")
    for ph in sorted(result.per_ph):
        fit = result.per_ph[ph]
        print(f"\n[pH {ph:g}]  n = {fit.n_points}")
        print(f"  eps_app (AU/µM):   {fit.eps_app:.6g}  (SE {fit.eps_app_stderr:.3g})")
        print(f"  R^2:               {fit.r_squared:.5f}")
        print(f"  blanks (row P):    "
              + ", ".join(f"c{c}={v:.3f}" for c, v in fit.blanks.items()))


if __name__ == "__main__":
    main()
