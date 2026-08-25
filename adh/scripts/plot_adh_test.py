"""Plot the 260722 ADH kinetic test plate and identify candidate protein-titration rows.

The BioTek export contains 384 kinetic traces. This script finds rows containing
substantial DCPIP absorbance loss, plots all active wells, and highlights the
strongest active row as the likely protein-titration row. The inference is
reported as a candidate only: confirm the well map before assigning protein
concentrations.

Usage
-----
    python plot_adh_test.py test/260722_adh_test.xlsx
    python plot_adh_test.py test/260722_adh_test.xlsx --row B
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

from kcat_plate import PLATE_COLS, PLATE_ROWS, parse_kinetic_plate


DEFAULT_OUTPUT_SUBDIR = "outputs"
ACTIVE_DROP_AU = 0.05

NATURE_RC = {
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
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
}
NATURE_PANEL_INCHES = (7.08, 5.20)


def _valid_trace(trace: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    indices = np.flatnonzero(np.isfinite(trace))
    return indices, trace[indices]


def _absorbance_drop(trace: np.ndarray) -> float:
    _, values = _valid_trace(trace)
    if len(values) < 6:
        return np.nan
    n = max(3, min(len(values) // 10, 20))
    return float(np.mean(values[:n]) - np.mean(values[-n:]))


def find_active_rows(traces, threshold: float = ACTIVE_DROP_AU) -> list[tuple[str, int, float]]:
    """Return rows with active wells: row, active-well count, mean absorbance drop."""
    summary = []
    for row in PLATE_ROWS:
        drops = np.asarray([_absorbance_drop(traces.traces[f"{row}{col}"]) for col in PLATE_COLS])
        active = drops[np.isfinite(drops) & (drops >= threshold)]
        if active.size:
            summary.append((row, int(active.size), float(np.mean(active))))
    return sorted(summary, key=lambda item: (item[1], item[2]), reverse=True)


def _time_minutes(traces) -> np.ndarray:
    """Build a monotonic time axis from the positive intervals in the export."""
    differences = np.diff(traces.times_s)
    positive = differences[np.isfinite(differences) & (differences > 0)]
    interval_s = float(np.median(positive)) if positive.size else 20.0
    return np.arange(len(traces.times_s), dtype=float) * interval_s / 60.0


def _linear_fit(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    """Return slope, intercept, and $R^2$ for an ordinary least-squares line."""
    if len(x) < 3:
        return np.nan, np.nan, np.nan
    slope, intercept = np.polyfit(x, y, 1)
    fitted = slope * x + intercept
    residual = float(np.sum((y - fitted) ** 2))
    total = float(np.sum((y - np.mean(y)) ** 2))
    return float(slope), float(intercept), 1 - residual / total if total > 0 else np.nan


def find_initial_linear_region(
    time_min: np.ndarray,
    values: np.ndarray,
    *,
    window_min: float = 5.0,
    min_r2: float = 0.90,
) -> tuple[float, float, float, float, float, int]:
    """Apply the same initial-window search used by ``kcat_plate.py``.

    The fit always begins at the first measurement and is limited to the shorter
    of ``window_min`` and half of the recorded trace. The longest window that
    meets the $R^2$ threshold is retained.
    """
    if len(values) < 5:
        return np.nan, np.nan, np.nan, np.nan, np.nan, 0
    fallback = None
    best = None
    for end in range(4, len(values) + 1):
        if time_min[end - 1] > time_min[0] + min(window_min, 0.5 * (time_min[-1] - time_min[0])):
            break
        slope, intercept, r_squared = _linear_fit(time_min[:end], values[:end])
        candidate = (time_min[0], time_min[end - 1], slope, intercept, r_squared, end)
        fallback = fallback or candidate
        if np.isfinite(r_squared) and r_squared >= min_r2:
            best = candidate
    return best or fallback


def plot_active_traces(
    traces,
    selected_row: str,
    out_path: Path,
    *,
    start_protein_uM: float = 100.0,
    dilution_factor: float = 0.5,
    fit_window_min: float = 5.0,
    min_r2: float = 0.90,
    title: str = "",
) -> list[dict[str, float | str]]:
    """Plot titration traces and their initial-fit quality statistics."""
    plt.rcParams.update(NATURE_RC)
    time_min = _time_minutes(traces)
    fig, axes = plt.subplots(1, 2, figsize=NATURE_PANEL_INCHES, layout="constrained")
    titration_ax, quality_ax = axes

    candidate_wells = [
        f"{selected_row}{col}" for col in PLATE_COLS
        if _absorbance_drop(traces.traces[f"{selected_row}{col}"]) >= ACTIVE_DROP_AU
    ]
    colours = plt.get_cmap("viridis")(np.linspace(0.15, 0.90, max(len(candidate_wells), 1)))
    linear_regions: list[dict[str, float | str]] = []
    for colour, well in zip(colours, candidate_wells):
        concentration = start_protein_uM * dilution_factor ** (int(well[1:]) - 1)
        indices, values = _valid_trace(traces.traces[well])
        well_time = time_min[indices]
        t0, t1, slope, intercept, r_squared, n_points = find_initial_linear_region(
            well_time, values, window_min=fit_window_min, min_r2=min_r2,
        )
        titration_ax.plot(time_min[indices], values, color=colour, lw=0.7,
                          label=f"{well}: {concentration:g} µM")
        if n_points:
            region_time = well_time[:n_points]
            region_values = values[:n_points]
            pad = max(np.ptp(region_values) * 0.08, 0.003)
            titration_ax.add_patch(Rectangle(
                (t0, float(np.min(region_values)) - pad), t1 - t0,
                float(np.ptp(region_values)) + 2 * pad,
                fill=False, edgecolor=colour, linewidth=0.7, linestyle="--",
            ))
            titration_ax.plot(region_time, slope * region_time + intercept,
                              color=colour, lw=0.7)
        linear_regions.append({
            "well": well, "protein_uM": concentration, "start_min": t0, "end_min": t1,
            "slope_A600_per_min": slope, "r_squared": r_squared, "n_points": n_points,
        })
    titration_ax.set_xlabel("time (min)")
    titration_ax.set_ylabel("A$_{600}$")
    titration_ax.set_title(f"Protein titration: row {selected_row}")
    titration_ax.legend(frameon=False, ncol=2, loc="upper right", fontsize=5)
    titration_ax.text(0.02, 0.02, "Dashed boxes: fitted initial regions",
                      transform=titration_ax.transAxes, ha="left", va="bottom", fontsize=5)

    concentrations = np.asarray([float(region["protein_uM"]) for region in linear_regions])
    r_squared = np.asarray([float(region["r_squared"]) for region in linear_regions])
    quality_ax.scatter(concentrations, r_squared, c=colours, s=16, edgecolors="none", zorder=3)
    for region in linear_regions:
        concentration_label = f"{float(region['protein_uM']):g} µM"
        quality_ax.annotate(
            f"{region['well']}\n{concentration_label}",
            (float(region["protein_uM"]), float(region["r_squared"])),
            xytext=(3, 3), textcoords="offset points", fontsize=5,
        )
    quality_ax.axhline(min_r2, color="0.3", lw=0.5, ls="--", label=f"R² threshold ({min_r2:.2f})")
    quality_ax.set_xscale("log")
    quality_ax.set_xlabel("protein concentration (µM; assumed)")
    quality_ax.set_ylabel("initial-fit R²")
    quality_ax.set_title("Linear-region fit quality")
    quality_ax.set_ylim(min(0.85, float(np.nanmin(r_squared)) - 0.01), 1.002)
    quality_ax.legend(frameon=False, loc="lower left", fontsize=5)

    if title:
        fig.suptitle(title, fontsize=7)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, format="pdf")
    fig.savefig(out_path.with_suffix(".png"), format="png", dpi=600)
    plt.close(fig)
    return linear_regions


def write_linear_regions(regions: list[dict[str, float | str]], out_path: Path) -> None:
    """Write fitted initial linear windows for each active well in the titration row."""
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "well", "protein_uM", "start_min", "end_min", "slope_A600_per_min", "r_squared", "n_points",
        ])
        writer.writeheader()
        writer.writerows(regions)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("xlsx", type=Path, help="BioTek kinetic .xlsx export.")
    parser.add_argument("--row", choices=PLATE_ROWS, default=None,
                        help="Known protein-titration row; defaults to strongest active row.")
    parser.add_argument("--start-protein-uM", type=float, default=100.0,
                        help="Protein concentration in column 1 of the titration row (default: 100 µM).")
    parser.add_argument("--fit-window-min", type=float, default=5.0,
                        help="Maximum initial fitting window in minutes (default: 5; matches kcat_plate.py).")
    parser.add_argument("--min-r2", type=float, default=0.90,
                        help="Minimum R² for an initial linear region (default: 0.90; matches kcat_plate.py).")
    parser.add_argument("--output", type=Path, default=None,
                        help="Output PDF (default: sibling outputs directory).")
    args = parser.parse_args()

    traces = parse_kinetic_plate(args.xlsx)
    active_rows = find_active_rows(traces)
    if not active_rows:
        raise SystemExit(f"No rows passed the {ACTIVE_DROP_AU:g} AU activity threshold.")
    selected_row = args.row or active_rows[0][0]
    output = args.output or (args.xlsx.resolve().parent.parent / DEFAULT_OUTPUT_SUBDIR /
                             f"{args.xlsx.stem}_active_traces.pdf")
    regions = plot_active_traces(traces, selected_row, output,
                                 start_protein_uM=args.start_protein_uM,
                                 fit_window_min=args.fit_window_min,
                                 min_r2=args.min_r2, title=args.xlsx.stem)
    regions_csv = output.with_name(f"{output.stem}_linear_regions.csv")
    write_linear_regions(regions, regions_csv)

    print("Active rows (ranked by active-well count, then mean absorbance drop):")
    for row, count, mean_drop in active_rows:
        print(f"  {row}: {count} wells; mean ΔA600 = {mean_drop:.3f}")
    print(f"Candidate protein-titration row: {selected_row} "
          f"(assumed {args.start_protein_uM:g} µM in col. 1; 1:1 serial dilution)")
    for region in regions:
        print(f"  {region['well']}: {region['start_min']:.1f}–{region['end_min']:.1f} min; "
              f"R² = {region['r_squared']:.3f}")
    print(f"Wrote: {output}")
    print(f"Wrote: {regions_csv}")


if __name__ == "__main__":
    main()
