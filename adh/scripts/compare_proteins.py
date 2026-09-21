"""Compare two proteins across (pH, metal, alcohol) conditions from kcat_plate CSVs.

Pipeline
--------
    1. Load one per-well k_cat CSV per (protein, pH) condition. Accepts either
       ``_wells.csv`` (single plate, column ``kcat_s``) or ``_avg_wells.csv``
       (averaged replicates, column ``kcat_mean_s``).
    2. Merge everything into a long-form table:
       protein, pH, well, row, col, metal, substrate, kcat, kcat_std,
    mw_da, aromatic.
    3. Dump summary CSVs marginalised by pH / metal / alcohol.
    4. Render Nature-style comparison figures:
        * pH effect (per protein, aggregated over metals + alcohols)
        * Metal effect (per protein, faceted by pH)
        * Alcohol size effect (k_cat vs molecular weight, linear fit)
        * Direct protein-vs-protein scatter (when exactly two proteins present)

Usage
-----
    # Manifest form (preferred):
    python compare_proteins.py --manifest manifest.json --outdir outputs/

    # CLI form (repeat --condition per file):
    python compare_proteins.py \\
        --condition ADH1 7  path/to/ADH1_pH7_avg_wells.csv \\
        --condition ADH1 8  path/to/ADH1_pH8_avg_wells.csv \\
        --condition PedE 7  path/to/PedE_pH7_avg_wells.csv \\
        --condition PedE 8  path/to/PedE_pH8_avg_wells.csv \\
        --outdir outputs/

Manifest JSON
-------------
    {
      "conditions": [
        {"protein": "ADH1", "pH": 7.0, "csv": "outputs/ADH1_pH7_avg_wells.csv"},
        {"protein": "ADH1", "pH": 8.0, "csv": "outputs/ADH1_pH8_avg_wells.csv"},
        {"protein": "PedE", "pH": 7.0, "csv": "outputs/PedE_pH7_avg_wells.csv"},
        {"protein": "PedE", "pH": 8.0, "csv": "outputs/PedE_pH8_avg_wells.csv"}
      ]
    }
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


# =============================================================================
# Editable configuration
# =============================================================================

# Metals for the screen plate, ordered by atomic number (must match kcat_plate.py).
METALS: list[str] = [
    "Al", "Ca", "Sc", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Y",
    "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho",
    "Er", "Tm", "Yb", "Lu",
]

METAL_MASS_DA = {
    "Al": 26.9815, "Ca": 40.078, "Sc": 44.9559, "Mn": 54.9380,
    "Fe": 55.845, "Co": 58.9332, "Ni": 58.6934, "Cu": 63.546,
    "Zn": 65.38, "Y": 88.9058, "La": 138.9055, "Ce": 140.116,
    "Pr": 140.9077, "Nd": 144.242, "Sm": 150.36, "Eu": 151.964,
    "Gd": 157.25, "Tb": 158.9254, "Dy": 162.500, "Ho": 164.9303,
    "Er": 167.259, "Tm": 168.9342, "Yb": 173.045, "Lu": 174.9668,
}

# Molecular weights are for the neutral substrates, in daltons (Da).
ALCOHOL_PROPS: dict[str, dict] = {
    "methanol":                {"mw_da": 32.04, "aromatic": False},
    "ethanol":                 {"mw_da": 46.07, "aromatic": False},
    "1-propanol":              {"mw_da": 60.10, "aromatic": False},
    "1-butanol":               {"mw_da": 74.12, "aromatic": False},
    "3-methyl-1-butanol":      {"mw_da": 88.15, "aromatic": False},
    "1-pentanol":              {"mw_da": 88.15, "aromatic": False},
    "2-methyl-1-butanol":      {"mw_da": 88.15, "aromatic": False},
    "1-hexanol":               {"mw_da": 102.17, "aromatic": False},
    "1,5-pentanediol":         {"mw_da": 104.15, "aromatic": False},
    "2-phenylethanol":         {"mw_da": 122.16, "aromatic": True},
    "4-hydroxybenzyl alcohol": {"mw_da": 124.14, "aromatic": True},
    "vanillin":                {"mw_da": 152.15, "aromatic": True},
    "protocatechuic acid":     {"mw_da": 154.12, "aromatic": True},
    "vanillyl alcohol":        {"mw_da": 154.16, "aromatic": True},
    "vanillic acid":           {"mw_da": 168.15, "aromatic": True},
}
ALCOHOLS: list[str] = list(ALCOHOL_PROPS.keys())

# Control markers (must match kcat_plate.py) that we exclude from analysis.
CONTROL_SUBSTRATES = {"water", "tcep", "no_substrate", ""}
CONTROL_METALS = {"no_metal", ""}

DEFAULT_OUTPUT_SUBDIR: str = "outputs"

# Nature-panel style (see .github/instructions/figure-style.instructions.md).
NATURE_RC = {
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"],
    "font.size": 7,
    "axes.labelsize": 7,
    "axes.titlesize": 7,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "legend.title_fontsize": 7,
    "figure.titlesize": 7,
    "figure.labelsize": 7,
    "mathtext.fontset": "custom",
    "mathtext.rm": "Helvetica",
    "mathtext.it": "Helvetica:italic",
    "mathtext.bf": "Helvetica:bold",
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size": 2.0,
    "ytick.major.size": 2.0,
    "axes.spines.top": False,
    "axes.spines.right": False,
}
NATURE_PANEL_INCHES = (7.08, 5.20)

# Two-protein palette; extra proteins fall back to matplotlib tab10.
PROTEIN_COLORS: dict[str, str] = {}
_TAB10 = plt.get_cmap("tab10").colors


# =============================================================================
# Data classes
# =============================================================================

@dataclass
class WellRecord:
    protein: str
    pH: float
    well: str
    row: str
    col: int
    metal: str
    substrate: str
    kcat: float
    kcat_std: float          # NaN if unknown (single-plate CSV)
    mw_da: float
    aromatic: bool


@dataclass(frozen=True)
class MetalProperties:
    charge: int
    ionic_radius_pm: float
    coordination_number: int = 6
    spin: str = ""
    source: str = "Shannon 1976; doi:10.1107/S0567739476001551"
    assumption: str = "Assumed salt charge and CN; verify before interpreting real data"


DEFAULT_METAL_PROPERTIES = {
    "Al": MetalProperties(3, 53.5),
    "Ca": MetalProperties(2, 100.0),
    "Sc": MetalProperties(3, 74.5),
    "Mn": MetalProperties(2, 83.0, spin="high"),
    "Fe": MetalProperties(3, 64.5, spin="high"),
    "Co": MetalProperties(2, 74.5, spin="high"),
    "Ni": MetalProperties(2, 69.0),
    "Cu": MetalProperties(2, 73.0),
    "Zn": MetalProperties(2, 74.0),
    "Y": MetalProperties(3, 90.0),
    "La": MetalProperties(3, 103.2),
    "Ce": MetalProperties(3, 101.0),
    "Pr": MetalProperties(3, 99.0),
    "Nd": MetalProperties(3, 98.3),
    "Sm": MetalProperties(3, 95.8),
    "Eu": MetalProperties(3, 94.7),
    "Gd": MetalProperties(3, 93.8),
    "Tb": MetalProperties(3, 92.3),
    "Dy": MetalProperties(3, 91.2),
    "Ho": MetalProperties(3, 90.1),
    "Er": MetalProperties(3, 89.0),
    "Tm": MetalProperties(3, 88.0),
    "Yb": MetalProperties(3, 86.8),
    "Lu": MetalProperties(3, 86.1),
}


def _load_metal_properties(path: Path | None) -> dict[str, MetalProperties]:
    """Read optional overrides of the assumed CN-6 Shannon effective radii."""
    properties = DEFAULT_METAL_PROPERTIES.copy()
    if path is None:
        return properties
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"metal", "charge", "ionic_radius_pm", "coordination_number", "spin", "source"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"{path}: expected columns {sorted(required)}")
        seen = set()
        for row in reader:
            metal = row["metal"].strip()
            charge = int(row["charge"])
            radius = float(row["ionic_radius_pm"])
            coordination = int(row["coordination_number"])
            if not metal or metal in seen or charge <= 0 or coordination <= 0 or not np.isfinite(radius) or radius <= 0:
                raise ValueError(f"{path}: invalid or duplicate metal descriptor: {row}")
            seen.add(metal)
            properties[metal] = MetalProperties(
                charge, radius, coordination, row["spin"].strip(), row["source"].strip(),
                row.get("assumption", "User-supplied descriptors"),
            )
    return properties


# =============================================================================
# CLI / main
# =============================================================================

def _build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--manifest", type=Path, default=None,
                   help="JSON manifest {conditions: [{protein, pH, csv}, ...]}.")
    p.add_argument("--condition", nargs=3, action="append", metavar=("PROTEIN", "PH", "CSV"),
                   help="Add a (protein, pH, csv) condition. May be repeated. "
                        "Ignored if --manifest is given.")
    p.add_argument("--outdir", type=Path, default=None,
                   help="Directory for figures + CSVs (default: manifest's ../outputs).")
    p.add_argument("--label", type=str, default="compare",
                   help="Filename prefix for outputs (default: 'compare').")
    p.add_argument("--kcat-upper", type=float, default=None,
                   help="Upper bound for y-axes / colour scales (s^-1). Default automatic.")
    p.add_argument("--include-aromatic", action="store_true",
                   help="Include aromatic substrates in the alcohol-size regression "
                        "(they are chemically distinct and off by default).")
    p.add_argument("--metal-properties", type=Path, default=None,
                   help="Override assumed CN-6 metal descriptors: CSV columns metal, charge, "
                        "ionic_radius_pm, coordination_number, spin, source. "
                        "Use the exported metal_properties.csv as a template.")
    return p


def main() -> None:
    parser = _build_argparser()
    args = parser.parse_args()
    if args.kcat_upper is not None and (not np.isfinite(args.kcat_upper) or args.kcat_upper <= 0):
        parser.error("--kcat-upper must be finite and greater than 0.")

    conditions = _resolve_conditions(args, parser)
    outdir = args.outdir or _default_outdir(conditions[0][2])
    outdir.mkdir(parents=True, exist_ok=True)

    records: list[WellRecord] = []
    for protein, pH, csv_path in conditions:
        records.extend(_load_wells_csv(csv_path, protein=protein, pH=pH))
    if not records:
        raise SystemExit("No usable enzyme wells found across the supplied CSVs.")

    proteins = sorted({r.protein for r in records})
    _assign_protein_colors(proteins)
    pH_values = sorted({r.pH for r in records})
    metal_properties = _load_metal_properties(args.metal_properties)
    unknown_metals = {record.metal for record in records} - metal_properties.keys()
    unknown_substrates = {record.substrate for record in records if not np.isfinite(record.mw_da)}
    if unknown_metals or unknown_substrates:
        parser.error(f"Missing descriptors: metals={sorted(unknown_metals)}, substrates={sorted(unknown_substrates)}")
    if any(metal_properties[record.metal].assumption.startswith("Assumed") for record in records):
        print("NOTE: metal charge/CN assumptions are unverified; see exported metal_properties.csv.")

    stem = args.label
    _dump_long_form(records, outdir / f"{stem}_wells_long.csv")
    _dump_summary_by_ph(records, outdir / f"{stem}_summary_by_ph.csv")
    _dump_summary_by_metal(records, outdir / f"{stem}_summary_by_metal.csv")
    _dump_summary_by_alcohol(records, outdir / f"{stem}_summary_by_alcohol.csv")
    _dump_factor_data(records, metal_properties, outdir, stem)

    plot_ph_effect(records, proteins, pH_values,
                   outdir / f"{stem}_ph_effect.pdf",
                   vmax=args.kcat_upper)
    plot_metal_effect(records, proteins, pH_values,
                      outdir / f"{stem}_metal_effect.pdf",
                      vmax=args.kcat_upper)
    plot_alcohol_size_effect(records, proteins, pH_values,
                             outdir / f"{stem}_alcohol_size.pdf",
                             include_aromatic=args.include_aromatic,
                             vmax=args.kcat_upper)
    if len(proteins) == 2:
        plot_protein_vs_protein(records, proteins, pH_values,
                                outdir / f"{stem}_protein_scatter.pdf",
                                vmax=args.kcat_upper)
    for factor in ("ph", "metal_radius", "metal_charge", "metal_id", "alcohol_mw"):
        plot_factor_scatter(records, proteins, pH_values, metal_properties, factor,
                            outdir / f"{stem}_{factor}_scatter.pdf", vmax=args.kcat_upper)
    plot_combined_factors(records, proteins, pH_values, metal_properties,
                          outdir / f"{stem}_combined_factors.pdf", vmax=args.kcat_upper)
    plot_combined_factors(records, proteins, pH_values, metal_properties,
                          outdir / f"{stem}_combined_factors_metal_mw.pdf",
                          vmax=args.kcat_upper, metal_order="mass")
    for metal_axis in ("radius", "mass"):
        plot_factors_3d(records, proteins, pH_values, metal_properties,
                        outdir / f"{stem}_3d_metal_{metal_axis}.pdf",
                        metal_axis=metal_axis, vmax=args.kcat_upper)

    print(f"Wrote {len(records)} records to {outdir}")
    print(f"  proteins: {proteins}")
    print(f"  pH:       {pH_values}")


# =============================================================================
# Input parsing
# =============================================================================

def _resolve_conditions(args, parser: argparse.ArgumentParser) -> list[tuple[str, float, Path]]:
    if args.manifest is not None:
        data = json.loads(args.manifest.read_text())
        raw = data.get("conditions", data)
        conditions: list[tuple[str, float, Path]] = []
        for entry in raw:
            try:
                protein = str(entry["protein"])
                pH = float(entry["pH"])
                csv_path = Path(entry["csv"])
            except (KeyError, TypeError, ValueError) as exc:
                parser.error(f"Malformed manifest entry {entry!r}: {exc}")
            if not csv_path.is_absolute():
                csv_path = (args.manifest.parent / csv_path).resolve()
            conditions.append((protein, pH, csv_path))
        if not conditions:
            parser.error("Manifest contained zero conditions.")
        return conditions
    if not args.condition:
        parser.error("Provide either --manifest or one or more --condition PROTEIN PH CSV.")
    conditions = []
    for protein, pH_str, csv_str in args.condition:
        try:
            pH = float(pH_str)
        except ValueError:
            parser.error(f"pH must be numeric, got {pH_str!r}.")
        conditions.append((protein, pH, Path(csv_str).resolve()))
    return conditions


def _default_outdir(sample_csv: Path) -> Path:
    return sample_csv.resolve().parent.parent / DEFAULT_OUTPUT_SUBDIR


def _load_wells_csv(csv_path: Path, *, protein: str, pH: float) -> list[WellRecord]:
    """Read a wells CSV (single-plate or averaged) and keep only enzyme wells."""
    with csv_path.open() as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
    if not rows:
        return []
    kcat_col = "kcat_mean_s" if "kcat_mean_s" in rows[0] else "kcat_s"
    std_col = "kcat_std_s" if "kcat_std_s" in rows[0] else None
    if kcat_col not in rows[0]:
        raise ValueError(f"{csv_path}: missing k_cat column (expected 'kcat_mean_s' or 'kcat_s').")
    records: list[WellRecord] = []
    for row in rows:
        substrate = row.get("substrate", "").strip()
        metal = row.get("metal", "").strip()
        if substrate in CONTROL_SUBSTRATES or metal in CONTROL_METALS:
            continue
        raw = row.get(kcat_col, "")
        if raw in ("", "nan", "NaN"):
            continue
        try:
            kcat = float(raw)
        except ValueError:
            continue
        if not np.isfinite(kcat):
            continue
        std_val = float("nan")
        if std_col and row.get(std_col, "") not in ("", "nan", "NaN"):
            try:
                std_val = float(row[std_col])
            except ValueError:
                std_val = float("nan")
        props = ALCOHOL_PROPS.get(substrate, {"mw_da": float("nan"), "aromatic": False})
        records.append(WellRecord(
            protein=protein,
            pH=pH,
            well=row.get("well", ""),
            row=row.get("row", ""),
            col=int(row.get("col", 0) or 0),
            metal=metal,
            substrate=substrate,
            kcat=kcat,
            kcat_std=std_val,
            mw_da=float(props["mw_da"]),
            aromatic=bool(props["aromatic"]),
        ))
    return records


# =============================================================================
# Palette helpers
# =============================================================================

def _assign_protein_colors(proteins: list[str]) -> None:
    for i, name in enumerate(proteins):
        if name not in PROTEIN_COLORS:
            PROTEIN_COLORS[name] = _TAB10[i % len(_TAB10)]


def _ph_alpha(pH: float, pH_values: list[float]) -> float:
    """Encode pH as an alpha in [0.35, 1.0]: low pH lighter, high pH darker."""
    if len(pH_values) <= 1:
        return 1.0
    lo, hi = min(pH_values), max(pH_values)
    if hi == lo:
        return 1.0
    return 0.35 + 0.65 * (pH - lo) / (hi - lo)


# =============================================================================
# Aggregation helpers
# =============================================================================

def _group_stats(values: list[float]) -> tuple[float, float, int]:
    arr = np.asarray([v for v in values if np.isfinite(v)], dtype=float)
    if arr.size == 0:
        return float("nan"), float("nan"), 0
    return float(np.mean(arr)), float(np.std(arr, ddof=1) if arr.size > 1 else 0.0), int(arr.size)


def _filter(records: list[WellRecord], **filters) -> list[WellRecord]:
    def match(rec: WellRecord) -> bool:
        for key, val in filters.items():
            if getattr(rec, key) != val:
                return False
        return True
    return [r for r in records if match(r)]


def _sem(std: float, n: int) -> float:
    if n <= 1 or not np.isfinite(std):
        return float("nan")
    return std / np.sqrt(n)


# =============================================================================
# CSV dumps
# =============================================================================

def _dump_long_form(records: list[WellRecord], out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["protein", "pH", "well", "row", "col", "metal", "substrate",
                    "kcat_s", "kcat_std_s", "mw_da", "aromatic"])
        for r in records:
            w.writerow([r.protein, r.pH, r.well, r.row, r.col, r.metal, r.substrate,
                        r.kcat,
                        "" if not np.isfinite(r.kcat_std) else r.kcat_std,
                        "" if not np.isfinite(r.mw_da) else r.mw_da,
                        int(r.aromatic)])


def _dump_summary_by_ph(records: list[WellRecord], out_csv: Path) -> None:
    proteins = sorted({r.protein for r in records})
    pH_values = sorted({r.pH for r in records})
    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["protein", "pH", "mean_kcat_s", "std_kcat_s", "sem_kcat_s", "n_wells"])
        for protein in proteins:
            for pH in pH_values:
                subset = _filter(records, protein=protein, pH=pH)
                mean, std, n = _group_stats([r.kcat for r in subset])
                w.writerow([protein, pH, mean, std, _sem(std, n), n])


def _dump_summary_by_metal(records: list[WellRecord], out_csv: Path) -> None:
    proteins = sorted({r.protein for r in records})
    pH_values = sorted({r.pH for r in records})
    metals_present = [m for m in METALS if any(r.metal == m for r in records)]
    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["protein", "pH", "metal", "mean_kcat_s", "std_kcat_s", "sem_kcat_s", "n_wells"])
        for protein in proteins:
            for pH in pH_values:
                for metal in metals_present:
                    subset = _filter(records, protein=protein, pH=pH, metal=metal)
                    mean, std, n = _group_stats([r.kcat for r in subset])
                    w.writerow([protein, pH, metal, mean, std, _sem(std, n), n])


def _dump_summary_by_alcohol(records: list[WellRecord], out_csv: Path) -> None:
    proteins = sorted({r.protein for r in records})
    pH_values = sorted({r.pH for r in records})
    alcohols_present = [a for a in ALCOHOLS if any(r.substrate == a for r in records)]
    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["protein", "pH", "alcohol", "mw_da", "aromatic",
                    "mean_kcat_s", "std_kcat_s", "sem_kcat_s", "n_wells"])
        for protein in proteins:
            for pH in pH_values:
                for alcohol in alcohols_present:
                    subset = _filter(records, protein=protein, pH=pH, substrate=alcohol)
                    mean, std, n = _group_stats([r.kcat for r in subset])
                    props = ALCOHOL_PROPS[alcohol]
                    w.writerow([protein, pH, alcohol, props["mw_da"],
                                int(props["aromatic"]), mean, std, _sem(std, n), n])


# =============================================================================
# Plots
# =============================================================================

def _condition_means(records: list[WellRecord]) -> list[WellRecord]:
    """Give each measured protein/pH/metal/substrate combination one plotting point."""
    grouped: dict[tuple[str, float, str, str], list[WellRecord]] = {}
    for record in records:
        key = (record.protein, record.pH, record.metal, record.substrate)
        grouped.setdefault(key, []).append(record)
    return [WellRecord(
        protein=key[0], pH=key[1], well="", row="", col=0, metal=key[2], substrate=key[3],
        kcat=float(np.mean([record.kcat for record in group])), kcat_std=float("nan"),
        mw_da=group[0].mw_da, aromatic=group[0].aromatic,
    ) for key, group in sorted(grouped.items())]


def _dump_factor_data(
    records: list[WellRecord], properties: dict[str, MetalProperties], outdir: Path, stem: str,
) -> None:
    with (outdir / f"{stem}_metal_properties.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metal", "charge", "ionic_radius_pm", "coordination_number", "spin", "source", "assumption"])
        for metal in sorted({record.metal for record in records}):
            prop = properties[metal]
            writer.writerow([metal, prop.charge, prop.ionic_radius_pm, prop.coordination_number,
                             prop.spin, prop.source, prop.assumption])
    with (outdir / f"{stem}_factor_conditions.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["protein", "pH", "metal", "charge", "ionic_radius_pm", "coordination_number",
                         "substrate", "mw_da", "aromatic", "mean_kcat_s"])
        for record in _condition_means(records):
            prop = properties[record.metal]
            writer.writerow([record.protein, record.pH, record.metal, prop.charge, prop.ionic_radius_pm,
                             prop.coordination_number, record.substrate, record.mw_da,
                             int(record.aromatic), record.kcat])


def plot_factor_scatter(
    records: list[WellRecord], proteins: list[str], pH_values: list[float],
    properties: dict[str, MetalProperties], factor: str, out_path: Path,
    *, vmax: float | None = None,
) -> None:
    """Plot condition means without pooling other factors into SEM error bars."""
    plt.rcParams.update(NATURE_RC)
    points = _condition_means(records)
    metals = [metal for metal in METALS if any(record.metal == metal for record in points)]
    metals += sorted({record.metal for record in points} - set(metals))
    labels = {
        "ph": "pH", "metal_radius": "metal ionic radius (pm)",
        "metal_charge": "metal charge", "metal_id": "metal",
        "alcohol_mw": "substrate MW (Da)",
    }
    if factor not in labels:
        raise ValueError(f"Unknown factor: {factor}")
    facets = [None] if factor == "ph" else pH_values
    fig, axes = plt.subplots(len(facets), 1, sharex=True, sharey=True, squeeze=False,
                             figsize=(7.08, min(2.1 * len(facets) + 0.35, 6.69)))
    markers = ["o", "^", "s", "D"]
    values = {
        "ph": lambda record: record.pH,
        "metal_radius": lambda record: properties[record.metal].ionic_radius_pm,
        "metal_charge": lambda record: properties[record.metal].charge,
        "metal_id": lambda record: metals.index(record.metal),
        "alcohol_mw": lambda record: record.mw_da,
    }
    for axis, ph in zip(axes.flat, facets):
        for protein_index, protein in enumerate(proteins):
            subset = [record for record in points if record.protein == protein and (ph is None or record.pH == ph)]
            axis.scatter([values[factor](record) for record in subset], [record.kcat for record in subset],
                         s=8, marker=markers[protein_index % len(markers)], alpha=0.55,
                         color=PROTEIN_COLORS[protein], linewidths=0, label=protein)
        if ph is not None:
            axis.set_title(f"pH {ph:g}", loc="left")
        if vmax is not None:
            axis.set_ylim(top=vmax)
        axis.margins(x=0.03)
    bottom_axis = axes[-1, 0]
    if factor == "ph":
        bottom_axis.set_xticks(pH_values, labels=[f"{ph:g}" for ph in pH_values])
        bottom_axis.set_xticks([], minor=True)
    elif factor == "metal_charge":
        charges = sorted({properties[record.metal].charge for record in points})
        bottom_axis.set_xticks(charges, labels=[f"+{charge}" for charge in charges])
        bottom_axis.set_xticks([], minor=True)
    elif factor == "metal_id":
        bottom_axis.set_xticks(range(len(metals)), labels=metals, rotation=90)
    bottom_axis.set_xlabel(labels[factor])
    fig.supylabel("k$_{cat}$ (s$^{-1}$)")
    axes[0, 0].legend(frameon=False, loc="upper right", ncol=len(proteins))
    fig.tight_layout()
    _finalize(fig, out_path)


def plot_combined_factors(
    records: list[WellRecord], proteins: list[str], pH_values: list[float],
    properties: dict[str, MetalProperties], out_path: Path, *, vmax: float | None = None,
    metal_order: str = "radius",
) -> None:
    """Maps ordered by metal radius or atomic mass retain separate chemical identities."""
    plt.rcParams.update(NATURE_RC)
    points = _condition_means(records)
    if metal_order == "mass":
        descriptors = METAL_MASS_DA
        metal_label = "Metal (atomic mass, Da); ordered by mass"
    elif metal_order == "radius":
        descriptors = {metal: prop.ionic_radius_pm for metal, prop in properties.items()}
        metal_label = "Metal (ionic radius, pm); ordered by radius"
    else:
        raise ValueError(f"Unknown metal order: {metal_order}")
    missing = {record.metal for record in points} - descriptors.keys()
    if missing:
        raise ValueError(f"Missing metal {metal_order} descriptors: {sorted(missing)}")
    metals = sorted({record.metal for record in points}, key=lambda metal: (descriptors[metal], metal))
    molecular_weights = {record.substrate: record.mw_da for record in points}
    substrates = sorted(molecular_weights, key=lambda substrate: (molecular_weights[substrate], substrate))
    lookup = {(record.protein, record.pH, record.metal, record.substrate): record.kcat for record in points}
    values = np.array(list(lookup.values()))
    lower = min(0.0, float(values.min()))
    upper = vmax if vmax is not None else max(float(values.max()), 0.001)
    cmap = plt.get_cmap("gray_r").copy()
    cmap.set_bad("#ff5252")
    fig, axes = plt.subplots(len(pH_values), len(proteins), squeeze=False, sharex=True, sharey=True,
                             figsize=(7.08, 6.69))
    for row_index, ph in enumerate(pH_values):
        for col_index, protein in enumerate(proteins):
            axis = axes[row_index, col_index]
            grid = np.array([[lookup.get((protein, ph, metal, substrate), np.nan)
                              for metal in metals] for substrate in substrates])
            image = axis.imshow(grid, cmap=cmap, aspect="auto", interpolation="nearest", vmin=lower, vmax=upper)
            axis.set_title(f"{protein}, pH {ph:g}")
            metal_labels = [f"{metal} ({descriptors[metal]:.2f})" if metal_order == "mass"
                            else f"{metal} ({descriptors[metal]:g})" for metal in metals]
            axis.set_xticks(range(len(metals)), labels=metal_labels, rotation=90)
            axis.set_yticks(range(len(substrates)), labels=[f"{substrate} ({molecular_weights[substrate]:.2f})"
                                                          for substrate in substrates])
            axis.tick_params(length=0, pad=2)
    fig.subplots_adjust(left=0.30, right=0.985, bottom=0.15, top=0.90, wspace=0.08, hspace=0.24)
    color_axis = fig.add_axes((0.38, 0.96, 0.45, 0.012))
    colorbar = fig.colorbar(image, cax=color_axis, orientation="horizontal")
    colorbar.set_label("k$_{cat}$ (s$^{-1}$)", labelpad=2)
    colorbar.ax.xaxis.set_label_position("top")
    colorbar.outline.set_linewidth(0.4)
    colorbar.ax.tick_params(labelsize=7, length=2, pad=1)
    fig.text(0.64, 0.015, metal_label, ha="center", fontsize=7)
    fig.text(0.01, 0.52, "Substrate (MW, Da); ordered by MW", rotation=90, va="center", fontsize=7)
    _finalize(fig, out_path)


def plot_factors_3d(
    records: list[WellRecord], proteins: list[str], pH_values: list[float],
    properties: dict[str, MetalProperties], out_path: Path,
    *, metal_axis: str = "radius", vmax: float | None = None,
) -> None:
    """Exact coordinates with activity-scaled area and protein color; overlaps remain.

    Area above a 0.5-point-squared baseline scales with nonnegative kcat.
    """
    plt.rcParams.update(NATURE_RC)
    _assign_protein_colors(proteins)
    points = _condition_means(records)
    if metal_axis == "radius":
        descriptors = {metal: prop.ionic_radius_pm for metal, prop in properties.items()}
        xlabel = "Metal ionic radius (pm)"
    elif metal_axis == "mass":
        descriptors = METAL_MASS_DA
        xlabel = "Metal atomic mass (Da)"
    else:
        raise ValueError(f"Unknown metal axis: {metal_axis}")
    activities = np.array([record.kcat for record in points])
    upper = vmax if vmax is not None else max(float(activities.max()), 0.001)

    def activity_fraction(values):
        return np.clip(np.asarray(values) / upper, 0, 1)

    def marker_area(values):
        return 0.5 + 55.0 * activity_fraction(values)

    fig = plt.figure(figsize=(7.08, 5.20))
    axis = fig.add_subplot(111, projection="3d", computed_zorder=False)
    handles = []
    for protein_index, protein in enumerate(proteins):
        subset = [record for record in points if record.protein == protein]
        protein_activities = np.array([record.kcat for record in subset])
        axis.scatter(
            [descriptors[record.metal] for record in subset],
            [record.mw_da for record in subset],
            [record.pH for record in subset],
            color=PROTEIN_COLORS[protein],
            marker="o", s=marker_area(protein_activities), alpha=0.55, depthshade=False,
            edgecolors="none", linewidths=0, zorder=protein_index + 2,
        )
        handles.append(Line2D([], [], marker="o", markersize=4, markeredgewidth=0,
                              color=PROTEIN_COLORS[protein], alpha=0.55, linestyle="none", label=protein))
    axis.set_xlabel(xlabel, labelpad=7)
    axis.set_ylabel("Substrate MW (Da)", labelpad=7)
    axis.set_zlabel("pH", labelpad=4)
    axis.set_zticks(pH_values, labels=[f"{ph:g}" for ph in pH_values])
    axis.set_zticks([], minor=True)
    axis.set_box_aspect((1.25, 1.25, 0.8))
    axis.view_init(elev=23, azim=-55)
    axis.tick_params(labelsize=7, pad=1)
    for dimension in (axis.xaxis, axis.yaxis, axis.zaxis):
        dimension.pane.fill = False
    axis.legend(handles=handles, loc="upper left", frameon=False)
    levels = np.linspace(0.0, upper, 10)
    size_handles = [Line2D([], [], marker="o", markersize=float(np.sqrt(marker_area(level))),
                           markeredgewidth=0, color="0.4", alpha=0.55,
                           linestyle="none", label=f"{level:.3g}")
                    for level in levels]
    fig.legend(handles=size_handles, title="Marker size\nk$_{cat}$ (s$^{-1}$)",
               loc="upper right", bbox_to_anchor=(0.98, 0.90),
               labelspacing=1.0, markerscale=1.0, frameon=False)
    fig.subplots_adjust(left=0.03, right=0.82, bottom=0.10, top=0.97)
    _finalize(fig, out_path)


def _finalize(fig, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, format="pdf")
    fig.savefig(out_path.with_suffix(".png"), format="png", dpi=600)
    plt.close(fig)


def plot_ph_effect(
    records: list[WellRecord],
    proteins: list[str],
    pH_values: list[float],
    out_path: Path,
    *,
    vmax: float | None = None,
) -> None:
    """Mean±SEM per-well k_cat vs pH, one line per protein (all metals + alcohols)."""
    plt.rcParams.update(NATURE_RC)
    fig, ax = plt.subplots(figsize=(3.46, NATURE_PANEL_INCHES[1]))
    for protein in proteins:
        color = PROTEIN_COLORS[protein]
        xs, ys, yerr = [], [], []
        for pH in pH_values:
            subset = _filter(records, protein=protein, pH=pH)
            mean, std, n = _group_stats([r.kcat for r in subset])
            if n == 0:
                continue
            xs.append(pH)
            ys.append(mean)
            yerr.append(_sem(std, n))
        ax.errorbar(xs, ys, yerr=yerr, color=color, marker="o", ms=3, lw=0.7,
                    capsize=1.5, elinewidth=0.4, label=protein)
    ax.set_xticks(pH_values, labels=[f"{ph:g}" for ph in pH_values])
    ax.set_xticks([], minor=True)
    ax.set_xlabel("pH")
    ax.set_ylabel("k$_{cat}$ (s$^{-1}$)")
    ax.set_title("pH effect (averaged over metals + alcohols)")
    ax.set_ylim(0, vmax)
    ax.legend(frameon=False, loc="best")
    fig.tight_layout()
    _finalize(fig, out_path)


def _grouped_bar(
    ax,
    categories: list[str],
    proteins: list[str],
    pH_values: list[float],
    means: dict[tuple[str, float, str], tuple[float, float, int]],
) -> None:
    """Draw a grouped bar plot: category on x, one bar cluster per (protein, pH)."""
    n_series = len(proteins) * len(pH_values)
    if n_series == 0:
        return
    total_width = 0.8
    bar_width = total_width / n_series
    x = np.arange(len(categories))
    idx = 0
    for protein in proteins:
        base_color = PROTEIN_COLORS[protein]
        for pH in pH_values:
            offsets = x - total_width / 2 + bar_width / 2 + idx * bar_width
            ys = np.array([means.get((protein, pH, cat), (np.nan, np.nan, 0))[0] for cat in categories])
            ns = np.array([means.get((protein, pH, cat), (np.nan, np.nan, 0))[2] for cat in categories])
            stds = np.array([means.get((protein, pH, cat), (np.nan, np.nan, 0))[1] for cat in categories])
            yerr = np.array([_sem(s, int(nn)) for s, nn in zip(stds, ns)])
            alpha = _ph_alpha(pH, pH_values)
            ax.bar(offsets, ys, width=bar_width * 0.95, color=base_color, alpha=alpha,
                   edgecolor="none", linewidth=0)
            ax.errorbar(offsets, ys, yerr=yerr, fmt="none",
                        ecolor="0.2", elinewidth=0.35, capsize=1.0)
            idx += 1
    ax.set_xticks(x)
    ax.set_xticklabels(categories)


def _legend_protein_ph(proteins: list[str], pH_values: list[float]) -> list[Patch]:
    handles: list[Patch] = []
    for protein in proteins:
        base = PROTEIN_COLORS[protein]
        for pH in pH_values:
            alpha = _ph_alpha(pH, pH_values)
            handles.append(Patch(facecolor=base, alpha=alpha,
                                 label=f"{protein}, pH {pH:g}"))
    return handles


def plot_metal_effect(
    records: list[WellRecord],
    proteins: list[str],
    pH_values: list[float],
    out_path: Path,
    *,
    vmax: float | None = None,
) -> None:
    """Grouped bar plot: k_cat vs metal, one cluster per (protein, pH); averaged over alcohols."""
    plt.rcParams.update(NATURE_RC)
    metals_present = [m for m in METALS if any(r.metal == m for r in records)]
    means: dict[tuple[str, float, str], tuple[float, float, int]] = {}
    for protein in proteins:
        for pH in pH_values:
            for metal in metals_present:
                subset = _filter(records, protein=protein, pH=pH, metal=metal)
                means[(protein, pH, metal)] = _group_stats([r.kcat for r in subset])

    fig, ax = plt.subplots(figsize=NATURE_PANEL_INCHES)
    _grouped_bar(ax, metals_present, proteins, pH_values, means)
    ax.set_xlabel("metal")
    ax.set_ylabel("k$_{cat}$ (s$^{-1}$)")
    ax.set_title("Metal effect (averaged over alcohols)")
    ax.set_ylim(0, vmax)
    ax.tick_params(axis="x", rotation=90)
    ax.legend(handles=_legend_protein_ph(proteins, pH_values),
              frameon=False, loc="upper right", ncol=len(proteins))
    fig.tight_layout()
    _finalize(fig, out_path)


def plot_alcohol_effect(
    records: list[WellRecord],
    proteins: list[str],
    pH_values: list[float],
    out_path: Path,
    *,
    vmax: float | None = None,
) -> None:
    """Grouped bar plot: k_cat vs alcohol, one cluster per (protein, pH); averaged over metals."""
    plt.rcParams.update(NATURE_RC)
    alcohols_present = [a for a in ALCOHOLS if any(r.substrate == a for r in records)]
    means: dict[tuple[str, float, str], tuple[float, float, int]] = {}
    for protein in proteins:
        for pH in pH_values:
            for alcohol in alcohols_present:
                subset = _filter(records, protein=protein, pH=pH, substrate=alcohol)
                means[(protein, pH, alcohol)] = _group_stats([r.kcat for r in subset])

    fig, ax = plt.subplots(figsize=NATURE_PANEL_INCHES)
    _grouped_bar(ax, alcohols_present, proteins, pH_values, means)
    ax.set_xlabel("alcohol")
    ax.set_ylabel("k$_{cat}$ (s$^{-1}$)")
    ax.set_title("Alcohol effect (averaged over metals)")
    ax.set_ylim(0, vmax)
    ax.tick_params(axis="x", rotation=45)
    for tick in ax.get_xticklabels():
        tick.set_horizontalalignment("right")
    ax.legend(handles=_legend_protein_ph(proteins, pH_values),
              frameon=False, loc="upper right", ncol=len(proteins))
    fig.tight_layout()
    _finalize(fig, out_path)


def plot_alcohol_size_effect(
    records: list[WellRecord],
    proteins: list[str],
    pH_values: list[float],
    out_path: Path,
    *,
    include_aromatic: bool = False,
    vmax: float | None = None,
) -> None:
    """Per-alcohol mean k_cat vs molecular weight, with per-(protein, pH) linear fit.

    Aromatic substrates are excluded by default because their electronic structure
    changes what "size" means to the enzyme; pass include_aromatic=True to keep them.
    """
    plt.rcParams.update(NATURE_RC)
    fig, ax = plt.subplots(figsize=NATURE_PANEL_INCHES)
    slopes: list[tuple[str, float, float, float, int]] = []  # (protein, pH, slope, intercept, n)

    for protein in proteins:
        base_color = PROTEIN_COLORS[protein]
        for pH in pH_values:
            xs, ys, yerr = [], [], []
            for alcohol, props in ALCOHOL_PROPS.items():
                if props["aromatic"] and not include_aromatic:
                    continue
                subset = _filter(records, protein=protein, pH=pH, substrate=alcohol)
                mean, std, n = _group_stats([r.kcat for r in subset])
                if n == 0:
                    continue
                xs.append(float(props["mw_da"]))
                ys.append(mean)
                yerr.append(_sem(std, n))
            if not xs:
                continue
            alpha = _ph_alpha(pH, pH_values)
            label = f"{protein}, pH {pH:g}"
            ax.errorbar(xs, ys, yerr=yerr, fmt="o", ms=3, mfc=base_color, mec="none",
                        alpha=alpha, ecolor=base_color, elinewidth=0.4, capsize=1.2,
                        linestyle="none", label=label)
            if len(set(xs)) >= 2:
                slope, intercept = np.polyfit(xs, ys, 1)
                slopes.append((protein, pH, float(slope), float(intercept), len(xs)))
                xf = np.array([min(xs), max(xs)], dtype=float)
                ax.plot(xf, slope * xf + intercept, color=base_color,
                        alpha=alpha, lw=0.7)

    ax.set_xlabel("substrate MW (Da)")
    ax.set_ylabel("k$_{cat}$ (s$^{-1}$)")
    aromatic_note = "including aromatic" if include_aromatic else "aliphatic only"
    ax.set_title(f"Alcohol size effect ({aromatic_note})")
    ax.set_ylim(0, vmax)
    ax.legend(frameon=False, loc="best")

    text_lines = [f"{p}, pH {ph:g}: slope = {m:.3g}, n = {n}"
                  for p, ph, m, _b, n in slopes]
    if text_lines:
        ax.text(0.02, 0.98, "\n".join(text_lines),
                transform=ax.transAxes, ha="left", va="top",
                fontsize=7, color="0.25")
    fig.tight_layout()
    _finalize(fig, out_path)


def plot_protein_vs_protein(
    records: list[WellRecord],
    proteins: list[str],
    pH_values: list[float],
    out_path: Path,
    *,
    vmax: float | None = None,
) -> None:
    """Scatter of protein A vs protein B k_cat at matched (pH, metal, alcohol) triples."""
    assert len(proteins) == 2, "protein-vs-protein scatter requires exactly two proteins"
    plt.rcParams.update(NATURE_RC)
    prot_a, prot_b = proteins
    n_panels = len(pH_values)
    ncols = min(2, n_panels)
    nrows = (n_panels + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(7.08 if ncols > 1 else 3.46,
                                      min(3.2 * nrows, 6.69)),
                             squeeze=False)

    for unused_ax in axes.flat[n_panels:]:
        unused_ax.set_visible(False)

    for ax, pH in zip(axes.flat, pH_values):
        pairs = _matched_pairs(records, prot_a, prot_b, pH)
        if not pairs:
            ax.text(0.5, 0.5, "no matched wells", transform=ax.transAxes,
                    ha="center", va="center", color="0.4")
            ax.set_title(f"pH {pH:g}")
            continue
        xa = np.array([p[0] for p in pairs])
        yb = np.array([p[1] for p in pairs])
        colors = [_TAB10[METALS.index(p[2]) % len(_TAB10)] if p[2] in METALS else "0.5"
                  for p in pairs]
        ax.scatter(xa, yb, c=colors, s=6, edgecolors="none", alpha=0.85)
        lim = max(float(np.nanmax(xa)), float(np.nanmax(yb)),
                  vmax if vmax is not None else 0.0)
        if lim <= 0 or not np.isfinite(lim):
            lim = 1.0
        ax.plot([0, lim], [0, lim], color="0.3", lw=0.4, ls="--")
        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel(f"{prot_a} k$_{{cat}}$ (s$^{{-1}}$)")
        ax.set_ylabel(f"{prot_b} k$_{{cat}}$ (s$^{{-1}}$)")
        ax.set_title(f"pH {pH:g} (n = {len(pairs)})")

        finite = np.isfinite(xa) & np.isfinite(yb)
        if finite.sum() >= 2:
            r = float(np.corrcoef(xa[finite], yb[finite])[0, 1])
            ax.text(0.02, 0.98, f"Pearson r = {r:.2f}",
                    transform=ax.transAxes, ha="left", va="top",
                    fontsize=7, color="0.25")

    fig.suptitle(f"{prot_a} vs {prot_b} (matched metal x alcohol)", fontsize=7)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    _finalize(fig, out_path)


def _matched_pairs(
    records: list[WellRecord],
    prot_a: str,
    prot_b: str,
    pH: float,
) -> list[tuple[float, float, str, str]]:
    """Aggregate to per-(metal, alcohol) means for each protein, then match."""
    def agg(protein: str) -> dict[tuple[str, str], float]:
        out: dict[tuple[str, str], float] = {}
        for metal in {r.metal for r in records if r.protein == protein and r.pH == pH}:
            for alcohol in {r.substrate for r in records
                            if r.protein == protein and r.pH == pH and r.metal == metal}:
                subset = _filter(records, protein=protein, pH=pH,
                                 metal=metal, substrate=alcohol)
                mean, _std, n = _group_stats([r.kcat for r in subset])
                if n:
                    out[(metal, alcohol)] = mean
        return out

    a = agg(prot_a)
    b = agg(prot_b)
    shared = sorted(set(a) & set(b))
    return [(a[key], b[key], key[0], key[1]) for key in shared]


# =============================================================================
if __name__ == "__main__":
    main()
