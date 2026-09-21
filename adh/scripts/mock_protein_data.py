"""Generate explicitly synthetic averaged PedE/PedH screening CSVs.

Run with python3 adh/scripts/mock_protein_data.py, then pass the generated
mock_manifest.json to compare_proteins.py with --label mock.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from kcat_plate import ALCOHOLS, PLATE_COLS, PLATE_ROWS, make_screen_plate_layout


PH_RATES = {"PedE": {6: 0.0, 7: 0.012, 8: 0.085},
            "PedH": {6: 0.45, 7: 0.32, 8: 0.22}}
LIGHT_LANTHANIDES = {"La": 0.90, "Ce": 1.0, "Pr": 0.95, "Nd": 0.88, "Sm": 0.80}
SUBSTRATE_FACTORS = {
    "methanol": 0.75,
    "ethanol": 1.0,
    "1-propanol": 0.95,
    "1-butanol": 0.90,
    "3-methyl-1-butanol": 0.80,
    "1-pentanol": 0.85,
    "2-methyl-1-butanol": 0.78,
    "1-hexanol": 0.75,
    "1,5-pentanediol": 0.82,
    "2-phenylethanol": 0.70,
    "vanillyl alcohol": 0.68,
}


def generate_mock_data(outdir: Path, seed: int = 20260921) -> Path:
    """Write six averaged CSVs with three simulated replicates per assay well."""
    outdir.mkdir(parents=True, exist_ok=True)
    random = np.random.default_rng(seed)
    layout = make_screen_plate_layout()
    shared_substrates = set(ALCOHOLS[: ALCOHOLS.index("2-phenylethanol") + 1])
    conditions = []
    header = ["well", "row", "col", "metal", "substrate", "kcat_mean_s",
              "kcat_std_s", "n_replicates", "kcat_mock_rep1_s", "kcat_mock_rep2_s",
              "kcat_mock_rep3_s", "data_origin", "mock_expected_active"]
    for protein, rates in PH_RATES.items():
        for ph, base_rate in rates.items():
            filename = f"mock_{protein}_pH{ph}_avg_wells.csv"
            with (outdir / filename).open("w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(header)
                for row_name in PLATE_ROWS:
                    for column in PLATE_COLS:
                        well = f"{row_name}{column}"
                        metal, substrate = layout.condition(well)
                        metadata = [well, row_name, column, metal, substrate]
                        if layout.is_water_control(well) or layout.is_tcep_control(well):
                            writer.writerow(metadata + ["", "", 0, "", "", "", "synthetic", ""])
                            continue
                        metal_factor = (float(metal == "Ca") if protein == "PedE"
                                        else LIGHT_LANTHANIDES.get(metal, 0.0))
                        substrate_active = (substrate in shared_substrates or
                                            (protein == "PedH" and substrate == "vanillyl alcohol"))
                        expected = base_rate * metal_factor * SUBSTRATE_FACTORS.get(substrate, 0.0)
                        active = substrate_active and expected > 0
                        if active:
                            replicates = expected * np.clip(random.normal(1.0, 0.05, 3), 0.85, 1.15)
                        else:
                            replicates = random.uniform(0.0, 0.0005, 3)
                        writer.writerow(metadata + [float(np.mean(replicates)),
                                        float(np.std(replicates, ddof=0)), 3,
                                        *replicates.tolist(), "synthetic", int(active)])
            conditions.append({"protein": protein, "pH": ph, "csv": filename})
    manifest = {
        "data_origin": "SYNTHETIC MOCK DATA - NOT EXPERIMENTAL RESULTS",
        "seed": seed,
        "units": "kcat in s^-1",
        "replicates": 3,
        "std_ddof": 0,
        "assumptions": [
            "Absolute rates and substrate preferences are illustrative, not fitted to real data.",
            "PedE is active only with Ca, at pH 8 and weakly at pH 7; inactive at pH 6.",
            "Light lanthanides are La, Ce, Pr, Nd, Sm; Pm is not on the source plate.",
            "PedH is active only with those light lanthanides, with pH 6 > 7 > 8.",
            "Active PedH combinations exceed the strongest active PedE combination.",
            "Both proteins accept the source substrate list through 2-phenylethanol inclusive.",
            "PedH additionally accepts vanillyl alcohol; all other substrates are inactive.",
            "Inactive assay wells have background between 0 and 0.0005 s^-1.",
            "Control rows retain the source layout but have blank kcat and zero replicates.",
            "Population standard deviations match kcat_plate.dump_averaged_wells.",
        ],
        "conditions": conditions,
    }
    manifest_path = outdir / "mock_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--outdir", type=Path,
                        default=Path(__file__).resolve().parents[1] / "outputs" / "mock_proteins")
    parser.add_argument("--seed", type=int, default=20260921)
    args = parser.parse_args()
    print(f"Synthetic data manifest: {generate_mock_data(args.outdir, args.seed)}")


if __name__ == "__main__":
    main()