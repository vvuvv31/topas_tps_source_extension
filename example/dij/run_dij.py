#!/usr/bin/env python3
"""Run one TOPAS simulation per spot (Dij beamlet) and index the outputs.

Each run simulates a single plan spot with DijMode (fixed histories per
spot) and writes one 3D DoseToMedium binary. Afterwards assemble_dij.py
combines the per-spot dose files into a sparse Dij matrix (.npz + MATLAB
compatible .mat).

The template parameter file must contain these tokens (rendered per spot):
  @FIRST_SPOT@, @LAST_SPOT@, @DIJ_HISTORIES@, @OUTPUT_FILE@, @SEED@

Plan-index convention mirrors TsPBSSpotPlan::Load: the first non-empty,
non-# line is the header; every following data line is one beamlet with
plan_index 0, 1, 2, ... Dij column k always corresponds to the plan spot
with index first + k.

Example (from example/dij):
  python3 run_dij.py --template pbs_dij_template.txt --spots ../spots.csv \\
      --outdir dij_out --histories-per-spot 100000
"""

import argparse
import csv
import subprocess
import sys
from pathlib import Path

REQUIRED_TOKENS = ("@FIRST_SPOT@", "@LAST_SPOT@", "@DIJ_HISTORIES@",
                   "@OUTPUT_FILE@", "@SEED@")


def read_plan_rows(spots_path):
    """Return data rows as dicts with plan_index, mirroring the C++ loader."""
    with open(spots_path, newline="") as handle:
        lines = [line for line in handle
                 if line.strip() and not line.lstrip().startswith("#")]
    if not lines:
        raise ValueError(f"spot plan has no header: {spots_path}")
    header = [h.strip().lower() for h in lines[0].split(",")]

    def col(*names):
        for name in names:
            if name in header:
                return header.index(name)
        return -1

    xcol = col("x", "x_mm")
    ycol = col("y", "y_mm")
    ecol = col("energy", "energy_mev")
    idcol = col("spot_id", "id")
    if min(xcol, ycol, ecol) < 0:
        raise ValueError(f"spot plan is missing x/y/energy columns: {spots_path}")
    rows = []
    for plan_index, line in enumerate(lines[1:]):
        fields = [f.strip() for f in line.split(",")]
        rows.append({
            "plan_index": plan_index,
            "spot_id": fields[idcol] if idcol >= 0 and idcol < len(fields) and fields[idcol] else str(plan_index),
            "x_mm": fields[xcol],
            "y_mm": fields[ycol],
            "energy_mev": fields[ecol],
        })
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template", required=True, help="parameter template with @TOKENS@")
    parser.add_argument("--spots", required=True, help="spot-plan CSV (for indexing)")
    parser.add_argument("--outdir", required=True, help="output directory")
    parser.add_argument("--topas", default="topas", help="TOPAS executable")
    parser.add_argument("--histories-per-spot", type=int, required=True)
    parser.add_argument("--first", type=int, default=0, help="first plan index (inclusive)")
    parser.add_argument("--last", type=int, default=-1, help="last plan index (inclusive, -1 = end)")
    parser.add_argument("--seed", type=int, default=1, help="base seed; run seed = base + column")
    parser.add_argument("--max-spots", type=int, default=-1, help="cap on number of spots (smoke test)")
    parser.add_argument("--dry-run", action="store_true", help="render files without running TOPAS")
    parser.add_argument("--no-skip-existing", action="store_true", help="re-run even if outputs exist")
    args = parser.parse_args()

    if args.histories_per_spot <= 0:
        parser.error("--histories-per-spot must be > 0")

    template = Path(args.template).read_text()
    missing = [t for t in REQUIRED_TOKENS if t not in template]
    if missing:
        parser.error(f"template is missing tokens: {', '.join(missing)}")

    rows = read_plan_rows(args.spots)
    last = len(rows) - 1 if args.last < 0 else args.last
    if not (0 <= args.first <= last < len(rows)):
        parser.error(f"invalid range first={args.first} last={last} for {len(rows)} spots")
    selected = rows[args.first:last + 1]
    if args.max_spots >= 0:
        selected = selected[:args.max_spots]

    outdir = Path(args.outdir)
    rundir = outdir / "runs"
    dosedir = outdir / "dose"
    rundir.mkdir(parents=True, exist_ok=True)
    dosedir.mkdir(parents=True, exist_ok=True)

    index_rows = []
    for column, row in enumerate(selected):
        tag = f"spot_{row['plan_index']:06d}"
        txt_path = rundir / f"{tag}.txt"
        out_name = f"dose/{tag}_dose"
        bin_path = outdir / (out_name + ".bin")
        header_path = outdir / (out_name + ".binheader")
        txt_path.write_text(
            template
            .replace("@FIRST_SPOT@", str(row["plan_index"]))
            .replace("@LAST_SPOT@", str(row["plan_index"]))
            .replace("@DIJ_HISTORIES@", str(args.histories_per_spot))
            .replace("@OUTPUT_FILE@", out_name)
            .replace("@SEED@", str(args.seed + column)))
        index_rows.append({
            "column": column,
            "plan_index": row["plan_index"],
            "spot_id": row["spot_id"],
            "x_mm": row["x_mm"],
            "y_mm": row["y_mm"],
            "energy_mev": row["energy_mev"],
            "histories": args.histories_per_spot,
            "bin_file": out_name + ".bin",
        })
        if args.dry_run:
            continue
        if bin_path.exists() and header_path.exists() and not args.no_skip_existing:
            print(f"[skip] {tag} outputs exist", flush=True)
            continue
        print(f"[run] {tag} plan_index={row['plan_index']} "
              f"seed={args.seed + column} histories={args.histories_per_spot}", flush=True)
        subprocess.run([args.topas, str(txt_path)], cwd=outdir, check=True)
        if not (bin_path.exists() and header_path.exists()):
            raise RuntimeError(f"{tag}: TOPAS finished but {out_name}.bin(.header) missing")

    index_path = outdir / "dij_index.csv"
    with open(index_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(index_rows[0].keys()))
        writer.writeheader()
        writer.writerows(index_rows)
    print(f"wrote {index_path} ({len(index_rows)} columns)")
    if args.dry_run:
        print("dry run: parameter files rendered, nothing executed")


if __name__ == "__main__":
    sys.exit(main())
