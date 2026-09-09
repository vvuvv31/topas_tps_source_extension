#!/usr/bin/env python3
"""Assemble per-spot TOPAS dose binaries into a sparse Dij matrix.

Reads dij_index.csv (written by run_dij.py) plus one DoseToMedium .bin /
.binheader pair per spot, thresholds each 3D dose, and stacks the spots as
columns of a sparse matrix:

    Dij : (nVoxels, nSpots) sparse CSC, float64, dose in Gy
          for DijHistoriesPerSpot histories per column.

Voxel ordering follows the TOPAS binary layout (X fastest):
    flat row r = (iz * ny + iy) * nx + ix   (0-based),
i.e. np.fromfile(path, dtype='<f8').reshape((nz, ny, nx)) in C order.

Outputs (next to --out, which should end in .mat or .npz):
  <out>.npz       sparse CSC matrix (scipy) + <out>_meta.json sidecar
  <out>.mat       MATLAB-compatible v5 file: Dij loads as a MATLAB sparse
                  matrix, plus spot table and grid metadata
  (both are written unless --no-npz / --no-mat is given)

Dose scaling for optimization: with optimizer particle counts x (one entry
per Dij column), the physical dose is Dij * (x / H), where H is the
histories-per-spot stored as histories_per_spot in the outputs.

Example:
  python3 assemble_dij.py --outdir dij_out --out dij_out/dij \\
      --threshold 1e-9 --relative-threshold 0.01
"""

import argparse
import csv
import json
import re
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.io as sio

HEADER_RE = re.compile(r"#\s*([XYZ])\s+in\s+(\d+)\s+bins\s+of\s+([\d.eE+-]+)\s*(\w*)")
UNIT_TO_MM = {"mm": 1.0, "cm": 10.0, "m": 1000.0, "": 10.0}


def parse_binheader(path):
    """Return (nx, ny, nz), (dx, dy, dz) in mm, plus raw meta lines."""
    dims = {}
    spacing = {}
    scorer = ""
    component = ""
    for line in Path(path).read_text().splitlines():
        match = HEADER_RE.match(line.strip())
        if match:
            axis, nbins, size, unit = match.groups()
            dims[axis] = int(nbins)
            spacing[axis] = float(size) * UNIT_TO_MM.get(unit.lower(), 10.0)
            continue
        if "Results for scorer" in line:
            scorer = line.split(":", 1)[-1].strip()
        if "Scored in component" in line:
            component = line.split(":", 1)[-1].strip()
    if set(dims) != {"X", "Y", "Z"}:
        raise ValueError(f"could not parse grid dims from {path}")
    return (dims["X"], dims["Y"], dims["Z"]), \
        (spacing["X"], spacing["Y"], spacing["Z"]), scorer, component


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", required=True, help="run_dij.py output dir")
    parser.add_argument("--out", required=True, help="output stem (dij -> dij.npz/dij.mat)")
    parser.add_argument("--threshold", type=float, default=0.0,
                        help="absolute cut in Gy: keep voxels with dose > threshold "
                             "(default 0 = all nonzero)")
    parser.add_argument("--relative-threshold", type=float, default=0.0,
                        help="per-column relative cut: keep voxels with dose > "
                             "relative-threshold * column max (0 = off). "
                             "Effective cut = max(threshold, "
                             "relative-threshold * column max).")
    parser.add_argument("--no-npz", action="store_true")
    parser.add_argument("--no-mat", action="store_true")
    args = parser.parse_args()

    if args.threshold < 0.0:
        parser.error("--threshold must be >= 0")
    if args.relative_threshold < 0.0:
        parser.error("--relative-threshold must be >= 0")

    outdir = Path(args.outdir)
    index_path = outdir / "dij_index.csv"
    with open(index_path, newline="") as handle:
        index_rows = list(csv.DictReader(handle))
    if not index_rows:
        raise ValueError(f"no columns in {index_path}")

    histories = {int(r["histories"]) for r in index_rows}
    if len(histories) != 1:
        raise ValueError(f"mixed histories per spot in {index_path}: {sorted(histories)}")
    histories_per_spot = histories.pop()

    grid = None
    data_rows, data_cols, data_vals = [], [], []
    col_max = []
    for row in index_rows:
        column = int(row["column"])
        bin_path = outdir / row["bin_file"]
        head_path = bin_path.with_suffix(".binheader")
        dims, spacing, scorer, component = parse_binheader(head_path)
        if grid is None:
            grid = (dims, spacing, scorer, component)
        elif dims != grid[0]:
            raise ValueError(f"grid change at {bin_path}: {dims} vs {grid[0]}")
        nx, ny, nz = dims
        expected = nx * ny * nz
        values = np.fromfile(bin_path, dtype="<f8")
        if values.size != expected:
            raise ValueError(f"{bin_path}: {values.size} values, expected {expected}")
        mask = values > args.threshold
        colmax = float(values.max(initial=0.0))
        cut = args.threshold
        if args.relative_threshold > 0.0:
            cut = max(cut, args.relative_threshold * colmax)
            mask = values > cut
        hit = np.flatnonzero(mask)
        data_rows.append(hit)
        data_cols.append(np.full(hit.shape, column, dtype=np.int64))
        data_vals.append(values[hit])
        col_max.append(colmax)
        print(f"[col {column}] {bin_path.name}: nnz={hit.size}/{expected} "
              f"max={col_max[-1]:.6g} Gy cut={cut:.6g} Gy", flush=True)

    n_vox = int(np.prod(grid[0]))
    n_spots = len(index_rows)
    Dij = sp.csc_matrix(
        (np.concatenate(data_vals),
         (np.concatenate(data_rows), np.concatenate(data_cols))),
        shape=(n_vox, n_spots), dtype=np.float64)
    Dij.sort_indices()
    print(f"Dij shape={Dij.shape} nnz={Dij.nnz} "
          f"density={Dij.nnz / Dij.shape[0] / Dij.shape[1]:.4g}")

    (nx, ny, nz), (dx, dy, dz), scorer, component = grid
    meta = {
        "nx": nx, "ny": ny, "nz": nz,
        "voxel_size_x_mm": dx, "voxel_size_y_mm": dy, "voxel_size_z_mm": dz,
        "voxel_order": "x_fastest: r=(iz*ny+iy)*nx+ix, 0-based",
        "dose_unit": "Gy",
        "histories_per_spot": histories_per_spot,
        "scorer": scorer, "component": component,
        "threshold_gy": args.threshold,
        "relative_threshold": args.relative_threshold,
        "threshold_note": "per-column cut = max(threshold_gy, "
                          "relative_threshold * column max); voxel kept if dose > cut",
        "dose_scale_note": "physical dose = Dij * (x / histories_per_spot), "
                           "x = optimizer particle counts per column",
    }
    spot_ids = [r["spot_id"] for r in index_rows]
    meta["columns"] = [
        {"column": int(r["column"]), "plan_index": int(r["plan_index"]),
         "spot_id": r["spot_id"], "x_mm": float(r["x_mm"]),
         "y_mm": float(r["y_mm"]), "energy_mev": float(r["energy_mev"])}
        for r in index_rows]

    out = Path(args.out)
    if out.suffix in (".mat", ".npz"):
        out = out.with_suffix("")
    if not args.no_npz:
        sp.save_npz(str(out.with_suffix(".npz")), Dij)
        out.with_suffix(".json").write_text(json.dumps(meta, indent=1))
        print(f"wrote {out.with_suffix('.npz')} + {out.with_suffix('.json')}")
    if not args.no_mat:
        # scipy sparse -> MATLAB sparse on load; plain v5 file keeps
        # compatibility with load() on all MATLAB versions.
        sio.savemat(str(out.with_suffix(".mat")), {
            "Dij": Dij,
            "spot_id": np.array(spot_ids, dtype=object),
            "spot_x_mm": np.array([c["x_mm"] for c in meta["columns"]]),
            "spot_y_mm": np.array([c["y_mm"] for c in meta["columns"]]),
            "spot_energy_mev": np.array([c["energy_mev"] for c in meta["columns"]]),
            "plan_index": np.array([c["plan_index"] for c in meta["columns"]]),
            "histories_per_spot": float(histories_per_spot),
            "nx": nx, "ny": ny, "nz": nz,
            "voxel_size_x_mm": dx, "voxel_size_y_mm": dy, "voxel_size_z_mm": dz,
            "dose_unit": "Gy",
            "voxel_order": meta["voxel_order"],
        }, format="5")
        print(f"wrote {out.with_suffix('.mat')} (MATLAB: load() -> sparse Dij)")


if __name__ == "__main__":
    sys.exit(main())
