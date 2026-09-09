# TOPAS TPS Source Extension

Pencil beam scanning (PBS) particle source for [TOPAS](https://github.com/OpenTOPAS/OpenTOPAS) / OpenTOPAS.

This is a TOPAS source extension (`TsSource` + `TsVGenerator` pair). It reads a compact spot-plan CSV plus an energy-dependent machine beam model, then generates primaries with scanning-magnet steering and BiGaussian emittance. Treatment-plan vendor formats and DICOM stay outside this source — convert them to CSV first.

Current version is MT-safe: each Geant4 event ID in `[0, N)` maps to exactly one history on one spot via a read-only prefix-sum table, so workers generate disjoint spots without shared mutable state.

## Layout

```text
include/
  TsSourcePencilBeamScanning.hh      # source, owns spot preparation + history->spot map
  TsGeneratorPencilBeamScanning.hh   # generator, samples one primary per event
  TsPBSSpotPlan.hh                   # spot-plan CSV loader
  TsPBSBeamModel.hh                  # beam-model CSV loader + lookup/interpolation
  TsPBSCoordinateModel.hh            # TPS / ComponentLocal ray geometry
src/
  TsSourcePencilBeamScanning.cc
  TsGeneratorPencilBeamScanning.cc
  TsPBSSpotPlan.cc
  TsPBSBeamModel.cc
  TsPBSCoordinateModel.cc
example/
  pbs_water.txt                # minimal runnable TOPAS parameter file (water phantom)
  spots.csv                    # 6-spot plan, 1200 histories total
  beam_model.csv               # machine optics for the spot energies
```

Drop these files into a TOPAS extension build (source + generator + extra classes) and register the `PencilBeamScanning` source type name.

## Input files

Both CSVs are case-insensitive on headers, tolerate surrounding whitespace, and skip blank lines and `#` comment lines. First non-comment line must be the header.

### Spot plan (`SpotPlanFile`)

Required columns: `x`, `y`, `energy`, `weight`. Optional: `spot_id` / `id`.

Aliases: `x_mm` for `x`, `y_mm` for `y`, `energy_mev` for `energy`.

| column | unit | meaning |
| --- | --- | --- |
| `x`, `y` | mm, isocenter plane | spot position |
| `energy` | MeV (per nucleon-equivalent as passed to TOPAS) | spot energy |
| `weight` | histories (see `WeightMode`) | primaries for this spot |
| `spot_id` / `id` | non-negative integer, optional | kept as `TsPBSSpot::id`; defaults to row index |

Example:

```csv
x,y,energy,weight
-40.0,12.0,150.2,1000
-30.0,12.0,150.2,1000
```

Validation: `energy > 0`, `weight >= 0`. Zero-weight rows are counted and skipped by default (`SkipZeroWeightSpots`). A plan with no usable spots aborts the session.

### Beam model (`BeamModelFile`)

Required columns: `energy`, `sigma_x_mm`, `sigma_xp_rad`, `corr_x`, `sigma_y_mm`, `sigma_yp_rad`, `corr_y`, `energy_spread_percent`.

Aliases: `energy_mev`, `sigmax`, `sigmaxprime`, `correlationx`, `sigmay`, `sigmayprime`, `correlationy`, `energyspread`.

Example:

```csv
energy,sigma_x_mm,sigma_xp_rad,corr_x,sigma_y_mm,sigma_yp_rad,corr_y,energy_spread_percent
150.2,3.1,0.0042,-0.55,3.3,0.0045,-0.60,0.35
200.5,2.6,0.0036,-0.50,2.8,0.0038,-0.55,0.30
```

Validation: `energy > 0`, sigmas `>= 0`, correlations in `[-1, 1]`, energy spread `>= 0`, no duplicate energies (relative tolerance `1e-6`). Rows are sorted by energy. Exact-match lookup always works; in-between energies require `InterpolateBeamModel = "True"` (linear interpolation), otherwise the run aborts. Out-of-range energies always abort.

## Parameters

```text
s:So/CarbonPBS/Type = "PencilBeamScanning"
s:So/CarbonPBS/Component = "PBSBeamFrame"
s:So/CarbonPBS/BeamParticle = "GenericIon(6,12,6)"

s:So/CarbonPBS/SpotPlanFile = "spots.csv"
s:So/CarbonPBS/BeamModelFile = "beam_model.csv"

d:So/CarbonPBS/VirtualScanningMagneticX = 6227.8 mm
d:So/CarbonPBS/VirtualScanningMagneticY = 7008.6 mm
d:So/CarbonPBS/VirtualSourceToIsocenterDistance = 450 mm

b:So/CarbonPBS/SkipZeroWeightSpots = "True"
s:So/CarbonPBS/WeightMode = "Histories"
u:So/CarbonPBS/HistoriesScale = 1.0
b:So/CarbonPBS/InterpolateBeamModel = "False"
i:So/CarbonPBS/FirstSpot = 0
i:So/CarbonPBS/LastSpot = -1
```

| parameter | type | default | meaning |
| --- | --- | --- | --- |
| `SpotPlanFile` | string | — (required) | spot-plan CSV path |
| `BeamModelFile` | string | — (required) | beam-model CSV path |
| `VirtualScanningMagneticX` (`VirtualSADX`) | Length | — (required, `> 0`) | virtual SAD for X steering |
| `VirtualScanningMagneticY` (`VirtualSADY`) | Length | — (required, `> 0`) | virtual SAD for Y steering |
| `VirtualSourceToIsocenterDistance` (`SAD`, `SourceToIsocenterDistance`) | Length | — (required, `> 0`) | source-plane to isocenter distance |
| `WeightMode` | string | `"Histories"` | only `Histories` is supported; CSV `weight` = number of primaries |
| `HistoriesScale` | unitless | `1.0` | per-spot histories = `llround(weight * HistoriesScale)`; must be `>= 0`; spots rounding to `<= 0` are skipped |
| `SkipZeroWeightSpots` | bool | `"True"` | skip `weight == 0` rows (still counted in the log) |
| `InterpolateBeamModel` | bool | `"False"` | allow linear interpolation between beam-model energies |
| `FirstSpot` / `LastSpot` | int | `0` / `-1` | inclusive sub-range of the loaded plan; `-1` means last spot |
| `SpotCoordinateConvention` (`CoordinateConvention`) | string | `"TPS"` | `TPS` (also accepts `IEC61217`) or `ComponentLocal` |

Notes:

- `NumberOfHistoriesInRun` is set automatically from the selected spots (`fNumberOfHistoriesInRun = total`). Do not set it manually and do not use Time Feature vectors for per-spot energy or optics.
- Total histories are capped at `1e9`; empty selections abort with an error.
- Only `WeightMode = Histories` exists in this version. There is no MU-to-particle calibration.

## Geometry

Hang `PBSBeamFrame` on the IEC gantry. The source only needs `s:So/.../Component = "PBSBeamFrame"`. Gantry / couch live on `IEC_G` / `IEC_S`.

```text
s:Ge/IEC_F/Parent = "World"
s:Ge/IEC_F/Type   = "Group"
s:Ge/IEC_G/Parent = "IEC_F"
s:Ge/IEC_G/Type   = "Group"

s:Ge/PBSBeamFrame/Parent = "IEC_G"
s:Ge/PBSBeamFrame/Type   = "Group"
```

- `PBSBeamFrame` sits at the isocenter. Do not translate or rotate it to place the nozzle — change `IEC_G` rotations for other gantry angles.
- The source plane is at `y = -VirtualSourceToIsocenterDistance` (TPS convention) and every spot is aimed at the isocenter origin.
- `VirtualScanningMagneticX/Y` are the scanning-magnet virtual distances used for spot steering: `thetaX = atan(xIso / VSADX)`, `thetaY = atan(yIso / VSADY)`.

### Gantry angle: TPS angle = −IEC_G

A TPS gantry angle `θ` is set as `Ge/IEC_G/RotZ = −θ`, because OpenTOPAS uses the IEC rotation convention `world = Rz(−RotZ) · local`:

- The local beam points along component +Y. Rotating `IEC_G` by `RotZ` moves the world-space beam to `(sin RotZ, cos RotZ, 0)`.
- The TPS convention for gantry `θ` is `w(θ) = (−sin θ, cos θ, 0)` (θ = 0 → +Y, θ = 90° → −X, θ = 270° → +X).
- Matching the two gives `RotZ = −θ`.

| TPS gantry `θ` | `d:Ge/IEC_G/RotZ` | world-space beam |
| --- | --- | --- |
| 0° | `0. deg` | +Y (`y-` → `y+`) |
| 90° | `270. deg` (= −90°) | −X |
| 180° | `180. deg` | −Y |
| 270° | `90. deg` | +X |

Example for TPS 90°:

```text
d:Ge/IEC_G/RotZ = 270. deg
```

WARNING: do not set `IEC_G/RotZ = 90` to mean "TPS 90°" — that points the beam at +X, which is TPS 270°. (The validated full-plan cases keep `IEC_G/RotZ = 0` and instead rotate the CT with `Ge/Patient/RotZ = θ`; that "rotate CT, not the beam" setup is equivalent for dose comparison but is a different way of expressing the same irradiation — pick one and keep it consistent.)

### Coordinate conventions

Default `TPS` (built-in IEC mapping, no user rotation needed):

- component origin = isocenter
- component +X = TPS scan X
- component +Y = beam (gantry 0: `y-` → `y+`)
- component +Z = TPS scan Y
- direction: `(tx, 1, ty) / sqrt(tx² + 1 + ty²)` with `tx = tan(thetaX)`, `ty = tan(thetaY)`

`ComponentLocal` leaves +Z as the beam (`(tx, ty, 1)` normalized) so a user-supplied component rotation can remap the axes. Only set `SpotCoordinateConvention` when overriding the default.

## Sampling

Per event, the generator looks up the owning spot from the event ID (`SpotForHistory`), then:

- Transverse phase space: correlated BiGaussian. With `ux, vx, uy, vy ~ N(0,1)`:
  `dx = sigmaX·ux`, `x' = sigmaXp·(corrX·ux + vx·sqrt(1-corrX²))`, same for Y.
- Direction: nominal steering angles plus sampled `x'`/`y'`, converted with `tan` and normalized (see convention above).
- Energy: `E ~ N(E_spot, E_spot·spread%/100)`; resampled while `E <= 0`. Zero spread gives monoenergetic spots.
- Start position: nominal source-plane point plus sampled offsets (`dx` along component X; `dy` along component Z in TPS mode, Y in `ComponentLocal` mode).
- `p.weight = 1`, `isNewHistory = true`, particle type from `BeamParticle`, then `TransformPrimaryForComponent`.

## Multithreading

`Ts/NumberOfThreads` may be greater than 1. The history-to-spot table (`fHistoryBegin` prefix sums) is built once in `ResolveParameters` and read-only afterwards; `SpotForHistory` is a binary search (`upper_bound`) with no mutable state. Geant4 thread-local RNGs sample optics independently per event. Time-ordered delivery / interplay is still not modeled — MT workers generate spots concurrently.

## Example

`example/` is a minimal runnable case (water phantom, no CT data needed):

- `pbs_water.txt` — main TOPAS parameter file. Run from inside `example/` with a TOPAS build that includes this extension: `topas pbs_water.txt`.
- `spots.csv` — 6 spots / 1200 histories (spot 4 has weight 0 to exercise `SkipZeroWeightSpots`).
- `beam_model.csv` — machine optics rows covering every spot energy, so the default `InterpolateBeamModel = "False"` works.

Provenance and adaptations: machine numbers (`VirtualScanningMagneticX/Y = 6227.8 / 7008.6 mm`, `SAD = 450 mm`, beam optics, spot positions/energies) come from the RT06423 full-plan case (`run_full_plan.txt`, `spots.csv`, `beam_model.csv`), which instead simulates ~1100 spots / ~15M histories on a DICOM patient with custom LET scorers. The example replaces the DICOM patient with a `G4_WATER` box, drops the custom `myHadronLET` scorers (they need a separate scorer extension), shrinks the plan to 6 spots with small weights, spells out the current-version parameters (`InterpolateBeamModel`, `FirstSpot`, `LastSpot`), and runs 4 threads to demonstrate the MT-safe mapping. Gantry is TPS 0° (`IEC_G/RotZ = 0`); see "Gantry angle" above for other angles.

## Limitations

- No DICOM / vendor plan parser (convert to CSV first)
- No MU-to-particle calibration (`HistoriesScale` is a plain multiplier)
- No delivery timing / interplay; no time-ordered spot delivery (MT is concurrent)
- BiGaussian emittance only
- `WeightMode = Histories` only

## License

MIT — see [LICENSE](LICENSE).
