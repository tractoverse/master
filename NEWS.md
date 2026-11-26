# master 0.0.0.9000

## Pyramidal-tract model (per-bundle somatotopy, D8)

New module `R/pyramid.R` implements the per-bundle pyramidal-tract model
recorded in `DECISIONS.qmd` (D8):

- `simulate_pyramidal_tract()` — one hemisphere's pyramidal tract as a
  `fiber::bundle_set` of five per-body-part bundles (`leg`, `trunk`, `arm`,
  `hand`, `face`).
- `pyt_defaults()` — the full D8 parameter table (stations, somatotopic
  `q`/`r` offsets, cross-section widths, twist, decussation and perturbation
  controls).
- `pyt_stations()` — the anatomical station grid (cortex → cord) with
  fractions and mm depths.
- `smoothstep_interp()` — exported C1 smoothstep interpolation helper,
  extracted from the centerline machinery for public reuse.

Geometry highlights (validated numerically):

- Parcels form a contiguous somatotopic fan at the cortex, funnel through the
  internal capsule, rotate at the crus–midbrain, and re-sort medial–lateral by
  the cord.
- Decussation is a per-streamline `+1 → −1` C1 sign flip (`cross` controls)
  with an optional posterior detour; crossed fibers land contralateral,
  uncrossed stay ipsilateral, and the two hemispheres interlace at the
  medullary crossing.
- The `face` parcel is corticobulbar: it terminates at the pons (fewer
  arc-length samples, `term_station = "pons"`) with the crossed branch
  terminating bilaterally.
- Per-stream metadata (`parcel`, `branch`, `hemisphere`) is lifted to the
  bundle level and flows into `rtists::plot3d()` tooltips.
- Axon spread is rescaled by the local section half-width, so bundles read as
  narrow contiguous clouds (~1.6 mm RMS at the cortex, ~0.3 mm at the cord)
  instead of hairline strokes; the zero-mean (ground-truth pathway) property is
  preserved.

`create_axons()` gained `...` so caller-supplied per-stream labels are stored
in `streamline_data` (`bundle_id`, plus any extras).
