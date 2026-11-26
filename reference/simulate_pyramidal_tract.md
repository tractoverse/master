# Simulate one pyramidal tract (one hemisphere) as a `fiber::bundle_set`

Generates the pyramidal tract of a single hemisphere as a
[`fiber::bundle_set`](https://tractoverse.github.io/fiber/reference/bundle_set.html):
five per-parcel bundles (`leg`, `trunk`, `arm`, `hand`, `face`) laid out
by the somatotopic table
[`pyt_defaults()`](https://tractoverse.github.io/master/reference/pyt_defaults.md).
Each parcel is a centreline swept along a common macro-axis with a
smooth (C1) moving frame, then perturbed into a bundle of axons with
[`create_axons()`](https://tractoverse.github.io/master/reference/create_axons.md).

## Usage

``` r
simulate_pyramidal_tract(
  parcels = NULL,
  macros = NULL,
  rot_st = NULL,
  cross = NULL,
  perturbation = NULL,
  hemisphere = "left",
  seed = NULL
)
```

## Arguments

- parcels:

  Optional character vector of parcel names to include (defaults to all
  five from
  [`pyt_defaults()`](https://tractoverse.github.io/master/reference/pyt_defaults.md)).

- macros:

  Optional list with `tract_length`, `midline_offset`, `bow`, `n_points`
  (defaults from
  [`pyt_defaults()`](https://tractoverse.github.io/master/reference/pyt_defaults.md)).

- rot_st:

  Optional numeric control points for the twist map `psi` (length equal
  to `nrow(pyt_defaults()$q_st)`).

- cross:

  Optional list with `cross_percent`, `cross_span`, `cross_detour`,
  `center_gap`.

- perturbation:

  Optional list with `variance`, `autocorrelation`, `n_basis`, `n_axons`
  (reused by
  [`create_axons()`](https://tractoverse.github.io/master/reference/create_axons.md)).

- hemisphere:

  Either `"left"` or `"right"`; the medial wall faces the midline (`+X`
  toward the neighbour hemisphere).

- seed:

  Optional integer seed for full reproducibility.

## Value

A
[`fiber::bundle_set`](https://tractoverse.github.io/fiber/reference/bundle_set.html)
of one bundle per parcel. Each stream's `streamline_data` carries
`parcel` and `branch` (`crossed`/`uncrossed`).

## Details

The `face` parcel terminates at the pons (corticobulbar); the limb/trunk
parcels continue through the decussation into the cord. Each parcel
contains a crossed and an uncrossed branch: the crossed branch's `q` is
multiplied by a smooth sign-flip `phi(s)` (default 85% of fibers cross),
nudged laterally by `center_gap` (mm, so the two copies interlace at the
decussation) and pushed posteriorly by `cross_detour` (mm, the
interdigitating X), without any streamline splitting (D1) or fusion
(C3).

## Examples

``` r
pyt <- simulate_pyramidal_tract(hemisphere = "left", seed = 1)
fiber::is_bundle_set(pyt)      # TRUE
#> [1] TRUE
```
