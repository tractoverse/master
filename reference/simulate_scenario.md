# Simulate a fascicle scenario as a `fiber::bundle_set`

Generates one of the master fascicle configurations and returns it as a
[`fiber::bundle_set`](https://tractoverse.github.io/fiber/reference/bundle_set.html)
ready for
[`rtists::plot3d()`](https://tractoverse.github.io/rtists/reference/plot3d.html).
Supported configurations:

## Usage

``` r
simulate_scenario(
  configuration = "bending",
  n_axons = 30,
  angle = 0.5,
  arm_length = 1,
  interaction_length = 2,
  variance = 0.02,
  autocorrelation = 2,
  n_basis = 8,
  n_points = 300L,
  separation = 0.3,
  seed = NULL
)
```

## Arguments

- configuration:

  One of `"bending"`, `"fanning"`, `"merging"`, `"kissing"`,
  `"crossing"`.

- n_axons:

  Number of streamlines (per bundle for the two-bundle configurations,
  or total for `"bending"` and `"fanning"`).

- angle:

  Bending angle (radians), shared geometry parameter.

- arm_length:

  Length of each bending arm, shared geometry parameter.

- interaction_length:

  Flat interaction length, shared geometry parameter.

- variance:

  Coefficient variance controlling bundle diameter.

- autocorrelation:

  Correlation rate (space between independent and coherent streamlines).

- n_basis:

  Number of B-spline coefficients per coordinate.

- n_points:

  Number of points per centreline after arc-length reparametrization.

- separation:

  Separation between mean pathways for the two-bundle configurations
  (`"merging"`, `"kissing"`, `"crossing"`).

- seed:

  Optional integer seed for full reproducibility.

## Value

A
[`fiber::bundle_set`](https://tractoverse.github.io/fiber/reference/bundle_set.html)
with one (bending, fanning) or two (merging, kissing, crossing) named
bundles.

## Details

- `"bending"`: a single bundle whose mean pathway bends
  ([`generate_centerline()`](https://tractoverse.github.io/master/reference/generate_centerline.md)).

- `"fanning"`: a single bundle that diverges into two branches after a
  fork; half the axons follow one branch, half the other (no axon
  splits).

- `"merging"`: two bundles that start apart and run together.

- `"kissing"`: two bundles that meet, travel together, then part.

- `"crossing"`: two bundles whose mean pathways touch at a point and
  continue; the streamlines interlace without fusion.

## Examples

``` r
st <- simulate_scenario("crossing", n_axons = 10, seed = 1)
fiber::is_bundle_set(st)      # TRUE
#> [1] TRUE
st@n_bundles                  # 2
#> [1] 2
```
