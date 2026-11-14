# Generate a pair of interacting centerlines

Builds two centreline mean-pathways for the multi-bundle configurations
(merging, kissing, crossing, fanning). Both are generated as single
bending centre-lines via
[`generate_centerline()`](https://tractoverse.github.io/master/reference/generate_centerline.md)
and then separated by `separation` along `Y` so they approach and
interact over the central region.

## Usage

``` r
generate_centerline_pair(
  angle = 0.5,
  arm_length = 1,
  interaction_length = 2,
  separation = 0.3,
  rotation = 0,
  n_points = 1000L
)
```

## Arguments

- angle:

  Bending angle (radians) of each centreline (see
  [`generate_centerline()`](https://tractoverse.github.io/master/reference/generate_centerline.md)).

- arm_length:

  Length of each bending arm (see
  [`generate_centerline()`](https://tractoverse.github.io/master/reference/generate_centerline.md)).

- interaction_length:

  Flat interaction length (see
  [`generate_centerline()`](https://tractoverse.github.io/master/reference/generate_centerline.md)).

- separation:

  Separation (units of `Y`) between the two mean-pathways at their
  closest approach. `separation = 0` with a zero interaction gives a
  pure crossing; small positive separations give merging/kissing.

- rotation:

  Optional rotation (radians) common to both centre-lines.

- n_points:

  Number of points per centreline after arc-length reparametrization.
  Defaults to `1000`.

## Value

A list of two point matrices (named `centerline1`, `centerline2`), each
a `n_points` x 3 matrix with columns `X`, `Y`, `Z`.

## Examples

``` r
cpair <- generate_centerline_pair(angle = 0.5, arm_length = 1,
                                  interaction_length = 2, separation = 0.3)
length(cpair)      # 2
#> [1] 2
nrow(cpair[[1]])   # 1000
#> [1] 1000
```
