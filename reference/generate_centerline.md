# Generate a smooth centerline with a flat interaction segment

Assembles a single reference mean-pathway by integrating a tangent-angle
field over arc length. The central interaction segment is exactly
straight for `interaction_length` units; the two equal arms bend
smoothly (C1) by `angle` (radians). The result is already parametrized
by curvilinear arc length and is discretized onto `n_points` equally
spaced points (via
[`fiber::reparametrize()`](https://tractoverse.github.io/fiber/reference/reparametrize.html)).

## Usage

``` r
generate_centerline(
  angle = 0.5,
  arm_length = 1,
  interaction_length = 2,
  rotation = 0,
  vertical_bend = 0,
  n_points = 1000L
)
```

## Arguments

- angle:

  Bending angle (radians) reached by each arm. Bounded to
  `(-pi/2, pi/2)`.

- arm_length:

  Length (in output units) of each straight run of the bending arms that
  lead up to the flat central segment. Together with
  `interaction_length` it sets how far the arms extend.

- interaction_length:

  Length of the flat central (straight) part of the centreline, in the
  same units as the output points.

- rotation:

  Rotation (radians) applied to the whole centreline about the vertical
  (`Y`) axis after generation. Defaults to `0`.

- vertical_bend:

  Optional bend angle (radians) applied out of the `xy`-plane (towards
  `z`), so the centreline can tile a truly 3-D fascicle. Defaults to `0`
  (a planar centreline).

- n_points:

  Number of equally spaced points on the grid (in curvilinear arc
  length). Defaults to `1000`.

## Value

A `n_points` x 3 numeric matrix with column names `X`, `Y`, `Z` giving
the arc-length-parametrized centreline.

## Examples

``` r
cl <- generate_centerline(angle = 0.6, arm_length = 1,
                          interaction_length = 2, n_points = 200)
colnames(cl) # c("X", "Y", "Z")
#> [1] "X" "Y" "Z"
nrow(cl)     # 200
#> [1] 200
```
