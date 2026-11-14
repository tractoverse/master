# Generate a bundle of axons around a centreline

Draws `n_axons` streamlines (axons) as the centreline plus a zero-mean
functional perturbation, and returns them as a
[`fiber::bundle`](https://tractoverse.github.io/fiber/reference/bundle.html).

## Usage

``` r
create_axons(
  centerline,
  n_axons = 30,
  n_basis = 8,
  variance = 0.02,
  autocorrelation = 2,
  radius = NULL,
  seed = NULL
)
```

## Arguments

- centerline:

  A numeric `N x 3` matrix (columns `X`, `Y`, `Z`) giving the mean
  pathway; see
  [`generate_centerline()`](https://tractoverse.github.io/master/reference/generate_centerline.md).

- n_axons:

  Number of streamlines (axons) in the bundle.

- n_basis:

  Number of B-spline coefficients per coordinate (smoothness of the
  perturbation).

- variance:

  Coefficient variance `sigma^2` controlling the bundle diameter.

- autocorrelation:

  Positive correlation rate `> 0` (space between independent and
  perfectly coherent streamlines).

- radius:

  Optional axonal radius (scalar or length-`N` vector).

- seed:

  Optional integer seed for reproducibility of the draws.

## Value

A
[`fiber::bundle`](https://tractoverse.github.io/fiber/reference/bundle.html)
of `n_axons` streamlines.

## Details

The perturbation of each coordinate is a linear combination of `n_basis`
B-spline functions over the centreline abscissa with zero-mean Gaussian
coefficients carrying an exponential autocorrelation. Because the field
is zero-mean, the sample mean of the bundle stays on the true
centreline.

Each axon carries per-point microstructure attributes as `@point_data`:
`curvature`, `torsion` and the unit tangent components `dir_x`, `dir_y`,
`dir_z`, computed from that axon's own geometry. If `radius` is a scalar
it is stored per-streamline; if it is a length-`N` vector it is stored
as a per-point `radius`.

## Examples

``` r
cl <- generate_centerline(angle = 0.5, arm_length = 1,
                          interaction_length = 2, n_points = 100)
b <- create_axons(cl, n_axons = 20, n_basis = 8, variance = 0.02,
                  autocorrelation = 2, seed = 1)
fiber::is_bundle(b)      # TRUE
#> [1] TRUE
b@n_streamlines          # 20
#> [1] 20
```
