# Draw zero-mean Gaussian B-spline coefficients

Draws one coefficient vector per coordinate (X, Y, Z) per axon from a
zero-mean multivariate normal with covariance
[`.bspline_covariance()`](https://tractoverse.github.io/master/reference/dot-bspline_covariance.md).
Returns one `n_basis` x `n_axons` matrix per coordinate.

## Usage

``` r
.draw_coefficients(n_basis, variance, autocorrelation, n_axons, seed = NULL)
```

## Arguments

- n_basis:

  Number of B-spline coefficients per coordinate.

- variance:

  Coefficient variance.

- autocorrelation:

  Correlation rate.

- n_axons:

  Number of independent draws.

- seed:

  Optional integer seed for reproducibility.

## Value

A list with elements `X`, `Y`, `Z`, each an `n_basis` x `n_axons`
matrix.
