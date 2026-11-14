# Build the exponential-autocorrelation covariance for B-spline coefficients

Returns the `n_basis` x `n_basis` zero-mean Gaussian covariance with
entries `Sigma[i, j] = variance * exp(-|i - j| / autocorrelation)`. This
is the covariance of the functional random field in the B-spline basis.

## Usage

``` r
.bspline_covariance(n_basis, variance, autocorrelation)
```

## Arguments

- n_basis:

  Number of B-spline coefficients (basis functions).

- variance:

  Variance `sigma^2` (amplitude) shared by all coefficients.

- autocorrelation:

  Positive correlation rate `rate > 0` (`rate -> 0` gives independent
  coefficients; larger values give smoother, more coherent streamlines).

## Value

An `n_basis` x `n_basis` numeric matrix.
