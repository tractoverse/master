# Reparametrize a point sequence onto uniform curvilinear arc length

Thin wrapper around
[`fiber::reparametrize()`](https://tractoverse.github.io/fiber/reference/reparametrize.html)
that takes a raw `n x 3` point matrix (columns `X`, `Y`, `Z`), resamples
it onto `n_points` equally spaced points in arc length, and returns the
new point matrix.

## Usage

``` r
reparametrize_3D_curve(points, n_points = NULL)
```

## Arguments

- points:

  A numeric matrix with columns `X`, `Y`, `Z`.

- n_points:

  Target number of equally spaced points. Defaults to the number of rows
  in `points`.

## Value

A numeric matrix with columns `X`, `Y`, `Z` and `n_points` rows.

## Examples

``` r
m <- cbind(X = 0:10, Y = (0:10)^2 / 10, Z = 0)
m2 <- reparametrize_3D_curve(m, n_points = 20)
nrow(m2) # 20
#> [1] 20
```
