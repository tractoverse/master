# Smooth C1 interpolation of control points over a station grid

Interpolates control values given at the station abscissae `stations`
onto the grid `s` using the smoothstep ramp `t^2 (3 - 2 t)` (C1: value
and slope agree at every station). Values before the first station equal
the first control value; values after the last equal the last. This is
the same ramp used by
[`generate_centerline()`](https://tractoverse.github.io/master/reference/generate_centerline.md),
extracted so station tables can drive any smoothly-varying field (twist
`psi`, widths `W`, or per-parcel positions).

## Usage

``` r
smoothstep_interp(s, stations, vals)
```

## Arguments

- s:

  Numeric grid of abscissae (typically arc-length fractions in
  `[0, 1]`).

- stations:

  Numeric station abscissae, strictly increasing, with `range(stations)`
  covering `range(s)` recommended.

- vals:

  Numeric control values, one per station (same length as `stations`).

## Value

Numeric vector, same length as `s`.

## Examples

``` r
s <- seq(0, 1, length.out = 101)
psi <- smoothstep_interp(s, c(0, 0.5, 1), c(0, pi / 2, pi))
psi[1]        # 0
#> [1] 0
psi[101]      # pi
#> [1] 3.141593
```
