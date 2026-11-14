# Per-point shape descriptors (orientation, curvature, torsion) of a curve

Estimates the unit tangent (orientation), local curvature and torsion at
each point of a 3-D curve from its finite differences. Returns scalar
per-point attributes suitable for `fiber::streamline@point_data` (each
must be a numeric vector of length `nrow(pts)`), so the 3-D orientation
is stored as its three components `dir_x`, `dir_y`, `dir_z`.

## Usage

``` r
.point_shape(pts)
```

## Arguments

- pts:

  A numeric matrix with columns `X`, `Y`, `Z`.

## Value

A list with numeric vectors `curvature`, `torsion`, `dir_x`, `dir_y`,
`dir_z`, each of length `nrow(pts)`.
