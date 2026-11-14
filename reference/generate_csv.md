# Export a `fiber::bundle_set` (or `bundle`) to a VTP-to-CSV-compatible file

Writes simulated tractography data to a plain-text CSV in the column
layout used by the master VTP-to-CSV tools. Every streamline becomes one
line, with the file name, a running point counter and the three
coordinates (`X`, `Y`, `Z`) of its points, plus an optional per-point
attribute column.

## Usage

``` r
generate_csv(x, file, attribute = NULL)
```

## Arguments

- x:

  A
  [`fiber::bundle_set`](https://tractoverse.github.io/fiber/reference/bundle_set.html),
  [`fiber::bundle`](https://tractoverse.github.io/fiber/reference/bundle.html),
  or
  [`fiber::streamline`](https://tractoverse.github.io/fiber/reference/streamline.html).

- file:

  Path to the output `.csv`.

- attribute:

  Optional name of a per-point attribute (e.g. `"radius"`,
  `"curvature"`) to append as an extra trailing column.

## Value

The path `file`, invisibly.

## Examples

``` r
cl <- generate_centerline(angle = 0.5, n_points = 20)
b <- create_axons(cl, n_axons = 3)
out <- generate_csv(b, tempfile(fileext = ".csv"))
head(read.csv(out))
#>                   File Point           X            Y            Z Tract
#> 1 file19d917ebd8bb.csv     0 -0.04994454  0.139913893  0.096364501     1
#> 2 file19d917ebd8bb.csv     1  0.07180139  0.093720189  0.069045253     1
#> 3 file19d917ebd8bb.csv     2  0.23657527 -0.004234815  0.030889684     1
#> 4 file19d917ebd8bb.csv     3  0.43312837 -0.109063975 -0.004054042     1
#> 5 file19d917ebd8bb.csv     4  0.64499874 -0.187065846 -0.021771059     1
#> 6 file19d917ebd8bb.csv     5  0.86121282 -0.227519520 -0.015306025     2
```
