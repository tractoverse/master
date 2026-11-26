# Anatomical stations of the pyramidal tract (cortex to cord)

Ordered stations along the tract, as a fraction of the total arc length
(`fraction`, in `[0, 1]`) and as a depth in millimetres measured from
the motor-cortex surface along the pathway (`depth_mm`, at the default
`tract_length = 75` mm).

## Usage

``` r
pyt_stations()
```

## Value

A data frame with one row per station and columns `station` (ordered
factor), `fraction`, `depth_mm` (at the default `tract_length = 70` mm).

## Examples

``` r
pyt_stations()
#>   station fraction depth_mm
#> 1  cortex     0.00      0.0
#> 2  corona     0.10      7.0
#> 3 centrum     0.20     14.0
#> 4      ic     0.35     24.5
#> 5    crus     0.50     35.0
#> 6    pons     0.60     42.0
#> 7 pyramid     0.72     50.4
#> 8     dec     0.82     57.4
#> 9    cord     1.00     70.0
```
