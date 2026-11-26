# Default pyramidal-tract parameter table

Returns the full D8 parameter table used by
[`simulate_pyramidal_tract()`](https://tractoverse.github.io/master/reference/simulate_pyramidal_tract.md)
when its arguments are `NULL`: the anatomical stations, the five
somatotopic parcels (`leg`, `trunk`, `arm`, `hand`, `face`), their
cross-section control points `q_st` (along `u`) and `r_st` (along `v`)
at each station, the section widths `width_u`/`width_v` (mm), the twist
map `twist`, per-parcel terminal stations (`face` ends at the `pons`,
everything else at the `cord`), calibers (axonal radius, mm), the
macro-axis `macros`, the decussation `cross` controls and the per-parcel
FDA `perturbation` controls (reused by
[`create_axons()`](https://tractoverse.github.io/master/reference/create_axons.md)).

## Usage

``` r
pyt_defaults()
```

## Value

A named list with elements `tract_length`, `stations`, `parcels`,
`q_st`, `r_st`, `width_u`, `width_v`, `twist`, `term_station`,
`caliber`, `macros`, `cross`, `perturbation`, as documented in DECISIONS
D8.

## Details

Coordinates follow the amended validation convention: at the cortex,
positive `u` points toward the medial wall of the parcel's hemisphere
and parcels are spaced by their relative sheet extents
(`leg -1.0, trunk -0.6, arm -0.15, hand +0.35, face +0.85`).

## Examples

``` r
d <- pyt_defaults()
d$parcels                       # five body-part parcels
#> [1] "leg"   "trunk" "arm"   "hand"  "face" 
d$term_station["face"]          # "pons"
#>   face 
#> "pons" 
```
