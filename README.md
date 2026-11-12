

<!-- README.md is generated from README.qmd. Please edit that file -->

# master <img src="man/figures/logo.png" align="right" height="138" alt="master package logo" />

<!-- badges: start -->

<!-- badges: end -->

`master` — the **Multivariate Axonal Simulator for Tractography
Evaluation** — generates realistic white-matter bundles (fascicles) with
microstructural attributes for evaluating diffusion-MRI tractography. It
implements the functional-data-analysis generative design of the MASTER
model: a reference centerline pathway (a flat central interaction
segment with smooth bending arms, assembled by tangent-angle integration
and reparametrized by curvilinear arc length) is perturbed by zero-mean
functional random fields of axons expressed in a B-spline basis with
exponentially autocorrelated Gaussian coefficients.

`master` builds on the `fiber` package (S7 `streamline` / `bundle` /
`bundle_set` classes) and emits `fiber::bundle` objects ready for
interactive visualization with `rtists`.

## Installation

You can install the development version of `master` from GitHub:

``` r
# install.packages("pak")
pak::pak("tractoverse/master")
```

Or from a local checkout:

``` r
pak::pak("local::/path/to/master")
```

## Example

Generate a single smooth centerline, then populate it with a bundle of
perturbed axons:

``` r
#| label: example
library(master)

# A straight central segment with two smoothly bending arms
cl <- generate_centerline(angle = 0.6, arm_length = 1,
                          interaction_length = 2, n_points = 200)
head(cl)          # N x 3 matrix, columns X, Y, Z

# A bundle of 20 axons coherently perturbed around that centerline
b <- create_axons(cl, n_axons = 20, n_basis = 8, variance = 0.02,
                  autocorrelation = 2, seed = 1)
fiber::is_bundle(b)   # TRUE
b@n_streamlines       # 20
```

Build a pre-defined interaction configuration as a `fiber::bundle_set`

``` r
#| label: scenario
st <- simulate_scenario("crossing", n_axons = 10, seed = 1)
fiber::is_bundle_set(st)   # TRUE
st@n_bundles               # 2

# Visualize with rtists (if installed)
# rtists::plot3d(st)

# Export to the master CSV layout for downstream tractography tools
generate_csv(st, tempfile(fileext = ".csv"))
```

## The generative model

Each bundle is grown from a reference centerline whose tangent angle
ramps smoothly from `-angle` to `+angle` with an exactly flat central
segment of length `interaction_length` (see \[generate_centerline()\]).
The centerline is reparametrized by curvilinear arc length, then each
axon is drawn as the centerline plus a smooth functional perturbation
modelled as a B-spline expansion of exponentially autocorrelated
Gaussian coefficients (see \[create_axons()\]). Pre-assembled
configurations — bending, fanning, merging, kissing, crossing — are
available through \[simulate_scenario()\].

## References

The MASTER design follows the functional-data generative approach
developed for diffusion-MRI tractography evaluation, building on the
`fiber` package and the `rtists` visualization stack.

You will still need to render `README.qmd` regularly, to keep
`README.md` up-to-date. `devtools::build_readme()` is handy for this.
