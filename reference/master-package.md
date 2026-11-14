# master: a Multivariate Axonal Simulator for Tractography Evaluation

`master` is a generative simulator for white-matter fascicles used to
build and evaluate tractography algorithms. It implements the FDA-based
generative design:

## Details

1.  a reference **centerline** (functional mean pathway) is assembled
    from a central straight interaction segment with smooth bending arms
    ([`generate_centerline()`](https://tractoverse.github.io/master/reference/generate_centerline.md));

2.  the centerline is discretized onto an **equally spaced grid in
    curvilinear arc length**
    ([`reparametrize_3D_curve()`](https://tractoverse.github.io/master/reference/reparametrize_3D_curve.md));

3.  individual **axons** are drawn as the centerline plus a zero-mean
    functional random field: each coordinate is expanded in a B-spline
    basis whose Gaussian coefficients carry an *exponential
    autocorrelation*
    ([`create_axons()`](https://tractoverse.github.io/master/reference/create_axons.md)).
    Zero mean keeps the sample mean of the bundle on the true
    centerline.

The emitted objects use the `fiber` model
([`fiber::streamline`](https://tractoverse.github.io/fiber/reference/streamline.html),
[`fiber::bundle`](https://tractoverse.github.io/fiber/reference/bundle.html),
[`fiber::bundle_set`](https://tractoverse.github.io/fiber/reference/bundle_set.html))
and drop straight into
[`rtists::plot3d()`](https://tractoverse.github.io/rtists/reference/plot3d.html)
for interactive visualization.

## See also

Useful links:

- <https://github.com/tractoverse/master>

- <https://tractoverse.github.io/master/>

- Report bugs at <https://github.com/tractoverse/master/issues>

## Author

**Maintainer**: Aymeric Stamm <aymeric.stamm@cnrs.fr>
([ORCID](https://orcid.org/0000-0002-8725-3654)) \[copyright holder\]

Authors:

- Aymeric Stamm <aymeric.stamm@cnrs.fr>
  ([ORCID](https://orcid.org/0000-0002-8725-3654)) \[copyright holder\]
