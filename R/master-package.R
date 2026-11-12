#' master: a Multivariate Axonal Simulator for Tractography Evaluation
#'
#' `master` is a generative simulator for white-matter fascicles used to build
#' and evaluate tractography algorithms. It implements the FDA-based generative design:
#'
#' 1. a reference **centerline** (functional mean pathway) is assembled from a
#'    central straight interaction segment with smooth bending arms
#'    ([generate_centerline()]);
#' 2. the centerline is discretized onto an **equally spaced grid in
#'    curvilinear arc length** ([reparametrize_3D_curve()]);
#' 3. individual **axons** are drawn as the centerline plus a zero-mean
#'    functional random field: each coordinate is expanded in a B-spline basis
#'    whose Gaussian coefficients carry an *exponential autocorrelation*
#'    ([create_axons()]). Zero mean keeps the sample mean of the bundle on the
#'    true centerline.
#'
#' The emitted objects use the `fiber` model (`fiber::streamline`,
#' `fiber::bundle`, `fiber::bundle_set`) and drop straight into
#' `rtists::plot3d()` for interactive visualization.
#'
#' @keywords internal
"_PACKAGE"
