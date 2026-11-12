# ---- axon (streamline) generation ------------------------------------------
#
# [PLAN Step III] Each axon is the centreline plus a zero-mean 3-D
# perturbation, so the sample mean of the bundle stays on the true centreline
# (the ground-truth mean pathway). Each coordinate of the perturbation is
# written in an order-`degree + 1` B-spline basis over the centreline
# abscissa, and the random coefficients are draws from a zero-mean multivariate
# Gaussian carrying an exponential autocorrelation:
#     Sigma[i, j] = variance * exp(-|i - j| / autocorrelation)
# `autocorrelation -> 0` gives independent points (free to mix); larger values
# produce more parallel, coherent streamlines. Zero mean is what guarantees the
# centreline remains the exact ground truth.

#' @importFrom stats rnorm
NULL

#' Build the exponential-autocorrelation covariance for B-spline coefficients
#'
#' Returns the `n_basis` x `n_basis` zero-mean Gaussian covariance with entries
#' `Sigma[i, j] = variance * exp(-|i - j| / autocorrelation)`. This is the
#' covariance of the functional random field in the B-spline basis.
#'
#' @param n_basis Number of B-spline coefficients (basis functions).
#' @param variance Variance `sigma^2` (amplitude) shared by all coefficients.
#' @param autocorrelation Positive correlation rate `rate > 0` (`rate -> 0`
#'   gives independent coefficients; larger values give smoother, more
#'   coherent streamlines).
#' @returns An `n_basis` x `n_basis` numeric matrix.
#' @keywords internal
.bspline_covariance <- function(n_basis, variance, autocorrelation) {
  r <- seq_len(n_basis)
  variance * exp(-abs(outer(r, r, "-")) / autocorrelation)
}

#' Draw zero-mean Gaussian B-spline coefficients
#'
#' Draws one coefficient vector per coordinate (X, Y, Z) per axon from a
#' zero-mean multivariate normal with covariance [.bspline_covariance()].
#' Returns one `n_basis` x `n_axons` matrix per coordinate.
#'
#' @param n_basis Number of B-spline coefficients per coordinate.
#' @param variance Coefficient variance.
#' @param autocorrelation Correlation rate.
#' @param n_axons Number of independent draws.
#' @param seed Optional integer seed for reproducibility.
#' @returns A list with elements `X`, `Y`, `Z`, each an `n_basis` x `n_axons`
#'   matrix.
#' @keywords internal
.draw_coefficients <- function(n_basis, variance, autocorrelation, n_axons, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  Sig <- .bspline_covariance(n_basis, variance, autocorrelation)
  L <- chol(Sig)
  draw <- function() {
    Z <- matrix(rnorm(n_basis * n_axons), n_basis, n_axons)
    crossprod(L, Z)
  }
  list(X = draw(), Y = draw(), Z = draw())
}

# Vectorised normalisation of rows to unit length (safe against zero rows).
.row_unit <- function(M) {
  n <- sqrt(rowSums(M * M))
  M / pmax(n, 1e-12)
}

# Three-component vector (cross) product, row-wise.
.cross <- function(u, v) {
  cbind(
    u[, 2] * v[, 3] - u[, 3] * v[, 2],
    u[, 3] * v[, 1] - u[, 1] * v[, 3],
    u[, 1] * v[, 2] - u[, 2] * v[, 1]
  )
}

#' Per-point shape descriptors (orientation, curvature, torsion) of a curve
#'
#' Estimates the unit tangent (orientation), local curvature and torsion at
#' each point of a 3-D curve from its finite differences. Returns scalar
#' per-point attributes suitable for `fiber::streamline@point_data` (each must
#' be a numeric vector of length `nrow(pts)`), so the 3-D orientation is stored
#' as its three components `dir_x`, `dir_y`, `dir_z`.
#'
#' @param pts A numeric matrix with columns `X`, `Y`, `Z`.
#' @returns A list with numeric vectors `curvature`, `torsion`, `dir_x`,
#'   `dir_y`, `dir_z`, each of length `nrow(pts)`.
#' @keywords internal
.point_shape <- function(pts) {
  n <- nrow(pts)
  d <- diff(pts)
  seg <- sqrt(rowSums(d * d))              # length n-1: per-segment arc length
  t_seg <- .row_unit(d)                    # unit tangent per segment
  # Average adjacent segment tangents onto the shared grid points.
  t_mid <- (rbind(t_seg[1, , drop = FALSE], t_seg) +
             rbind(t_seg, t_seg[n - 1, , drop = FALSE])) / 2
  t_pt <- .row_unit(t_mid)                 # unit tangent at each point

  # Curvature at each interface = angle between adjacent point tangents /
  # the travelled arc length (the shared segment).
  ang <- acos(pmax(pmin(rowSums(t_pt[-n, , drop = FALSE] * t_pt[-1, , drop = FALSE]), 1), -1))
  curv_if <- ang / pmax(seg, 1e-12)        # length n-1
  curv <- c(
    curv_if[1],
    (curv_if[-1] + curv_if[-length(curv_if)]) / 2,
    curv_if[length(curv_if)]
  )

  # Frenet frame from the turning of the tangent: normal is the component of
  # acceleration perpendicular to the tangent; binormal = T x N.
  acc <- cbind(
    pts[c(1, 1:(n - 1)), "X"] - 2 * pts[, "X"] + pts[c(2:n, n), "X"],
    pts[c(1, 1:(n - 1)), "Y"] - 2 * pts[, "Y"] + pts[c(2:n, n), "Y"],
    pts[c(1, 1:(n - 1)), "Z"] - 2 * pts[, "Z"] + pts[c(2:n, n), "Z"]
  )
  normal <- .row_unit(acc - rowSums(acc * t_pt) * t_pt)
  binormal <- .row_unit(.cross(t_pt, normal))

  # Torsion = -dot(dB/ds, N), from adjacent binormals over the shared segment.
  db <- diff(binormal)
  mid_normal <- (normal[-n, , drop = FALSE] + normal[-1, , drop = FALSE]) / 2
  tor_if <- -rowSums(db * mid_normal) / pmax(seg, 1e-12)
  tor <- c(
    tor_if[1],
    (tor_if[-1] + tor_if[-length(tor_if)]) / 2,
    tor_if[length(tor_if)]
  )

  list(
    curvature = as.numeric(curv),
    torsion = as.numeric(tor),
    dir_x = t_pt[, 1],
    dir_y = t_pt[, 2],
    dir_z = t_pt[, 3]
  )
}

#' Generate a bundle of axons around a centreline
#'
#' Draws `n_axons` streamlines (axons) as the centreline plus a zero-mean
#' functional perturbation, and returns them as a `fiber::bundle`.
#'
#' The perturbation of each coordinate is a linear combination of `n_basis`
#' B-spline functions over the centreline abscissa with zero-mean Gaussian
#' coefficients carrying an exponential autocorrelation. Because the field is
#' zero-mean, the sample mean of the bundle stays on the true centreline.
#'
#' Each axon carries per-point microstructure attributes as `@point_data`:
#' `curvature`, `torsion` and the unit tangent components `dir_x`, `dir_y`,
#' `dir_z`, computed from that axon's own geometry. If `radius` is a scalar it
#' is stored per-streamline; if it is a length-`N` vector it is stored as a
#' per-point `radius`.
#'
#' @param centerline A numeric `N x 3` matrix (columns `X`, `Y`, `Z`) giving
#'   the mean pathway; see [generate_centerline()].
#' @param n_axons Number of streamlines (axons) in the bundle.
#' @param n_basis Number of B-spline coefficients per coordinate (smoothness of
#'   the perturbation).
#' @param variance Coefficient variance `sigma^2` controlling the bundle
#'   diameter.
#' @param autocorrelation Positive correlation rate `> 0` (space between
#'   independent and perfectly coherent streamlines).
#' @param radius Optional axonal radius (scalar or length-`N` vector).
#' @param seed Optional integer seed for reproducibility of the draws.
#' @returns A `fiber::bundle` of `n_axons` streamlines.
#' @export
#' @examples
#' cl <- generate_centerline(angle = 0.5, arm_length = 1,
#'                           interaction_length = 2, n_points = 100)
#' b <- create_axons(cl, n_axons = 20, n_basis = 8, variance = 0.02,
#'                   autocorrelation = 2, seed = 1)
#' fiber::is_bundle(b)      # TRUE
#' b@n_streamlines          # 20
create_axons <- function(
    centerline,
    n_axons = 30,
    n_basis = 8,
    variance = 0.02,
    autocorrelation = 2,
    radius = NULL,
    seed = NULL
) {
  centerline <- as.matrix(centerline)
  n <- nrow(centerline)
  if (ncol(centerline) != 3L) {
    cli::cli_abort("{.arg centerline} must be an N x 3 point matrix.")
  }
  colnames(centerline) <- c("X", "Y", "Z")
  if (n_axons < 1L || n_basis < 1L || n_basis > n) {
    cli::cli_abort("{.arg n_basis} must be in [1, nrow(centerline)].")
  }

  # One B-spline basis over the (arc-length) abscissa, shared by all axons.
  sgrid <- seq(0, 1, length.out = n)
  B <- splines::bs(sgrid, df = n_basis, intercept = TRUE)

  coefs <- .draw_coefficients(n_basis, variance, autocorrelation, n_axons, seed)

  sls <- vector("list", n_axons)
  for (j in seq_len(n_axons)) {
    eps <- cbind(
      X = as.numeric(B %*% coefs$X[, j]),
      Y = as.numeric(B %*% coefs$Y[, j]),
      Z = as.numeric(B %*% coefs$Z[, j])
    )
    pts <- centerline + eps
    sh <- .point_shape(pts)

    point_data <- list(
      curvature = sh$curvature,
      torsion = sh$torsion,
      dir_x = sh$dir_x,
      dir_y = sh$dir_y,
      dir_z = sh$dir_z
    )
    streamline_data <- list(bundle_id = "bend")
    if (!is.null(radius)) {
      if (length(radius) == 1L) {
        streamline_data$radius <- radius
      } else if (length(radius) == n) {
        point_data$radius <- as.numeric(radius)
      } else {
        cli::cli_abort(c(
          "{.arg radius} must be a scalar or a vector of length nrow(centerline).",
          "i" = "Got length {length(radius)} for {n} centreline points."
        ))
      }
    }

    sls[[j]] <- fiber::streamline(
      points = pts,
      point_data = point_data,
      streamline_data = streamline_data
    )
  }

  fiber::bundle(streamlines = sls)
}
