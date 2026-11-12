# ---- centerline generation --------------------------------------------------
#
# A centreline is a smooth reference mean-pathway assembled by integrating a
# prescribed tangent-angle field over arc length (the PLAN's tangent-angle
# approach). Walking along +x from u = 0, the tangent angle theta(u) is:
#
#   * left arm  [0, A]        : theta ramps -alpha -> 0   (smooth step)
#   * central   [A, A+L]      : theta stays 0             (exactly flat, length L)
#   * right arm [A+L, A+L+A]  : theta ramps  0  -> +alpha (smooth step)
#
# Integrating (cos theta, sin theta) over u yields a symmetric bow: a flat
# central interaction segment of exact length `interaction_length` whose two
# arms bend smoothly (C1, singularity-free) up to +-angle. Because each
# integrand step has unit length, the result is already parametrized by arc
# length, matching [PLAN Step II] directly.

# Mirror smoothstep ramp R(t) = t^2 (3 - 2 t), t in [0, 1]. C1: value and
# slope agree at both endpoints (R'(0) = R'(1) = 0).
.smoothstep <- function(t) t^2 * (3 - 2 * t)

#' Generate a smooth centerline with a flat interaction segment
#'
#' Assembles a single reference mean-pathway by integrating a tangent-angle
#' field over arc length. The central interaction segment is
#' exactly straight for `interaction_length` units; the two equal arms bend
#' smoothly (C1) by `angle` (radians). The result is already parametrized by
#' curvilinear arc length and is discretized onto `n_points` equally spaced
#' points (via [fiber::reparametrize()]).
#'
#' @param angle Bending angle (radians) reached by each arm. Bounded to
#'   `(-pi/2, pi/2)`.
#' @param arm_length Length (in output units) of each straight run of the
#'   bending arms that lead up to the flat central segment. Together with
#'   `interaction_length` it sets how far the arms extend.
#' @param interaction_length Length of the flat central (straight) part of the
#'   centreline, in the same units as the output points.
#' @param rotation Rotation (radians) applied to the whole centreline about
#'   the vertical (`Y`) axis after generation. Defaults to `0`.
#' @param vertical_bend Optional bend angle (radians) applied out of the
#'   `xy`-plane (towards `z`), so the centreline can tile a truly 3-D
#'   fascicle. Defaults to `0` (a planar centreline).
#' @param n_points Number of equally spaced points on the grid (in
#'   curvilinear arc length). Defaults to `1000`.
#' @returns A `n_points` x 3 numeric matrix with column names `X`, `Y`, `Z`
#'   giving the arc-length-parametrized centreline.
#' @export
#' @examples
#' cl <- generate_centerline(angle = 0.6, arm_length = 1,
#'                           interaction_length = 2, n_points = 200)
#' colnames(cl) # c("X", "Y", "Z")
#' nrow(cl)     # 200
generate_centerline <- function(
    angle = 0.5,
    arm_length = 1,
    interaction_length = 2,
    rotation = 0,
    vertical_bend = 0,
    n_points = 1000L
) {
  if (!(n_points >= 2L)) {
    cli::cli_abort("{.arg n_points} must be at least 2.")
  }
  if (interaction_length <= 0) {
    cli::cli_abort("{.arg interaction_length} must be strictly positive.")
  }
  if (arm_length <= 0) {
    cli::cli_abort("{.arg arm_length} must be strictly positive.")
  }
  angle <- as.numeric(angle)
  if (abs(angle) >= pi / 2) {
    cli::cli_abort("{.arg angle} must lie inside (-pi/2, pi/2).")
  }
  arm_length <- as.numeric(arm_length)
  interaction_length <- as.numeric(interaction_length)

  A <- arm_length
  L <- interaction_length
  T <- L + 2 * A            # total arc length
  alpha <- angle

  # Dense integration over arc length on a uniform grid.
  N <- 6000L
  u <- seq(0, T, length.out = N)
  du <- T / (N - 1)

  theta <- numeric(N)
  left <- u <= A
  right <- u > A + L
  central <- !left & !right
  t_left <- u[left] / A
  t_right <- (u[right] - (A + L)) / A
  theta[left] <- -alpha * (1 - .smoothstep(t_left))
  theta[right] <- alpha * .smoothstep(t_right)
  # central stays theta = 0 -> exactly straight.

  # Integrate unit tangent: points are already arc-length-parametrized.
  x <- cumsum(cos(theta) * du)
  y <- cumsum(sin(theta) * du)

  # Vertical (out-of-plane) bend: linear ramp over z so the centreline stays C1.
  z <- if (vertical_bend == 0) {
    rep(0, N)
  } else {
    tan(vertical_bend) * u / T
  }

  pts <- cbind(X = x, Y = y, Z = z)

  # Apply in-plane rotation about the vertical (Y = up) axis.
  if (rotation != 0) {
    c_r <- cos(rotation)
    s_r <- sin(rotation)
    x_new <- pts[, "X"] * c_r + pts[, "Z"] * s_r
    z_new <- -pts[, "X"] * s_r + pts[, "Z"] * c_r
    pts[, "X"] <- x_new
    pts[, "Z"] <- z_new
  }

  # Discretize onto a uniform grid in curvilinear arc length.
  sl <- fiber::streamline(points = pts)
  out <- fiber::reparametrize(sl, n_points = n_points)@points
  rownames(out) <- NULL
  out
}

#' Generate a pair of interacting centerlines
#'
#' Builds two centreline mean-pathways for the multi-bundle configurations
#' (merging, kissing, crossing, fanning). Both are generated as single bending
#' centre-lines via [generate_centerline()] and then separated by `separation`
#' along `Y` so they approach and interact over the central region.
#'
#' @param angle Bending angle (radians) of each centreline (see
#'   [generate_centerline()]).
#' @param arm_length Length of each bending arm (see
#'   [generate_centerline()]).
#' @param interaction_length Flat interaction length (see
#'   [generate_centerline()]).
#' @param separation Separation (units of `Y`) between the two mean-pathways
#'   at their closest approach. `separation = 0` with a zero interaction gives
#'   a pure crossing; small positive separations give merging/kissing.
#' @param rotation Optional rotation (radians) common to both centre-lines.
#' @param n_points Number of points per centreline after arc-length
#'   reparametrization. Defaults to `1000`.
#' @returns A list of two point matrices (named `centerline1`, `centerline2`),
#'   each a `n_points` x 3 matrix with columns `X`, `Y`, `Z`.
#' @export
#' @examples
#' cpair <- generate_centerline_pair(angle = 0.5, arm_length = 1,
#'                                   interaction_length = 2, separation = 0.3)
#' length(cpair)      # 2
#' nrow(cpair[[1]])   # 1000
generate_centerline_pair <- function(
    angle = 0.5,
    arm_length = 1,
    interaction_length = 2,
    separation = 0.3,
    rotation = 0,
    n_points = 1000L
) {
  cl1 <- generate_centerline(
    angle = angle, arm_length = arm_length,
    interaction_length = interaction_length,
    rotation = rotation, n_points = n_points
  )
  cl2 <- generate_centerline(
    angle = angle, arm_length = arm_length,
    interaction_length = interaction_length,
    rotation = rotation, n_points = n_points
  )
  cl2[, "Y"] <- cl2[, "Y"] + separation
  list(centerline1 = cl1, centerline2 = cl2)
}

#' Reparametrize a point sequence onto uniform curvilinear arc length
#'
#' Thin wrapper around [fiber::reparametrize()] that takes a raw `n x 3` point
#' matrix (columns `X`, `Y`, `Z`), resamples it onto `n_points` equally spaced
#' points in arc length, and returns the new point matrix.
#'
#' @param points A numeric matrix with columns `X`, `Y`, `Z`.
#' @param n_points Target number of equally spaced points. Defaults to the
#'   number of rows in `points`.
#' @returns A numeric matrix with columns `X`, `Y`, `Z` and `n_points` rows.
#' @export
#' @examples
#' m <- cbind(X = 0:10, Y = (0:10)^2 / 10, Z = 0)
#' m2 <- reparametrize_3D_curve(m, n_points = 20)
#' nrow(m2) # 20
reparametrize_3D_curve <- function(points, n_points = NULL) {
  sl <- fiber::streamline(points = points)
  out <- fiber::reparametrize(sl, n_points = n_points)@points
  rownames(out) <- NULL
  out
}
