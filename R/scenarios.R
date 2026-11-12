# ---- fascicle configurations ------------------------------------------------
#
# master models fascicle geometry from a small vocabulary of configurations.
# Two rules hold across all of them:
#   * a streamline never splits in two: each axon is one continuous single
#     path from end to end (configurations that look like "splitting" are
#     really fanning);
#   * at a crossing streamlines interlace but no streamline fuses with or
#     enters another.
# All multi-bundle configurations are produced from a single centreline-pair
# generator plus one axon generator.

#' Simulate a fascicle scenario as a `fiber::bundle_set`
#'
#' Generates one of the master fascicle configurations and returns it as a
#' `fiber::bundle_set` ready for [rtists::plot3d()]. Supported configurations:
#'
#' - `"bending"`: a single bundle whose mean pathway bends ([generate_centerline()]).
#' - `"fanning"`: a single bundle that diverges into two branches after a fork;
#'   half the axons follow one branch, half the other (no axon splits).
#' - `"merging"`: two bundles that start apart and run together.
#' - `"kissing"`: two bundles that meet, travel together, then part.
#' - `"crossing"`: two bundles whose mean pathways touch at a point and
#'   continue; the streamlines interlace without fusion.
#'
#' @param configuration One of `"bending"`, `"fanning"`, `"merging"`,
#'   `"kissing"`, `"crossing"`.
#' @param n_axons Number of streamlines (per bundle for the two-bundle
#'   configurations, or total for `"bending"` and `"fanning"`).
#' @param angle Bending angle (radians), shared geometry parameter.
#' @param arm_length Length of each bending arm, shared geometry parameter.
#' @param interaction_length Flat interaction length, shared geometry
#'   parameter.
#' @param variance Coefficient variance controlling bundle diameter.
#' @param autocorrelation Correlation rate (space between independent and
#'   coherent streamlines).
#' @param n_basis Number of B-spline coefficients per coordinate.
#' @param n_points Number of points per centreline after arc-length
#'   reparametrization.
#' @param separation Separation between mean pathways for the two-bundle
#'   configurations (`"merging"`, `"kissing"`, `"crossing"`).
#' @param seed Optional integer seed for full reproducibility.
#' @returns A `fiber::bundle_set` with one (bending, fanning) or two
#'   (merging, kissing, crossing) named bundles.
#' @export
#' @examples
#' st <- simulate_scenario("crossing", n_axons = 10, seed = 1)
#' fiber::is_bundle_set(st)      # TRUE
#' st@n_bundles                  # 2
simulate_scenario <- function(
    configuration = "bending",
    n_axons = 30,
    angle = 0.5,
    arm_length = 1,
    interaction_length = 2,
    variance = 0.02,
    autocorrelation = 2,
    n_basis = 8,
    n_points = 300L,
    separation = 0.3,
    seed = NULL
) {
  configuration <- match.arg(
    configuration,
    c("bending", "fanning", "merging", "kissing", "crossing")
  )
  if (!is.null(seed)) {
    set.seed(seed)
  }

  one_bundle <- function(cl, name, n_ax) {
    fiber::bundle(
      streamlines = create_axons(
        cl, n_axons = n_ax, n_basis = n_basis, variance = variance,
        autocorrelation = autocorrelation
      )@streamlines,
      bundle_data = list(role = name)
    )
  }

  if (configuration == "bending") {
    cl <- generate_centerline(angle, arm_length, interaction_length, n_points = n_points)
    b <- one_bundle(cl, "bending", n_axons)
    return(fiber::bundle_set(bundles = list("bending" = b)))
  }

  if (configuration == "fanning") {
    # Two branches from a common flat trunk: branch A bends up, branch B bends
    # down; half the axons follow each branch.
    n_half <- max(1L, floor(n_axons / 2))
    bA <- generate_centerline(angle, arm_length, interaction_length, n_points = n_points)
    bB <- generate_centerline(angle, arm_length, interaction_length, n_points = n_points)
    # Mirror branch B across the flat axis so it fans the other way.
    bB[, "Y"] <- -bB[, "Y"]
    slA <- create_axons(bA, n_axons = n_half, n_basis = n_basis,
                        variance = variance, autocorrelation = autocorrelation)@streamlines
    slB <- create_axons(bB, n_axons = n_half, n_basis = n_basis,
                        variance = variance, autocorrelation = autocorrelation)@streamlines
    b <- fiber::bundle(
      streamlines = c(slA, slB),
      bundle_data = list(role = "fanning")
    )
    return(fiber::bundle_set(bundles = list("fanning" = b)))
  }

  # Two-bundle configurations share one code path. A pure crossing has the
  # mean pathways touching (separation 0); merging/kissing keep a gap.
  eff_sep <- if (configuration == "crossing") 0 else separation
  cpair <- generate_centerline_pair(
    angle, arm_length, interaction_length,
    separation = eff_sep, n_points = n_points
  )
  b1 <- one_bundle(cpair$centerline1, configuration, n_axons)
  b2 <- one_bundle(cpair$centerline2, configuration, n_axons)
  fiber::bundle_set(bundles = list("bundle1" = b1, "bundle2" = b2))
}
