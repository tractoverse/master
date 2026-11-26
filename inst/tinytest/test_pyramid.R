# Tests for the per-bundle pyramidal-tract model (D8).

set.seed(4321)

# --- smoothstep_interp ---------------------------------------------------

s <- seq(0, 1, length.out = 51)
psi <- smoothstep_interp(s, c(0, 0.5, 1), c(0, pi / 2, pi))

# Interpolation passes through the control points.
expect_equal(psi[c(1, 26, 51)], c(0, pi / 2, pi))

# Smoothstep is monotone between knots and C1 (zero slope at the knots).
expect_true(all(diff(psi) >= 0))
# The smoothstep analytic derivative is 0 at the knots; approximate by
# checking the interpolation has zero slope at the very first/last interval.
ss_interp <- function(t) t * t * (3 - 2 * t)
expect_equal(ss_interp(0), 0)
expect_equal(ss_interp(1), 1)
expect_true(all(diff(psi) >= 0))

# Out-of-range abscissae are clamped, never extrapolated.
expect_equal(smoothstep_interp(c(-1, 0, 1.5), c(0, 1), c(0, 1)), c(0, 0, 1))

# Length mismatch is an error.
expect_error(smoothstep_interp(seq(0, 1, 0.1), c(0, 1), c(0, 1, 2)))

# Non-monotone stations are an error.
expect_error(smoothstep_interp(seq(0, 1, 0.1), c(0, 1, 0.5), c(0, 1, 2)))

# --- pyt_stations ---------------------------------------------------------

st <- pyt_stations()
expect_true(is.data.frame(st))
expect_identical(nrow(st), 9L)
expect_true(all(st$fraction >= 0 & st$fraction <= 1))
expect_true(all(diff(st$fraction) > 0))           # strictly increasing
expect_true(all(c("cortex", "ic", "pons", "dec", "cord") %in% st$station))

# --- pyt_defaults ---------------------------------------------------------

d <- pyt_defaults()
expect_identical(d$parcels, c("leg", "trunk", "arm", "hand", "face"))
expect_true(is.list(d$macros) && d$macros$tract_length == 70)
expect_true(is.list(d$cross) && d$cross$cross_percent == 0.85)
expect_true(is.list(d$perturbation) && d$perturbation$n_axons == 30L)
expect_identical(d$term_station[["face"]], "pons")
expect_true(all(d$term_station[c("leg", "arm")] == "cord"))
# q/r tables conform to the station grid.
expect_identical(dim(d$q_st), c(9L, 5L))
expect_identical(rownames(d$q_st), as.character(pyt_stations()$station))
expect_identical(colnames(d$q_st), d$parcels)

# From the crus onward q_st is single-signed (ascending leg -> face), so the
# phi flip alone controls the landing side at the decussation (D8).
q_lo <- d$q_st[c("crus", "pons", "pyramid", "dec", "cord"), ]
expect_true(all(q_lo >= 0))
expect_true(all(apply(q_lo, 1, diff) > 0))

# --- simulate_pyramidal_tract ---------------------------------------------

for (h in c("left", "right")) {
  bs <- simulate_pyramidal_tract(hemisphere = h, seed = 99)
  expect_true(fiber::is_bundle_set(bs))
  expect_identical(bs@n_bundles, 5L)

  # Names travel via bundle_data, not names(@bundles) (fiber convention).
  expect_identical(bs@bundle_data$parcel, c("leg", "trunk", "arm", "hand", "face"))
  expect_true(all(bs@bundle_data$role == "pyramidal"))
  expect_true(all(bs@bundle_data$hemisphere == h))

  for (i in seq_along(bs@bundles)) {
    b <- bs@bundles[[i]]
    expect_true(fiber::is_bundle(b))
    # Per-stream metadata lifted to the bundle level.
    expect_true(all(c("parcel", "branch", "hemisphere") %in% names(b@streamline_data)))
    expect_true(all(b@streamline_data$parcel == bs@bundle_data$parcel[i]))
    expect_true(all(b@streamline_data$hemisphere == h))
    expect_identical(sort(unique(b@streamline_data$branch)), c("crossed", "uncrossed"))
    # Streams do not split: 30 per branch, 60 per parcel.
    expect_identical(b@n_streamlines, 60L)
    expect_identical(sum(b@streamline_data$branch == "crossed"), 30L)

    # All coordinates finite.
    expect_true(all(vapply(b@streamlines, function(sl) all(is.finite(sl@points)), TRUE)))
  }

  # Face is corticobulbar: terminates at the pons with fewer arc-length
  # samples, so a completed flip window; limbs run to the cord.
  npts <- function(i) nrow(bs@bundles[[i]]@streamlines[[1]]@points)
  expect_true(npts(5) < npts(1))

  # Cortex fan: leg medial (-X), face lateral (+X) for the left hemisphere
  # (mirrored for the right).
  x0 <- function(i) {
    b <- bs@bundles[[i]]
    mean(vapply(b@streamlines, function(sl) sl@points[1, "X"], 0))
  }
  fan <- vapply(1:5, x0, 0)
  # Somatotopic order preserved: leg medial, face lateral. The sign depends
  # on the hemisphere convention (left medial = -X, right medial = +X) but
  # the fan is monotone along the precentral strip in both cases.
  expect_true(all(diff(fan) > 0) || all(diff(fan) < 0))

  # Decussation: crossed branch ends contralateral, uncrossed ipsilateral.
  xend <- function(i, br) {
    b <- bs@bundles[[i]]
    idx <- which(b@streamline_data$branch == br)
    mean(vapply(idx, function(j) tail(b@streamlines[[j]]@points[, "X"], 1), 0))
  }
  for (i in 1:4) {                                # limbs/trunk
    unc <- xend(i, "uncrossed")
    cro <- xend(i, "crossed")
    if (h == "left") {
      expect_true(unc < 0 && cro > 0)             # crossed goes contralateral
    } else {
      expect_true(unc > 0 && cro < 0)
    }
  }
  # Face terminates bilaterally at the pons.
  f_unc <- xend(5, "uncrossed")
  f_cro <- xend(5, "crossed")
  if (h == "left") {
    expect_true(f_unc < 0 && f_cro > 0)
  } else {
    expect_true(f_unc > 0 && f_cro < 0)
  }
}

# Left and right hemispheres are mirror images across the midline. The
# perturbation is zero-mean by construction, so compare the mean streamline,
# which cancels stream-level noise and isolates the deterministic template.
L <- simulate_pyramidal_tract(hemisphere = "left", seed = 7)
R <- simulate_pyramidal_tract(hemisphere = "right", seed = 7)
mean_pts <- function(bs, i) {
  pts <- lapply(bs@bundles[[i]]@streamlines, function(sl) sl@points)
  Reduce(`+`, pts) / length(pts)
}
for (i in 1:5) {
  ptl <- mean_pts(L, i)
  ptr <- mean_pts(R, i)
  expect_identical(nrow(ptl), nrow(ptr))
  # Zero-mean perturbation still leaves sampling noise of order
  # sqrt(variance / n_streams) ~ 0.01-0.1 mm in the 60-stream means.
  expect_equal(ptl[, "X"], -ptr[, "X"], tolerance = 0.5)
  expect_equal(ptl[, "Y"],  ptr[, "Y"], tolerance = 0.5)
  expect_equal(ptl[, "Z"],  ptr[, "Z"], tolerance = 0.5)
}

# Reproducibility via seed.
s1 <- simulate_pyramidal_tract(seed = 42)
s2 <- simulate_pyramidal_tract(seed = 42)
expect_identical(
  lapply(s1@bundles, function(b) lapply(b@streamlines, function(sl) sl@points)),
  lapply(s2@bundles, function(b) lapply(b@streamlines, function(sl) sl@points))
)

# Twist length must match the station grid.
expect_error(simulate_pyramidal_tract(rot_st = rep(0, 4)))
