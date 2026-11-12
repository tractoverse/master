# Tests for axon generation around a centreline.

set.seed(99)
cl <- generate_centerline(angle = 0.6, arm_length = 1, interaction_length = 2, n_points = 150)

b <- create_axons(cl, n_axons = 20, n_basis = 8, variance = 0.02,
                  autocorrelation = 2, seed = 42)

expect_true(fiber::is_bundle(b))
expect_identical(b@n_streamlines, 20L)

# Every streamline has the per-point microstructure attributes.
sl <- b[[1]]
expect_true(all(c("curvature", "torsion", "dir_x", "dir_y", "dir_z") %in% sl@point_attributes))
expect_identical(length(sl@point_data$curvature), 150L)
expect_identical(length(sl@point_data$dir_x), 150L)

# Zero-mean property: the sample mean of the bundle stays on the centreline
# (within sampling noise for a modest number of axons).
pts <- do.call(rbind, lapply(b@streamlines, function(s) s@points))
mn <- colMeans(pts)
expect_true(max(abs(mn - cl[75, ])) < 0.5)

# Reproducibility via the seed.
b2 <- create_axons(cl, n_axons = 20, n_basis = 8, variance = 0.02,
                   autocorrelation = 2, seed = 42)
pts2 <- do.call(rbind, lapply(b2@streamlines, function(s) s@points))
expect_equal(pts, pts2)

# Radius handled as metadata: a scalar is lifted to the bundle-level
# `streamline_data` (one entry per streamline); a per-point vector stays in
# each streamline's `point_data`. fiber lifts common streamline_data to the
# bundle, so the scalar radius is read back from `b_r@streamline_data`.
b_r <- create_axons(cl, n_axons = 3, radius = 0.1, seed = 1)
expect_true("radius" %in% names(b_r@streamline_data))
expect_true(all(b_r@streamline_data$radius == 0.1))
# The accessor re-pushes the bundle-level radius back into each streamline.
expect_identical(b_r[[1]]@streamline_data$radius, 0.1)
bv <- create_axons(cl, n_axons = 3, radius = seq_len(150) * 0.001, seed = 1)
expect_identical(bv[[1]]@point_data$radius[1], 0.001)

# Exponential autocorrelation of the coefficient covariance: entries decay
# exponentially with lag and are symmetric/PD.
Sig <- master:::.bspline_covariance(10, variance = 1, autocorrelation = 2)
expect_true(isSymmetric(Sig))
expect_true(all(eigen(Sig, symmetric = TRUE)$values > 0))
expect_equal(Sig[1, 2], exp(-1 / 2))
