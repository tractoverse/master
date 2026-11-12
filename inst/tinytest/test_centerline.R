# Tests for centerline generation and arc-length reparametrization.

cl <- generate_centerline(
  angle = 0.6, arm_length = 1, interaction_length = 2, n_points = 200
)

expect_true(is.matrix(cl))
expect_identical(dim(cl), c(200L, 3L))
expect_identical(colnames(cl), c("X", "Y", "Z"))
expect_true(is.numeric(cl))

# Points should be ordered in X (monotonic one-way sweep).
expect_true(all(diff(cl[, "X"]) > 0))

# The centreline should be flat (y ~ 0) at the interaction and bend away.
mid <- cl[100, "Y"]
expect_true(abs(mid) < 0.5)

# Arc-length uniformity: successive inter-point distances are almost equal
# because reparametrize() lands points on a uniform curvilinear grid.
d <- sqrt(rowSums(diff(cl)^2))
expect_true(diff(range(d)) < 0.05 * max(d))

# Reproducible (deterministic geometry).
cl2 <- generate_centerline(angle = 0.6, arm_length = 1, interaction_length = 2, n_points = 200)
expect_equal(cl, cl2)

# reparametrize_3D_curve wrapper.
m <- cbind(X = 0:10, Y = (0:10)^2 / 10, Z = 0)
m2 <- reparametrize_3D_curve(m, n_points = 20)
expect_identical(dim(m2), c(20L, 3L))
expect_identical(colnames(m2), c("X", "Y", "Z"))
