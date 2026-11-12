# Tests for CSV export.

cl <- generate_centerline(angle = 0.5, n_points = 20)
b <- create_axons(cl, n_axons = 3, radius = 0.1, seed = 1)

f <- tempfile(fileext = ".csv")
expect_identical(generate_csv(b, f), f)
d <- read.csv(f)
# 3 streamlines x 20 points each.
expect_identical(nrow(d), 60L)
expect_true(all(c("X", "Y", "Z") %in% colnames(d)))

f2 <- tempfile(fileext = ".csv")
generate_csv(b, f2, attribute = "curvature")
d2 <- read.csv(f2)
expect_true("curvature" %in% colnames(d2))
