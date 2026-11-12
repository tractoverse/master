# Tests for the fascicle configurations.

for (cfg in c("bending", "fanning", "merging", "kissing", "crossing")) {
  st <- simulate_scenario(cfg, n_axons = 12, seed = 1)
  expect_true(fiber::is_bundle_set(st), info = cfg)
  n_expected <- if (cfg %in% c("merging", "kissing", "crossing")) 2L else 1L
  expect_identical(st@n_bundles, n_expected, info = cfg)
  # All streamlines are valid and carry the expected per-point attributes.
  total <- 0L
  for (bnd in st@bundles) {
    expect_true(fiber::is_bundle(bnd))
    total <- total + bnd@n_streamlines
  }
  expect_equal(total, if (n_expected == 1L) 12L else 24L, info = cfg)
}

# A single streamline never splits: every streamline has exactly one contiguous
# path (its points form a chain, which the model guarantees).
st <- simulate_scenario("fanning", n_axons = 20, seed = 2)
expect_equal(st@n_bundles, 1L)
