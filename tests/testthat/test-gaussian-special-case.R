test_that("package Z-scores agree with (y-mu)/sigma", {
  y <- c(-1, 0, 0.5, 2)
  mu <- 0.25
  sigma <- 1.5
  d <- norm_dist("gaussian", location = mu, scale = sigma)
  sc <- as_scores(d, y)
  expect_equal(sc$z, (y - mu) / sigma, tolerance = 1e-10)
  expect_equal(sc$centile, stats::pnorm(y, mu, sigma), tolerance = 1e-10)
})

test_that("as_scores has no abnormality column", {
  d <- norm_dist("gaussian", location = 0, scale = 1)
  sc <- as_scores(d, 0)
  expect_false(any(grepl("abnormal", names(sc), ignore.case = TRUE)))
})
