test_that("flag reports exceedances rather than abnormalities", {
  d <- norm_dist("gaussian", location = 0, scale = 1)
  sc <- as_scores(d, c(0, 3, -2.5))
  sc$.in_sample <- FALSE
  fl <- flag(sc, threshold = 2)
  expect_equal(sum(fl$exceedance), 2)
  expect_false(any(grepl("abnormal", names(fl), ignore.case = TRUE)))
})

test_that("as_wide and filter_scores work", {
  d <- norm_dist("gaussian", location = c(0, 0), scale = 1)
  sc <- as_scores(d, c(1, -1))
  sc$.row <- 1:2
  sc$.id <- c("A", "B")
  sc$.outcome <- "y"
  sc$.in_sample <- FALSE
  w <- as_wide(sc, value = "z")
  expect_true("y" %in% names(w))
  f <- filter_scores(sc, abs(z) > 0.5)
  expect_equal(nrow(f), 2)
})
