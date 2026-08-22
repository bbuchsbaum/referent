test_that("flag reports exceedances rather than abnormalities", {
  d <- distributional::dist_normal(0, 1)
  sc <- as_scores(d, c(0, 3, -2.5))
  sc$.in_sample <- FALSE
  fl <- flag(sc, threshold = 2)
  expect_equal(sum(fl$exceedance), 2)
  expect_false(any(grepl("abnormal", names(fl), ignore.case = TRUE)))
})

test_that("scores_matrix pivots a long table and keeps ids", {
  d <- distributional::dist_normal(c(0, 0), 1)
  sc <- as_scores(d, c(1, -1))
  sc$.row <- 1:2
  sc$.id <- c("A", "B")
  sc$.outcome <- "y"
  w <- scores_matrix(sc, value = "z")
  expect_equal(colnames(w$matrix), "y")
  expect_equal(w$.id, c("A", "B"))
  expect_equal(as.numeric(w$matrix[, "y"]), sc$z)
  expect_error(scores_matrix(sc, value = "nope"), "Unknown score column")
})
