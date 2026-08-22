test_that("norm_flag reports exceedances and FDR rather than abnormalities", {
  d <- distributional::dist_normal(0, 1)
  y <- c(0, 3, -2.5, 0.5, 1.9)
  sc <- as_scores(d, y)
  sc$.in_sample <- FALSE
  fl <- norm_flag(sc, threshold = 2)
  expect_s3_class(fl, "norm_flags")
  expect_equal(fl$exceedance, abs(y) > 2)
  expect_equal(attr(fl, "observed_exceedances"), 2L)
  expect_equal(attr(fl, "expected_exceedances"), 2 * stats::pnorm(-2) * 5)
  expect_equal(fl$fdr, stats::p.adjust(sc$tail_prob, method = "fdr"))
  expect_false(any(grepl("abnormal", names(fl), ignore.case = TRUE)))
  body <- utils::capture.output(msg <- cli::cli_fmt(print(fl)))
  expect_match(paste(msg, collapse = "\n"), "2 exceedances")
  expect_true(any(grepl("exceedance", body)))
  # in-sample scores warn
  sc_in <- sc
  sc_in$.in_sample <- TRUE
  expect_warning(norm_flag(sc_in), "in-sample")
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
