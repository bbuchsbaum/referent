simple_spec <- function(family = "gaussian", scale = FALSE) {
  loc <- ~ s(age, k = 5) + sex
  sc <- if (isTRUE(scale)) ~ s(age, k = 4) else ~ 1
  if (identical(family, "shash")) {
    norm_spec(
      family = norm_shash(),
      location = loc,
      scale = sc,
      skew = ~1,
      tail = ~1
    )
  } else {
    norm_spec(family = norm_gaussian(), location = loc, scale = sc)
  }
}

expect_near <- function(x, y, tol = 1e-6) {
  expect_true(all(abs(x - y) < tol, na.rm = TRUE))
}
