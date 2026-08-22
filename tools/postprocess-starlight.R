reference_dir <- file.path("starlight", "src", "content", "docs", "reference")

if (!dir.exists(reference_dir)) {
  stop("Run starlightdown::build_site() before this postprocessor.", call. = FALSE)
}

# starlightdown currently emits empty labels for figures evaluated from Rd
# examples. These two reference plots are meaningful, so supply concrete
# alternatives after generation and fail if any empty label remains.
autoplot_page <- file.path(reference_dir, "autoplot.ref_fit.md")
if (!file.exists(autoplot_page)) {
  stop("Missing generated page: ", autoplot_page, call. = FALSE)
}

lines <- readLines(autoplot_page, warn = FALSE, encoding = "UTF-8")
replacements <- c(
  "![](figures/autoplot.ref_fit-1.png)" =
    "![Fitted median and centile bands across age, faceted by sex, with observations overlaid.](figures/autoplot.ref_fit-1.png)",
  "![](figures/autoplot.ref_fit-2.png)" =
    "![Heatmap of individual deviation scores across observations and outcomes.](figures/autoplot.ref_fit-2.png)"
)

for (source in names(replacements)) {
  lines[lines == source] <- unname(replacements[[source]])
}

if (any(grepl("^!\\[\\]", lines))) {
  stop("An example image still has empty alternative text in ", autoplot_page,
       call. = FALSE)
}

writeLines(lines, autoplot_page, useBytes = TRUE)
message("Verified generated reference-plot alt text.")
