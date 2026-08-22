manifest_path <- file.path("starlight", ".starlightdown", "site.json")
dist_dir <- file.path("starlight", "dist")

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

if (!file.exists(manifest_path) || !dir.exists(dist_dir)) {
  stop("Build the Starlight site before running this check.", call. = FALSE)
}

site <- jsonlite::read_json(manifest_path, simplifyVector = TRUE)
routes <- site$routes
base <- sub("/$", "", site$site$base %||% "")

failures <- character()

for (i in seq_len(nrow(routes))) {
  route <- routes$route[[i]]
  id <- routes$id[[i]]
  kind <- routes$kind[[i]]
  relative <- sub("^/", "", sub("/$", "", route))
  page <- if (nzchar(relative)) {
    file.path(dist_dir, relative, "index.html")
  } else {
    file.path(dist_dir, "index.html")
  }

  if (!file.exists(page)) {
    failures <- c(failures, paste0(route, ": missing ", page))
    next
  }

  doc <- xml2::read_html(page)
  title <- xml2::xml_text(xml2::xml_find_first(doc, "//title"))
  h1 <- trimws(xml2::xml_text(xml2::xml_find_first(doc, "//h1")))
  expected_h1 <- switch(
    kind,
    home = site$package$name,
    article = routes$title[[i]],
    `reference-index` = routes$title[[i]],
    reference = basename(id),
    routes$title[[i]]
  )

  if (is.na(h1) || !identical(h1, expected_h1)) {
    failures <- c(
      failures,
      paste0(route, ": expected H1 '", expected_h1, "', found '", h1, "'")
    )
  }
  if (startsWith(title, "Redirecting to ") ||
      length(xml2::xml_find_all(doc, "//meta[translate(@http-equiv, 'REFSH', 'refsh')='refresh']"))) {
    failures <- c(failures, paste0(route, ": resolved to a redirect document"))
  }

  canonical <- xml2::xml_attr(
    xml2::xml_find_first(doc, "//link[@rel='canonical']"),
    "href"
  )
  canonical_path <- sub("^https?://[^/]+", "", canonical)
  expected_path <- paste0(base, route)
  if (is.na(canonical) || !identical(canonical_path, expected_path)) {
    failures <- c(
      failures,
      paste0(route, ": canonical path '", canonical_path,
             "' does not equal '", expected_path, "'")
    )
  }

  images <- xml2::xml_find_all(doc, "//main//img")
  if (length(images)) {
    alt <- trimws(xml2::xml_attr(images, "alt"))
    bad <- is.na(alt) | !nzchar(alt)
    if (any(bad)) {
      failures <- c(
        failures,
        paste0(route, ": ", sum(bad), " content image(s) have empty alt text")
      )
    }
  }
}

if (length(failures)) {
  stop(
    paste(c("Starlight route identity check failed:", paste0("- ", failures)),
          collapse = "\n"),
    call. = FALSE
  )
}

message(
  "Verified ", nrow(routes),
  " advertised routes by file, H1 identity, canonical URL, redirect state, and image alt text."
)
