## Polyester (<=1.38.0) add_platform_error hard-stops for reads >101 bp, which
## blocks custom error models for 150-bp reads. This patch removes only that
## guard: deparse -> drop the 3-line guard -> eval -> assignInNamespace.
## Usage: library(polyester); then source() this file.
patch_polyester_150 <- function() {
  src <- deparse(polyester:::add_platform_error)
  gs <- grep("> 101", src, fixed = TRUE)
  stopifnot(length(gs) == 1)
  stopifnot(grepl("102bp", src[gs + 1], fixed = TRUE))
  stopifnot(src[gs + 2] == "    }")
  fixed_src <- src[-(gs:(gs + 2))]
  fun <- eval(parse(text = paste(fixed_src, collapse = "\n")))
  environment(fun) <- asNamespace("polyester")
  assignInNamespace("add_platform_error", fun, ns = "polyester")
  invisible(TRUE)
}
patch_polyester_150()
