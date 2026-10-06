#!/usr/bin/env Rscript
# PESTFLY: isolated Hudson FST and mutual information implementation checks.
# Run from the repository or extracted S6 root. This uses base R only.
# Parse the original scripts and evaluate only the named function definitions.
# No option parsing, sequence processing or genomic analysis is executed.

load_functions <- function(path, wanted, into) {
  code <- parse(file = path)
  loaded <- character()
  for (expr in code) {
    if (!is.call(expr) || length(expr) != 3L ||
        !identical(expr[[1]], as.name("<-")) || !is.symbol(expr[[2]])) next
    name <- as.character(expr[[2]])
    if (!(name %in% wanted)) next
    value <- expr[[3]]
    if (!is.call(value) || !identical(value[[1]], as.name("function"))) next
    eval(expr, envir = into)
    loaded <- c(loaded, name)
  }
  if (!setequal(loaded, wanted)) {
    stop("Required implementation functions not found: ",
         paste(setdiff(wanted, loaded), collapse = ", "))
  }
}

fst <- new.env(parent = baseenv())
load_functions("steps/02_snp_panel_discovery_and_ranking/run.R",
               c("hudson_fst", "validate_hudson_fst"), fst)
fst$validate_hudson_fst()
stopifnot(is.na(fst$hudson_fst(0.5, 0.5, 1L, 4L)),
          is.na(fst$hudson_fst(NA_real_, 0.5, 4L, 4L)),
          is.na(fst$hudson_fst(Inf, 0.5, 4L, 4L)))
cat("Hudson FST checks passed, including undefined and insufficient call cases.\n")

mi <- new.env(parent = baseenv())
load_functions("steps/03_mutual_information_snp_diagnostics/run.R",
               c("entropy_bits", "mutual_information_bits", "run_mi_self_test"), mi)
mi$run_mi_self_test()
stopifnot(is.na(mi$mutual_information_bits(c(NA, NA), c("A", "B"))))
cat("Implementation checks completed.\n")
