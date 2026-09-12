#' Run an AMOVA with a significance test
#'
#' Wraps [poppr::poppr.amova()] plus [ade4::randtest()] so the analysis and its
#' permutation test are produced together, with the seed recorded, rather than
#' as two loosely related objects.
#'
#' @param gen A `genind` or `genclone` object with strata assigned.
#' @param hier A hierarchical formula, e.g. `~Pop` or `~Pop/Subpop`. The terms
#'   must exist in `strata(gen)`.
#' @param nrepet Number of permutations for the significance test.
#'   Default 999.
#' @param seed Random seed, so the test is reproducible. Default 1999.
#' @param ... Passed through to [poppr::poppr.amova()].
#'
#' @return A list of class `popgen_amova` with elements `amova`, `test`,
#'   `hierarchy`, `nrepet` and `seed`.
#' @export
#'
#' @examples
#' \dontrun{
#' gen <- read_genalex_csv("data/nj.csv")
#' res <- run_amova(gen, ~Pop)
#' res$amova
#' plot(res$test)
#' }
run_amova <- function(gen, hier = ~Pop, nrepet = 999, seed = 1999, ...) {
  stopifnot(inherits(gen, "genind"))
  if (!inherits(hier, "formula")) {
    stop("`hier` must be a formula such as ~Pop or ~Pop/Subpop.", call. = FALSE)
  }

  strata_df <- poppr::strata(gen)
  if (is.null(strata_df)) {
    stop("`gen` has no strata. Assign them with poppr::strata(gen) <- ... ",
         "before running AMOVA.", call. = FALSE)
  }

  wanted  <- all.vars(hier)
  missing <- setdiff(wanted, names(strata_df))
  if (length(missing)) {
    stop("Strata not found in `gen`: ", paste(missing, collapse = ", "),
         ".\n  Available: ", paste(names(strata_df), collapse = ", "),
         call. = FALSE)
  }

  amova_res <- poppr::poppr.amova(gen, hier, ...)

  set.seed(seed)
  test_res <- ade4::randtest(amova_res, nrepet = nrepet)

  structure(
    list(
      amova     = amova_res,
      test      = test_res,
      hierarchy = deparse(hier),
      nrepet    = nrepet,
      seed      = seed
    ),
    class = "popgen_amova"
  )
}

#' @export
print.popgen_amova <- function(x, ...) {
  cat("AMOVA:", x$hierarchy, "\n")
  cat("Permutations:", x$nrepet, " Seed:", x$seed, "\n\n")
  print(x$amova)
  cat("\nSignificance test:\n")
  print(x$test)
  invisible(x)
}
