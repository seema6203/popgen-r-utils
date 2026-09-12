#' Per-population diversity summary
#'
#' Collects the statistics a diversity table usually reports -- number of
#' individuals, number of multilocus genotypes, observed and expected
#' heterozygosity, and F_IS -- into a single tidy data frame, instead of
#' reading them off three different printed objects.
#'
#' @param gen A `genind` object with populations assigned.
#'
#' @return A data frame with one row per population and an `Overall` row.
#'   Columns: `population`, `n`, `mlg`, `Ho`, `Hs`, `Fis`.
#' @export
#'
#' @examples
#' \dontrun{
#' gen <- read_genalex_csv("data/nj.csv")
#' diversity_summary(gen)
#' }
diversity_summary <- function(gen) {
  stopifnot(inherits(gen, "genind"))
  if (is.null(adegenet::pop(gen))) {
    stop("`gen` has no population factor; assign one with adegenet::pop().",
         call. = FALSE)
  }

  hf    <- hierfstat::genind2hierfstat(gen)
  stats <- hierfstat::basic.stats(hf)

  # basic.stats returns loci x population matrices; average over loci.
  ho <- colMeans(stats$Ho,  na.rm = TRUE)
  hs <- colMeans(stats$Hs,  na.rm = TRUE)
  fis <- colMeans(stats$Fis, na.rm = TRUE)

  pops <- names(ho)
  n    <- as.integer(table(adegenet::pop(gen))[pops])

  mlg <- vapply(
    pops,
    function(p) {
      sub <- gen[adegenet::pop(gen) == p, ]
      nrow(unique(adegenet::tab(sub)))
    },
    integer(1)
  )

  per_pop <- data.frame(
    population = pops,
    n          = n,
    mlg        = as.integer(mlg),
    Ho         = round(ho,  4),
    Hs         = round(hs,  4),
    Fis        = round(fis, 4),
    row.names  = NULL,
    stringsAsFactors = FALSE
  )

  overall <- data.frame(
    population = "Overall",
    n          = adegenet::nInd(gen),
    mlg        = nrow(unique(adegenet::tab(gen))),
    Ho         = round(unname(stats$overall["Ho"]),  4),
    Hs         = round(unname(stats$overall["Hs"]),  4),
    Fis        = round(unname(stats$overall["Fis"]), 4),
    row.names  = NULL,
    stringsAsFactors = FALSE
  )

  rbind(per_pop, overall)
}

#' Pairwise F_ST between populations
#'
#' @param gen A `genind` object with populations assigned.
#' @param method `"WC"` for Weir & Cockerham (1984) or `"nei"` for Nei (1987).
#'   Default `"WC"`.
#'
#' @return A symmetric matrix of pairwise F_ST values, with population names as
#'   dimnames. The diagonal is zero.
#' @export
pairwise_fst <- function(gen, method = c("WC", "nei")) {
  method <- match.arg(method)
  stopifnot(inherits(gen, "genind"))
  if (is.null(adegenet::pop(gen))) {
    stop("`gen` has no population factor; assign one with adegenet::pop().",
         call. = FALSE)
  }
  if (adegenet::nPop(gen) < 2) {
    stop("Need at least two populations to compute pairwise F_ST.",
         call. = FALSE)
  }

  hf <- hierfstat::genind2hierfstat(gen)

  m <- switch(
    method,
    WC  = hierfstat::pairwise.WCfst(hf),
    nei = hierfstat::pairwise.neifst(hf)
  )

  m <- as.matrix(m)
  m[is.na(m)] <- 0
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  diag(m) <- 0
  m
}

#' Genetic distance between individuals
#'
#' One entry point for the distance measures used across these analyses, so a
#' script can switch measure without switching package.
#'
#' @param gen A `genind` object.
#' @param method `"nei"` (Nei's standard distance), `"euclidean"` (on the
#'   allele-frequency table) or `"bruvo"` (stepwise-mutation distance for
#'   microsatellites). Default `"nei"`.
#' @param replen Repeat lengths per locus. Required for `method = "bruvo"`.
#'
#' @return A `dist` object over individuals.
#' @export
genetic_distance <- function(gen,
                             method = c("nei", "euclidean", "bruvo"),
                             replen = NULL) {
  method <- match.arg(method)
  stopifnot(inherits(gen, "genind"))

  d <- switch(
    method,
    nei       = poppr::nei.dist(gen),
    euclidean = stats::dist(adegenet::tab(gen, NA.method = "mean")),
    bruvo     = {
      if (is.null(replen)) {
        stop("method = \"bruvo\" needs `replen`, the repeat length of each ",
             "locus (e.g. rep(2, nLoc(gen))).", call. = FALSE)
      }
      poppr::bruvo.dist(gen, replen = replen)
    }
  )

  stats::as.dist(d)
}
