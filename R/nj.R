#' Bootstrapped neighbour-joining tree
#'
#' Builds an NJ tree and its bootstrap support values in one call, and drops
#' support below a cutoff so the plotted tree stays readable.
#'
#' @param gen A `genind` object.
#' @param method Distance measure, passed to [genetic_distance()].
#'   Default `"nei"`.
#' @param bootstrap Number of bootstrap replicates. Default 1000. Set to 0 to
#'   skip bootstrapping, which is much faster while you are still iterating on
#'   a figure.
#' @param support_cutoff Support values below this percentage are set to `NA`
#'   and therefore not drawn. Default 70.
#' @param replen Repeat lengths per locus, required for `method = "bruvo"`.
#'
#' @return A list of class `popgen_nj` with elements `tree` (a `phylo`),
#'   `support` (percentages, `NA` below the cutoff), `raw_support`, `dist`,
#'   `method`, `bootstrap` and `support_cutoff`.
#' @export
#'
#' @examples
#' \dontrun{
#' gen  <- read_genalex_csv("data/nj.csv")
#' # start with bootstrap = 0 while tuning the figure, then raise it
#' fit  <- nj_bootstrap(gen, bootstrap = 100)
#' plot_nj(fit, groups = adegenet::pop(gen))
#' }
nj_bootstrap <- function(gen,
                         method = "nei",
                         bootstrap = 1000,
                         support_cutoff = 70,
                         replen = NULL) {
  stopifnot(inherits(gen, "genind"))

  d   <- genetic_distance(gen, method = method, replen = replen)
  tre <- ape::nj(d)

  support <- NULL
  if (bootstrap > 0) {
    tab <- adegenet::tab(gen, NA.method = "mean")
    raw <- ape::boot.phylo(
      tre,
      tab,
      FUN = function(e) ape::nj(stats::dist(e)),
      B = bootstrap,
      quiet = TRUE
    )
    support <- 100 * raw / bootstrap
  }

  trimmed <- support
  if (!is.null(trimmed)) {
    trimmed <- round(trimmed)
    trimmed[trimmed < support_cutoff] <- NA_real_
  }

  structure(
    list(
      tree           = tre,
      support        = trimmed,
      raw_support    = support,
      dist           = d,
      method         = method,
      bootstrap      = bootstrap,
      support_cutoff = support_cutoff
    ),
    class = "popgen_nj"
  )
}

#' Plot a bootstrapped neighbour-joining tree
#'
#' @param fit A `popgen_nj` object from [nj_bootstrap()].
#' @param groups Factor giving each tip's group, used for tip colours.
#'   Usually `adegenet::pop(gen)`. If `NULL`, tips are drawn in black.
#' @param palette Passed to [cluster_colours()].
#' @param type Tree layout: `"unrooted"`, `"phylogram"`, `"fan"` or `"radial"`.
#'   Default `"unrooted"`.
#' @param cex Tip label size. Default 0.5.
#' @param show_legend Draw a legend for the groups? Default `TRUE`.
#' @param ... Passed to [ape::plot.phylo()].
#'
#' @return `fit`, invisibly. Called for the plot it draws.
#' @export
plot_nj <- function(fit,
                    groups = NULL,
                    palette = NULL,
                    type = "unrooted",
                    cex = 0.5,
                    show_legend = TRUE,
                    ...) {
  stopifnot(inherits(fit, "popgen_nj"))

  if (is.null(groups)) {
    cols   <- rep("black", length(fit$tree$tip.label))
    legend_map <- NULL
  } else {
    if (length(groups) != length(fit$tree$tip.label)) {
      stop("`groups` has ", length(groups), " entries but the tree has ",
           length(fit$tree$tip.label), " tips.", call. = FALSE)
    }
    cols       <- cluster_colours(groups, palette = palette)
    legend_map <- attr(cols, "legend")
  }

  ape::plot.phylo(
    fit$tree,
    type         = type,
    cex          = cex,
    tip.color    = cols,
    edge.width   = 1,
    label.offset = 0.05,
    show.tip.label = TRUE,
    ...
  )

  if (!is.null(fit$support)) {
    ape::nodelabels(
      fit$support,
      frame = "c", cex = 0.4, bg = "black", col = "white"
    )
  }

  ape::tiplabels(col = cols, frame = "circle", pch = 16, cex = 0.7)

  if (show_legend && !is.null(legend_map)) {
    graphics::legend(
      "topleft",
      legend = names(legend_map),
      fill   = unname(legend_map),
      cex    = 0.8,
      bty    = "n",
      border = unname(legend_map)
    )
  }

  ape::add.scale.bar("bottomleft", cex = 0.7, font = 2, col = "black")

  invisible(fit)
}

#' @export
print.popgen_nj <- function(x, ...) {
  cat("Neighbour-joining tree\n")
  cat("  tips:      ", length(x$tree$tip.label), "\n")
  cat("  distance:  ", x$method, "\n")
  cat("  bootstrap: ", x$bootstrap, " replicates\n")
  if (!is.null(x$raw_support)) {
    kept <- sum(!is.na(x$support))
    cat("  nodes with support >= ", x$support_cutoff, "%: ",
        kept, " of ", length(x$support), "\n", sep = "")
  }
  invisible(x)
}
