#' Principal coordinates analysis
#'
#' Runs PCoA on a genetic distance matrix and, crucially, returns the variance
#' explained by each axis alongside the coordinates -- so axis labels can be
#' generated from the data rather than typed in by hand.
#'
#' Nei's distance is not guaranteed to be Euclidean, which makes
#' [ade4::dudi.pco()] produce negative eigenvalues. `correction` applies a
#' standard fix; the default `"cailliez"` is the usual choice.
#'
#' @param gen A `genind` object.
#' @param method Distance measure, passed to [genetic_distance()].
#'   Default `"nei"`.
#' @param nf Number of axes to keep. Default 3.
#' @param correction `"cailliez"`, `"lingoes"` or `"none"`.
#'   Default `"cailliez"`.
#' @param replen Repeat lengths per locus, required for `method = "bruvo"`.
#'
#' @return A list of class `popgen_pcoa` with elements `pco` (the `dudi`
#'   object), `coords` (a data frame of axis scores with individual names),
#'   `variance` (percentage explained per axis), `dist`, `method` and
#'   `correction`.
#' @export
#'
#' @examples
#' \dontrun{
#' gen <- read_genalex_csv("data/pcoa.csv")
#' pc  <- run_pcoa(gen)
#' pc$variance[1:2]
#' plot_pcoa(pc, groups = adegenet::pop(gen))
#' }
run_pcoa <- function(gen,
                     method = "nei",
                     nf = 3,
                     correction = c("cailliez", "lingoes", "none"),
                     replen = NULL) {
  correction <- match.arg(correction)
  stopifnot(inherits(gen, "genind"))

  d <- genetic_distance(gen, method = method, replen = replen)

  d_corrected <- switch(
    correction,
    cailliez = ade4::cailliez(d, print = FALSE),
    lingoes  = ade4::lingoes(d, print = FALSE),
    none     = d
  )

  pco <- ade4::dudi.pco(d_corrected, scannf = FALSE, nf = nf)

  positive <- pco$eig[pco$eig > 0]
  variance <- round(100 * positive / sum(positive), 2)
  names(variance) <- paste0("PCo", seq_along(variance))

  coords <- as.data.frame(pco$li)
  names(coords) <- paste0("PCo", seq_len(ncol(coords)))
  coords$individual <- adegenet::indNames(gen)

  structure(
    list(
      pco        = pco,
      coords     = coords,
      variance   = variance,
      dist       = d,
      method     = method,
      correction = correction
    ),
    class = "popgen_pcoa"
  )
}

#' Plot a PCoA ordination
#'
#' Draws two PCoA axes with ggplot2, labelling each axis with the variance it
#' actually explains. Built as a plain ggplot, so it can be extended with `+`
#' in the normal way rather than patched through `ggplot_build()` internals.
#'
#' @param pc A `popgen_pcoa` object from [run_pcoa()].
#' @param groups Factor giving each individual's group, used for point colour.
#'   If `NULL`, all points are drawn in one colour.
#' @param axes Length-2 integer vector: which axes to plot. Default `c(1, 2)`.
#' @param palette Passed to [cluster_colours()].
#' @param label_points Draw individual names next to the points?
#'   Default `FALSE` -- useful for checking outliers, noisy for publication.
#' @param point_size Point size. Default 2.
#' @param ellipse Draw a normal-probability ellipse per group? Default `TRUE`.
#'
#' @return A `ggplot` object.
#' @export
plot_pcoa <- function(pc,
                      groups = NULL,
                      axes = c(1, 2),
                      palette = NULL,
                      label_points = FALSE,
                      point_size = 2,
                      ellipse = TRUE) {
  stopifnot(inherits(pc, "popgen_pcoa"))
  if (length(axes) != 2) {
    stop("`axes` must be a vector of two axis numbers, e.g. c(1, 2).",
         call. = FALSE)
  }
  if (max(axes) > ncol(pc$coords) - 1) {
    stop("Only ", ncol(pc$coords) - 1, " axes were kept. ",
         "Re-run run_pcoa() with a larger `nf`.", call. = FALSE)
  }

  df <- pc$coords
  x_col <- paste0("PCo", axes[1])
  y_col <- paste0("PCo", axes[2])

  if (is.null(groups)) {
    df$group <- factor("all")
    lookup   <- c(all = "#0072B2")
  } else {
    if (length(groups) != nrow(df)) {
      stop("`groups` has ", length(groups), " entries but there are ",
           nrow(df), " individuals.", call. = FALSE)
    }
    df$group <- as.factor(groups)
    lookup   <- attr(cluster_colours(df$group, palette = palette), "legend")
  }

  x_lab <- sprintf("PCo%d (%.2f%%)", axes[1], pc$variance[axes[1]])
  y_lab <- sprintf("PCo%d (%.2f%%)", axes[2], pc$variance[axes[2]])

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], colour = .data$group)
  )

  if (ellipse && nlevels(df$group) > 1) {
    p <- p + ggplot2::stat_ellipse(type = "norm", linewidth = 0.4,
                                   show.legend = FALSE)
  }

  p <- p +
    ggplot2::geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::scale_colour_manual(values = lookup, name = NULL) +
    ggplot2::labs(x = x_lab, y = y_lab) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey92",
                                               linewidth = 0.2),
      legend.background = ggplot2::element_rect(fill = "transparent")
    )

  if (label_points) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = .data$individual),
      size = 1.9, vjust = -0.9, show.legend = FALSE
    )
  }

  p
}

#' @export
print.popgen_pcoa <- function(x, ...) {
  cat("Principal coordinates analysis\n")
  cat("  distance:   ", x$method, "\n")
  cat("  correction: ", x$correction, "\n")
  cat("  individuals:", nrow(x$coords), "\n")
  cat("  variance explained (first axes):\n")
  print(utils::head(x$variance, 5))
  invisible(x)
}
