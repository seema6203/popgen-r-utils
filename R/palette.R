#' Colours for a grouping factor
#'
#' Returns one colour per individual, looked up from the group each individual
#' belongs to.
#'
#' This exists to replace the positional pattern
#' `c(rep("purple", 2), rep("red", 39), ...)`, which silently produces the wrong
#' figure the moment a sample is added, removed or reordered. Deriving colours
#' from the grouping factor instead keeps the figure correct by construction.
#'
#' @param groups A factor (or something coercible to one) giving each
#'   individual's group -- typically `adegenet::pop(gen)`.
#' @param palette Either a character vector of colours, recycled over the levels
#'   of `groups`, or a named vector mapping level names to colours. Defaults to
#'   a colour-blind-safe eight-colour palette.
#'
#' @return A character vector of colours, the same length as `groups`, with a
#'   `"legend"` attribute holding the level-to-colour mapping.
#' @export
#'
#' @examples
#' groups <- factor(c("A", "A", "B", "C", "B"))
#' cols <- cluster_colours(groups)
#' attr(cols, "legend")
cluster_colours <- function(groups, palette = NULL) {
  groups <- as.factor(groups)
  levs   <- levels(groups)

  if (is.null(palette)) {
    # Okabe-Ito: distinguishable under the common forms of colour blindness.
    palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                 "#0072B2", "#D55E00", "#CC79A7", "#000000")
  }

  if (!is.null(names(palette))) {
    missing <- setdiff(levs, names(palette))
    if (length(missing)) {
      stop("No colour given for group(s): ", paste(missing, collapse = ", "),
           call. = FALSE)
    }
    lookup <- palette[levs]
  } else {
    if (length(palette) < length(levs)) {
      warning("Palette has ", length(palette), " colours for ", length(levs),
              " groups; recycling.", call. = FALSE)
    }
    lookup <- stats::setNames(rep_len(palette, length(levs)), levs)
  }

  out <- unname(lookup[as.character(groups)])
  attr(out, "legend") <- lookup
  out
}
