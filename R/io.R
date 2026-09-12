#' Read a GenAlEx-formatted CSV
#'
#' Thin wrapper around [poppr::read.genalex()] that fails with a useful message
#' instead of a cryptic parse error when the file is not actually in GenAlEx
#' format.
#'
#' A GenAlEx CSV has three header rows -- number of loci, number of samples,
#' population sizes -- followed by one row per individual with two columns per
#' locus for a diploid.
#'
#' @param path Path to the CSV file.
#' @param ploidy Ploidy of the organism. Default 2.
#' @param genclone Return a `genclone` object rather than a `genind`?
#'   Default `FALSE`, which is what most downstream functions here expect.
#' @param ... Passed through to [poppr::read.genalex()].
#'
#' @return A `genind` (or `genclone`) object.
#' @export
#'
#' @examples
#' \dontrun{
#' gen <- read_genalex_csv("data/nj.csv")
#' gen
#' }
read_genalex_csv <- function(path, ploidy = 2, genclone = FALSE, ...) {
  if (!file.exists(path)) {
    stop("File not found: ", path, call. = FALSE)
  }

  gen <- tryCatch(
    poppr::read.genalex(path, ploidy = ploidy, genclone = genclone, ...),
    error = function(e) {
      stop(
        "Could not read '", basename(path), "' as a GenAlEx file.\n",
        "  A GenAlEx CSV needs three header rows (loci, samples, population\n",
        "  sizes) before the genotype rows. Original error:\n  ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  if (is.null(adegenet::pop(gen))) {
    warning(
      "No population factor was found. AMOVA and pairwise F-statistics ",
      "will not work until populations are assigned with adegenet::pop().",
      call. = FALSE
    )
  }

  gen
}

#' Write a plot to a file
#'
#' Convenience wrapper so scripts do not have to repeat the
#' `pdf(); plot(); dev.off()` dance, and so the device is always closed even if
#' plotting fails.
#'
#' @param filename Output path. The extension decides the device: `.pdf`,
#'   `.png`, `.tiff` or `.svg`.
#' @param plot_expr An expression that draws the plot, or a ggplot object.
#' @param width,height Dimensions in inches. Default 8 x 8.
#' @param res Resolution in dpi, used by raster devices only. Default 300.
#'
#' @return `filename`, invisibly.
#' @export
save_figure <- function(filename, plot_expr, width = 8, height = 8, res = 300) {
  ext <- tolower(tools::file_ext(filename))

  switch(
    ext,
    pdf  = grDevices::pdf(filename, width = width, height = height),
    png  = grDevices::png(filename, width = width, height = height,
                          units = "in", res = res),
    tiff = grDevices::tiff(filename, width = width, height = height,
                           units = "in", res = res),
    svg  = grDevices::svg(filename, width = width, height = height),
    stop("Unsupported file extension '", ext,
         "'. Use pdf, png, tiff or svg.", call. = FALSE)
  )

  on.exit(grDevices::dev.off(), add = TRUE)

  if (inherits(plot_expr, "ggplot")) {
    print(plot_expr)
  } else {
    force(plot_expr)
  }

  invisible(filename)
}
