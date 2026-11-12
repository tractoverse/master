# ---- export -----------------------------------------------------------------

#' Export a `fiber::bundle_set` (or `bundle`) to a VTP-to-CSV-compatible file
#'
#' Writes simulated tractography data to a plain-text CSV in the column layout
#' used by the master VTP-to-CSV tools. Every streamline becomes one line, with
#' the file name, a running point counter and the three coordinates
#' (`X`, `Y`, `Z`) of its points, plus an optional per-point attribute column.
#'
#' @param x A `fiber::bundle_set`, `fiber::bundle`, or `fiber::streamline`.
#' @param file Path to the output `.csv`.
#' @param attribute Optional name of a per-point attribute (e.g. `"radius"`,
#'   `"curvature"`) to append as an extra trailing column.
#' @returns The path `file`, invisibly.
#' @export
#' @examples
#' cl <- generate_centerline(angle = 0.5, n_points = 20)
#' b <- create_axons(cl, n_axons = 3)
#' out <- generate_csv(b, tempfile(fileext = ".csv"))
#' head(read.csv(out))
generate_csv <- function(x, file, attribute = NULL) {
  streams <- if (fiber::is_bundle_set(x)) {
    unlist(lapply(x@bundles, function(b) b@streamlines), recursive = FALSE)
  } else if (fiber::is_bundle(x)) {
    x@streamlines
  } else if (fiber::is_streamline(x)) {
    list(x)
  } else {
    cli::cli_abort("{.arg x} must be a fiber bundle_set, bundle or streamline.")
  }

  rows <- lapply(seq_along(streams), function(i) {
    sl <- streams[[i]]
    pts <- sl@points
    n <- nrow(pts)
    data.frame(
      File = rep(basename(file), n),
      Point = seq.int(n) - 1L,
      X = pts[, "X"],
      Y = pts[, "Y"],
      Z = pts[, "Z"]
    )
  })
  df <- do.call(rbind, rows)
  df$Tract <- rep(seq_along(streams), lengths(rows))

  if (!is.null(attribute)) {
    vals <- lapply(streams, function(sl) as.numeric(sl@point_data[[attribute]]))
    df[[attribute]] <- unlist(vals)
  }

  utils::write.csv(df, file, row.names = FALSE)
  invisible(file)
}
