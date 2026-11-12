# ---- plot dispatch ----------------------------------------------------------
#
# plot(<fiber object>) redirects to rtists::plot3d() when rtists is installed,
# so a simulated scenario is one call away from an interactive 3-D view.

plot_fiber_object <- function(x, ...) {
  if (requireNamespace("rtists", quietly = TRUE)) {
    rtists::plot3d(x, ...)
  } else {
    cat("Install the 'rtists' package to visualize ",
        sprintf("<fiber::%s>", class(x)[1]), " interactively.\n", sep = "")
    invisible(x)
  }
}

S7::method(plot, fiber::streamline) <- function(x, ...) plot_fiber_object(x, ...)
S7::method(plot, fiber::bundle) <- function(x, ...) plot_fiber_object(x, ...)
S7::method(plot, fiber::bundle_set) <- function(x, ...) plot_fiber_object(x, ...)
