# ---- pyramidal tract (per-bundle somatotopic model) -------------------------
#
# [DECISIONS D8] Each body-part parcel (leg, trunk, arm, hand, face) becomes
# one bundle whose centreline is a *mean pathway* swept along a common
# macro-axis:
#
#     c_k(s) = A(s) + W_u(s) q_k(s) e_u(s) + W_v(s) r_k(s) e_v(s)
#
# with (e_u, e_v, T) a moving orthonormal frame (v = u x T) twisted by psi(s),
# and q_k, r_k the somatotopic cross-section coordinates, W_u, W_v the section
# widths. The decussation is a smooth per-streamline sign flip phi(s) of q
# over a short interval, with an optional posterior detour bump(s); because it
# happens on single streamlines, no streamline ever splits (D1) and nothing
# fuses (C3) - interlacing is just displacing crossed/uncrossed copies by a
# centroid shift.

#' Smooth C1 interpolation of control points over a station grid
#'
#' Interpolates control values given at the station abscissae `stations` onto
#' the grid `s` using the smoothstep ramp `t^2 (3 - 2 t)` (C1: value and slope
#' agree at every station). Values before the first station equal the first
#' control value; values after the last equal the last. This is the same ramp
#' used by [generate_centerline()], extracted so station tables can drive any
#' smoothly-varying field (twist `psi`, widths `W`, or per-parcel positions).
#'
#' @param s Numeric grid of abscissae (typically arc-length fractions in
#'   `[0, 1]`).
#' @param stations Numeric station abscissae, strictly increasing, with
#'   `range(stations)` covering `range(s)` recommended.
#' @param vals Numeric control values, one per station (same length as
#'   `stations`).
#' @returns Numeric vector, same length as `s`.
#' @export
#' @examples
#' s <- seq(0, 1, length.out = 101)
#' psi <- smoothstep_interp(s, c(0, 0.5, 1), c(0, pi / 2, pi))
#' psi[1]        # 0
#' psi[101]      # pi
smoothstep_interp <- function(s, stations, vals) {
  stations <- as.numeric(stations)
  vals <- as.numeric(vals)
  if (length(stations) != length(vals)) {
    cli::cli_abort(c(
      "{.arg stations} and {.arg vals} must have the same length.",
      "i" = "Got {length(stations)} stations and {length(vals)} values."
    ))
  }
  if (any(diff(stations) <= 0)) {
    cli::cli_abort("{.arg stations} must be strictly increasing.")
  }

  out <- numeric(length(s))
  for (i in seq_len(length(stations) - 1L)) {
    a <- stations[i]
    b <- stations[i + 1L]
    inb <- s >= a & s <= b
    t <- (s[inb] - a) / (b - a)
    out[inb] <- vals[i] + (vals[i + 1L] - vals[i]) * (t * t * (3 - 2 * t))
  }
  # Clamp outside the last station (no extrapolation).
  out[s > stations[length(stations)]] <- vals[length(stations)]
  out
}

#' Anatomical stations of the pyramidal tract (cortex to cord)
#'
#' Ordered stations along the tract, as a fraction of the total arc length
#' (`fraction`, in `[0, 1]`) and as a depth in millimetres measured from the
#' motor-cortex surface along the pathway (`depth_mm`, at the default
#' `tract_length = 75` mm).
#'
#' @returns A data frame with one row per station and columns `station`
#'   (ordered factor), `fraction`, `depth_mm` (at the default
#'   `tract_length = 70` mm).
#' @export
#' @examples
#' pyt_stations()
pyt_stations <- function() {
  stations <- c(
    "cortex", "corona", "centrum", "ic", "crus",
    "pons", "pyramid", "dec", "cord"
  )
  fraction <- c(0, 0.10, 0.20, 0.35, 0.50, 0.60, 0.72, 0.82, 1.00)
  data.frame(
    station = factor(stations, levels = stations),
    fraction = fraction,
    depth_mm = round(70 * fraction, 1)
  )
}

#' Default pyramidal-tract parameter table
#'
#' Returns the full D8 parameter table used by
#' [simulate_pyramidal_tract()] when its arguments are `NULL`: the anatomical
#' stations, the five somatotopic parcels (`leg`, `trunk`, `arm`, `hand`,
#' `face`), their cross-section control points `q_st` (along `u`) and `r_st`
#' (along `v`) at each station, the section widths `width_u`/`width_v` (mm),
#' the twist map `twist`, per-parcel terminal stations (`face` ends at the
#' `pons`, everything else at the `cord`), calibers (axonal radius, mm), the
#' macro-axis `macros`, the decussation `cross` controls and the per-parcel
#' FDA `perturbation` controls (reused by [create_axons()]).
#'
#' Coordinates follow the amended validation convention: at the cortex,
#' positive `u` points toward the medial wall of the parcel's hemisphere and
#' parcels are spaced by their relative sheet extents (`leg -1.0, trunk -0.6,
#' arm -0.15, hand +0.35, face +0.85`).
#'
#' @returns A named list with elements `tract_length`, `stations`, `parcels`,
#'   `q_st`, `r_st`, `width_u`, `width_v`, `twist`, `term_station`, `caliber`,
#'   `macros`, `cross`, `perturbation`, as documented in DECISIONS D8.
#' @export
#' @examples
#' d <- pyt_defaults()
#' d$parcels                       # five body-part parcels
#' d$term_station["face"]          # "pons"
pyt_defaults <- function() {
  parcels <- c("leg", "trunk", "arm", "hand", "face")

  # From the crus onward `q_st` is single-signed (ascending leg -> face) so the
  # per-streamline sign flip `phi` alone decides the landing side: crossed
  # fibres land contralateral, uncrossed fibres stay ipsilateral, and within
  # each branch leg is most medial, face most lateral (D8 convention).
  q_st <- rbind(
    cortex  = c(-1.00, -0.60, -0.15,  0.35,  0.85),
    corona  = c(-1.00, -0.60, -0.15,  0.35,  0.85),
    centrum = c(-0.90, -0.50, -0.10,  0.38,  0.80),
    ic      = c( 0.90,  0.40, -0.10, -0.45, -0.85),
    crus    = c( 0.25,  0.35,  0.45,  0.55,  0.65),
    pons    = c( 0.25,  0.35,  0.45,  0.55,  0.75),
    pyramid = c( 0.30,  0.38,  0.46,  0.54,  0.68),
    dec     = c( 0.35,  0.45,  0.55,  0.65,  0.72),
    cord    = c( 0.35,  0.45,  0.55,  0.65,  0.75)
  )
  r_st <- rbind(
    cortex  = c( 0.55,  0.25,  0.00, -0.25, -0.55),
    corona  = c( 0.45,  0.20,  0.00, -0.20, -0.45),
    centrum = c( 0.22,  0.10,  0.00, -0.10, -0.22),
    ic      = rep(0, 5),
    crus    = c( 0.12,  0.05,  0.00, -0.05, -0.15),
    pons    = c( 0.12,  0.05,  0.00, -0.05, -0.15),
    pyramid = rep(0, 5),
    dec     = rep(0, 5),
    cord    = rep(0, 5)
  )
  colnames(q_st) <- parcels
  colnames(r_st) <- parcels

  list(
    tract_length = 70,                         # mm, cortex to cord
    stations     = pyt_stations(),
    parcels      = parcels,
    q_st         = q_st,
    r_st         = r_st,
    width_u      = setNames(c(40, 35, 22, 7, 10, 12, 6, 10, 6), rownames(q_st)),
    width_v      = setNames(c(10, 9, 8, 1.5, 3, 5, 3, 6, 4), rownames(q_st)),
    twist        = setNames(c(0, 0, 0.1, pi / 2, pi, pi, 0, 0, pi), rownames(q_st)),
    term_station = c(leg = "cord", trunk = "cord", arm = "cord",
                     hand = "cord", face = "pons"),
    caliber      = setNames(rep(0.002, 5), parcels),   # mm (1-4 um thin axons)
    macros       = list(tract_length = 70, midline_offset = 0,
                        bow = 1.5, n_points = 300L),
    cross        = list(cross_percent = 0.85, cross_span = 0.06,
                        cross_detour = 1.8, center_gap = 0.5),
    perturbation = list(variance = 0.02, autocorrelation = 2,
                        n_basis = 8L, n_axons = 30L)
  )
}

# Build the full macro-axis point grid for one hemisphere. X is the
# medial-lateral axis (0 = midline), Y the anterior axis (cortex is bow-anterior
# at the surface), Z the depth (cortex = 0, cord = -tract_length).
.macro_axis_grid <- function(macros, hemisphere = "left") {
  tl <- macros$tract_length       # mm cortex -> cord
  ml <- macros$midline_offset     # mm, lateral position of the tract centre
  bow <- macros$bow               # mm, anterior bow amplitude
  n  <- macros$n_points
  side <- if (hemisphere == "left") -1 else 1

  s <- seq(0, 1, length.out = n)
  list(
    s = s,
    X = rep(side * ml, n),        # tract centre: left at -ml, right at +ml
    Y = bow * sin(pi * s),        # anterior bow (cortex/cord at 0)
    Z = -tl * s                   # depth, cortex = 0
  )
}

# Moving frame: S_u = fixed X axis (u positive = medial for the left
# hemisphere), v = u x T. The right hemisphere mirrors across the YZ plane, so
# its u/v flip sign in X.
.moving_frame <- function(s, psi, hemisphere = "left") {
  m <- if (hemisphere == "left") 1 else -1
  list(
    ux = m * cos(psi),
    uy = sin(psi),
    uz = rep(0, length(psi)),
    vx = -m * sin(psi),
    vy = cos(psi),
    vz = rep(0, length(psi))
  )
}

#' Simulate one pyramidal tract (one hemisphere) as a `fiber::bundle_set`
#'
#' Generates the pyramidal tract of a single hemisphere as a
#' `fiber::bundle_set`: five per-parcel bundles (`leg`, `trunk`, `arm`,
#' `hand`, `face`) laid out by the somatotopic table [pyt_defaults()]. Each
#' parcel is a centreline swept along a common macro-axis with a smooth (C1)
#' moving frame, then perturbed into a bundle of axons with [create_axons()].
#'
#' The `face` parcel terminates at the pons (corticobulbar); the limb/trunk
#' parcels continue through the decussation into the cord. Each parcel
#' contains a crossed and an uncrossed branch: the crossed branch's `q` is
#' multiplied by a smooth sign-flip `phi(s)` (default 85% of fibers cross),
#' nudged laterally by `center_gap` (mm, so the two copies interlace at the
#' decussation) and pushed posteriorly by `cross_detour` (mm, the
#' interdigitating X), without any streamline splitting (D1) or fusion (C3).
#'
#' @param parcels Optional character vector of parcel names to include
#'   (defaults to all five from [pyt_defaults()]).
#' @param macros Optional list with `tract_length`, `midline_offset`, `bow`,
#'   `n_points` (defaults from [pyt_defaults()]).
#' @param rot_st Optional numeric control points for the twist map `psi`
#'   (length equal to `nrow(pyt_defaults()$q_st)`).
#' @param cross Optional list with `cross_percent`, `cross_span`,
#'   `cross_detour`, `center_gap`.
#' @param perturbation Optional list with `variance`, `autocorrelation`,
#'   `n_basis`, `n_axons` (reused by [create_axons()]).
#' @param hemisphere Either `"left"` or `"right"`; the medial wall faces the
#'   midline (`+X` toward the neighbour hemisphere).
#' @param seed Optional integer seed for full reproducibility.
#' @returns A `fiber::bundle_set` of one bundle per parcel. Each stream's
#'   `streamline_data` carries `parcel` and `branch` (`crossed`/`uncrossed`).
#' @export
#' @examples
#' pyt <- simulate_pyramidal_tract(hemisphere = "left", seed = 1)
#' fiber::is_bundle_set(pyt)      # TRUE
simulate_pyramidal_tract <- function(
    parcels = NULL,
    macros = NULL,
    rot_st = NULL,
    cross = NULL,
    perturbation = NULL,
    hemisphere = "left",
    seed = NULL
) {
  hemisphere <- match.arg(hemisphere, c("left", "right"))
  d <- pyt_defaults()
  if (is.null(parcels))      parcels      <- d$parcels
  if (is.null(macros))       macros       <- d$macros
  if (is.null(rot_st))       rot_st       <- d$twist
  if (is.null(cross))        cross        <- d$cross
  if (is.null(perturbation)) perturbation <- d$perturbation

  if (!is.null(seed)) set.seed(seed)
  n_st <- nrow(d$q_st)
  if (length(rot_st) != n_st) {
    cli::cli_abort(c(
      "{.arg rot_st} must have one twist value per station (got {length(rot_st)} for {n_st}).",
      "i" = "See {.fn pyt_defaults} for the station grid."
    ))
  }
  q_st <- d$q_st[, parcels, drop = FALSE]
  r_st <- d$r_st[, parcels, drop = FALSE]
  term_station <- d$term_station[parcels]

  axis <- .macro_axis_grid(macros, hemisphere = hemisphere)
  s <- axis$s
  n <- length(s)

  psi <- smoothstep_interp(s, d$stations$fraction, rot_st)
  Wu  <- smoothstep_interp(s, d$stations$fraction, d$width_u)
  Wv  <- smoothstep_interp(s, d$stations$fraction, d$width_v)
  fr  <- .moving_frame(s, psi, hemisphere = hemisphere)

  # C1 sign-flip for the crossed branch. Limbs/trunk flip just before the
  # `dec` station (contralateral cord course). The face is the corticobulbar
  # parcel: it flips right before the pons so both branches terminate there
  # bilaterally (D8), and does not take the posterior detour.
  dec_frac <- d$stations$fraction[d$stations$station == "dec"]
  pons_frac <- d$stations$fraction[d$stations$station == "pons"]
  span <- cross$cross_span

  branches <- c(uncrossed = FALSE, crossed = TRUE)
  sls <- vector("list", length(parcels) * 2L)
  names(sls) <- as.vector(outer(parcels, names(branches), paste, sep = "."))

  for (k in seq_along(parcels)) {
    q <- smoothstep_interp(s, d$stations$fraction, q_st[, k])
    r <- smoothstep_interp(s, d$stations$fraction, r_st[, k])
    keep <- s <= d$stations$fraction[d$stations$station == term_station[k]]

    centre <- if (term_station[k] == "pons") pons_frac else dec_frac
    # For the face the flip must *complete at* the terminal pons station, so
    # the window is [pons - span, pons]; elsewhere it straddles the centre.
    hi <- if (term_station[k] == "pons") pons_frac else centre + span / 2
    lo <- if (term_station[k] == "pons") pons_frac - span else centre - span / 2
    phi <- rep(1, n)
    inb <- s >= lo & s <= hi
    tt <- (s[inb] - lo) / (hi - lo)
    phi[inb] <- cos(pi * .smoothstep(tt))  # C1 taper +1 -> -1
    phi[s > hi] <- -1                      # crossed fibers stay contralateral

    # Posterior detour (interdigitating X): smooth bump centred at the flip.
    bump <- rep(0, n)
    if (term_station[k] != "pons") {
      bump[inb] <- cross$cross_detour * sin(pi * .smoothstep(tt))
    }

    for (bb in seq_along(branches)) {
      qq <- q
      rr <- r
      if (branches[bb]) {
        qq <- qq * phi                        # per-streamline sign flip
        rr <- rr + bump                       # posterior detour (X interlace)
        # Nudge the crossed copy laterally across the midline at the
        # decussation so the two copies visibly interlace (D8 centroid shift).
        qq <- qq + cross$center_gap / pmax(Wu, 1e-6)
      }
      pts <- cbind(
        X = axis$X[keep] + Wu[keep] * qq[keep] * fr$ux[keep] +
            Wv[keep] * rr[keep] * fr$vx[keep],
        Y = axis$Y[keep] + Wu[keep] * qq[keep] * fr$uy[keep] +
            Wv[keep] * rr[keep] * fr$vy[keep],
        Z = axis$Z[keep] + Wv[keep] * rr[keep] * fr$vz[keep]
      )
      name <- paste(parcels[k], names(branches)[bb], sep = ".")
      # Keep the branch bundle so per-stream metadata stays lifted (see
      # creation note: fiber::bundle() drops per-stream copies on wrap).
      sls[[name]] <- create_axons(
        pts, n_axons = perturbation$n_axons,
        n_basis = perturbation$n_basis,
        variance = perturbation$variance,
        autocorrelation = perturbation$autocorrelation,
        parcel = parcels[k],
        branch = names(branches)[bb],
        hemisphere = hemisphere
      )
      # Spread that visibly fills each parcel: rescale each axon's
      # displacement from its centreline by the local mean cross-section
      # half-width (W_u + W_v) / 4. With the default coefficient variance
      # (0.02), this yields a within-branch RMS deviation of ~0.3-1.6 mm
      # (cortex/corona wide, pyramid narrow), so bundles read as narrow but
      # contiguous clouds instead of hairline strokes.
      scale <- (smoothstep_interp(s, d$stations$fraction, d$width_u) +
                smoothstep_interp(s, d$stations$fraction, d$width_v)) / 4
      scale <- scale[keep]
      for (j in seq_len(sls[[name]]@n_streamlines)) {
        sls[[name]]@streamlines[[j]]@points <-
          cbind(
            X = pts[, "X"] + (sls[[name]]@streamlines[[j]]@points[, "X"] - pts[, "X"]) * scale,
            Y = pts[, "Y"] + (sls[[name]]@streamlines[[j]]@points[, "Y"] - pts[, "Y"]) * scale,
            Z = pts[, "Z"] + (sls[[name]]@streamlines[[j]]@points[, "Z"] - pts[, "Z"]) * scale
          )
      }
    }
  }

  bundles <- lapply(parcels, function(p) {
    ub <- sls[[paste(p, "uncrossed", sep = ".")]]
    cb <- sls[[paste(p, "crossed", sep = ".")]]
    sd_nms <- union(names(ub@streamline_data), names(cb@streamline_data))
    sd <- lapply(sd_nms, function(nm) {
      c(ub@streamline_data[[nm]], cb@streamline_data[[nm]])
    })
    names(sd) <- sd_nms
    fiber::bundle(
      streamlines = c(ub@streamlines, cb@streamlines),
      streamline_data = sd,
      bundle_data = list(role = "pyramidal", parcel = p, hemisphere = hemisphere)
    )
  })
  names(bundles) <- parcels
  fiber::bundle_set(bundles = bundles)
}
