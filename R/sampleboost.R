#' @title Boosted Sampling Optimization
#'
#' @description
#' Multi-objective optimization of field sampling design. Generates \code{perm}
#' random configurations of \code{nplot} non-overlapping circular plots within
#' a study area, scores each on up to three objectives, and returns the
#' best-scoring configuration with field-ready outputs.
#'
#' The three objectives, all maximized, are:
#' \itemize{
#'   \item environmental heterogeneity, as between-plot NDVI variance;
#'   \item floristic ignorance, as the mean MRFI value across plots
#'     (D'Antraccoli et al., 2022);
#'   \item spatial dispersion, as the mean nearest-neighbour distance between
#'     plot centres (Clark & Evans, 1954).
#' }
#'
#' Each objective can be switched off by setting its weight to 0.
#'
#' @param ndvi A \code{terra} \code{SpatRaster} holding NDVI or another
#'   environmental index.
#' @param ignorance A \code{terra} \code{SpatRaster} holding the Map of
#'   Relative Floristic Ignorance (MRFI), as returned by
#'   \code{\link{ignorance_map}}.
#' @param site An \code{sf} or \code{sfc} polygon defining the study area.
#'   Multiple features are dissolved before sampling.
#' @param excl_areas Optional \code{sf} or \code{sfc} polygon delimiting areas
#'   unsuitable for sampling (water bodies, inaccessible terrain). If no CRS is
#'   set, EPSG:4326 is assumed.
#' @param CRS.new EPSG code of a projected CRS in metres. If \code{NULL} (the
#'   default) the CRS of \code{ndvi} is used when projected, avoiding
#'   unnecessary reprojection.
#' @param nplot Integer, number of plots per configuration (minimum 2).
#' @param plot_radius Numeric, plot radius in metres.
#' @param perm Integer, number of configurations to generate and score.
#' @param ndvi.weight,igno.weight,dist.weight Non-negative numeric weights for
#'   the three objectives (defaults 1). Set a weight to 0 to disable the
#'   corresponding objective. At least one weight must be positive.
#' @param seed Optional integer seed for reproducibility. Recorded in
#'   \code{statistics}.
#' @param block_size Integer, configurations drawn per batch (default 250).
#'   Larger values are faster but use more memory.
#' @param verbose Logical, print progress messages (default \code{TRUE}).
#' @param output_dir Optional directory for CSV and PDF output. If \code{NULL}
#'   (default) no files are written.
#' @param output_prefix Character prefix for output filenames (default
#'   \code{"SampleBoost"}).
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{best_solution}}{Data frame of the winning plot centres with
#'     NDVI and ignorance values.}
#'   \item{\code{best_solution_sf}}{The same as an \code{sf} point object.}
#'   \item{\code{best_scores}}{One-row data frame of the winning
#'     configuration's raw and normalized scores.}
#'   \item{\code{scores}}{Data frame of all configurations.}
#'   \item{\code{field_sheet}}{Table with WGS84 coordinates for field use.}
#'   \item{\code{plots}}{Named list of \code{ggplot} objects.}
#'   \item{\code{statistics}}{Run parameters and summary.}
#' }
#'
#' @details
#' \strong{Extraction.} Each plot is represented by the raster cell containing
#' its centre. A warning is issued when plot area exceeds cell area.
#'
#' \strong{Batching.} Points for \code{block_size} configurations are drawn in
#' a single \code{sf::st_sample()} call and extracted in a single
#' \code{terra::extract()} call per raster. This removes per-call overhead that
#' otherwise dominates runtime.
#'
#' \strong{Non-overlap.} Configurations with any pair of centres closer than
#' \code{2 * plot_radius} are rejected. A warning is issued when fewer than
#' 10\% of configurations are valid.
#'
#' \strong{Boundary constraint.} Plot centres are drawn from the study area
#' eroded by \code{plot_radius}, so every plot lies entirely inside the site.
#'
#' \strong{Spatial dispersion.} The dispersion objective is the mean
#' nearest-neighbour distance (NND) across plot centres. NND rewards even
#' spacing across the study area: every plot must be well separated from its
#' closest neighbour, which prevents local redundancy. Under the distance-decay
#' of community similarity (Nekola & White, 1999; Soininen et al., 2007),
#' ensuring adequate separation at every point in the network maximises species
#' turnover across the sampling design. The use of NND as a dispersion metric
#' follows Clark & Evans (1954).
#'
#' \strong{Score semantics.} Each objective is min-max normalized across the
#' configurations in this run. \code{final_score} ranks configurations within a
#' run and is not comparable across runs or values of \code{perm}. When a
#' weight is 0, the corresponding objective contributes nothing to the final
#' score. Raw values are returned alongside normalized ones.
#'
#' @references
#' Clark, P.J. & Evans, F.C. (1954). Distance to nearest neighbour as a
#' measure of spatial relationships in populations. \emph{Ecology}, 35(4),
#' 445--453.
#' 
#' D’Antraccoli, M., Bedini, G. & Peruzzi, L. (2022). Maps of relative 
#' floristic ignorance and virtual floristic lists: An R package to incorporate 
#' uncertainty in mapping and analysing biodiversity data. \emph{Ecological 
#' Informatics}, 67, 101512.
#'
#' Nekola, J.C. & White, P.S. (1999). The distance decay of similarity in
#' biogeography and ecology. \emph{Journal of Biogeography}, 26, 867--878.
#'
#' Soininen, J., McDonald, R. & Hillebrand, H. (2007). The distance decay of
#' similarity in ecological communities. \emph{Ecography}, 30, 3--12.
#'
#' @importFrom sf st_area st_as_sf st_bbox st_buffer st_cast st_coordinates
#'   st_crs st_difference st_intersection st_intersects st_is_empty
#'   st_is_longlat st_make_valid st_sample st_transform st_union
#' @importFrom terra crop crs extract mask project res vect
#' @importFrom tidyterra geom_spatraster
#' @importFrom ggplot2 aes coord_flip element_blank element_text geom_col
#'   geom_point geom_sf geom_text ggplot ggtitle labs scale_color_distiller
#'   scale_fill_distiller scale_y_continuous theme theme_classic theme_minimal
#'   xlab ylab
#' @importFrom rlang .data
#' @importFrom stats dist var
#' @importFrom utils setTxtProgressBar txtProgressBar write.csv
#' @importFrom grDevices pdf dev.off
#' @importFrom grid grid.draw gpar textGrob
#' @importFrom gridExtra grid.arrange tableGrob
#' @export
#'
#' @examples
#' \dontrun{
#' ndvi <- load_ndvi_example()
#' mrfi <- load_mrfi_example()
#' data(park)
#'
#' # Basic usage (equal weights)
#' res <- sampleboost(
#'   ndvi = ndvi, ignorance = mrfi, site = park,
#'   nplot = 50, plot_radius = 5.64, perm = 1000, seed = 1
#' )
#'
#' res$best_scores
#' res$plots$ndvi
#'
#' # NDVI + distance only (ignorance off)
#' res2 <- sampleboost(
#'   ndvi = ndvi, ignorance = mrfi, site = park,
#'   nplot = 50, plot_radius = 5.64, perm = 1000, seed = 1,
#'   igno.weight = 0
#' )
#' }
sampleboost <- function(ndvi, ignorance, site,
                        excl_areas = NULL,
                        CRS.new = NULL,
                        nplot, plot_radius, perm,
                        ndvi.weight = 1, igno.weight = 1, dist.weight = 1,
                        seed = NULL,
                        block_size = 250L,
                        verbose = TRUE,
                        output_dir = NULL,
                        output_prefix = "SampleBoost") {
  
  msg <- function(...) if (verbose) message(...)
  start_time <- Sys.time()
  
  # --- 1. INPUT VALIDATION ---------------------------------------------------
  
  if (missing(ndvi) || missing(ignorance) || missing(site) ||
      missing(nplot) || missing(plot_radius) || missing(perm)) {
    stop("Missing required arguments: ndvi, ignorance, site, nplot, plot_radius, perm.")
  }
  if (!inherits(ndvi, "SpatRaster")) stop("'ndvi' must be a terra SpatRaster.")
  if (!inherits(ignorance, "SpatRaster")) stop("'ignorance' must be a terra SpatRaster.")
  
  if (!is.numeric(nplot) || length(nplot) != 1L || nplot < 2)
    stop("'nplot' must be a single number >= 2.")
  if (!is.numeric(plot_radius) || length(plot_radius) != 1L || plot_radius <= 0)
    stop("'plot_radius' must be a single positive number.")
  if (!is.numeric(perm) || length(perm) != 1L || perm < 1)
    stop("'perm' must be a single number >= 1.")
  if (!is.numeric(block_size) || length(block_size) != 1L || block_size < 1)
    stop("'block_size' must be a single number >= 1.")
  
  weights <- c(ndvi = ndvi.weight, igno = igno.weight, dist = dist.weight)
  if (!is.numeric(weights) || anyNA(weights) || any(weights < 0))
    stop("Weights must be non-negative numbers.")
  if (sum(weights) == 0)
    stop("At least one weight must be positive.")
  
  nplot <- as.integer(nplot)
  perm <- as.integer(perm)
  block_size <- min(as.integer(block_size), perm)
  
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L) stop("'seed' must be a single number.")
    set.seed(seed)
  }
  
  areaplot <- pi * plot_radius^2
  msg("Inputs validated.")
  
  # --- 2. WORKING CRS --------------------------------------------------------
  
  ndvi_crs <- sf::st_crs(terra::crs(ndvi))
  
  if (is.null(CRS.new)) {
    if (is.na(ndvi_crs))
      stop("'ndvi' has no CRS. Supply 'CRS.new' explicitly.")
    if (isTRUE(sf::st_is_longlat(ndvi_crs)))
      stop("'ndvi' is in geographic coordinates. Supply a projected 'CRS.new'.")
    crs_sf <- ndvi_crs
    msg(paste0("Working CRS from NDVI: ",
               if (is.na(crs_sf$epsg)) "(no EPSG)" else paste0("EPSG:", crs_sf$epsg), "."))
  } else {
    if (!is.numeric(CRS.new) || length(CRS.new) != 1L || is.na(CRS.new) || CRS.new <= 0)
      stop("'CRS.new' must be NULL or a positive numeric EPSG code.")
    crs_sf <- sf::st_crs(CRS.new)
    if (is.na(crs_sf)) stop("EPSG:", CRS.new, " not recognised.")
    if (isTRUE(sf::st_is_longlat(crs_sf)))
      stop("'CRS.new' must be a projected CRS in metres.")
    msg(paste0("Working CRS: EPSG:", CRS.new, "."))
  }
  
  # Reproject rasters if needed
  if (is.na(ndvi_crs) || ndvi_crs != crs_sf) {
    msg("  Reprojecting NDVI ...")
    ndvi <- terra::project(ndvi, crs_sf$wkt)
  }
  igno_crs <- sf::st_crs(terra::crs(ignorance))
  if (is.na(igno_crs) || igno_crs != crs_sf) {
    msg("  Reprojecting ignorance ...")
    ignorance <- terra::project(ignorance, crs_sf$wkt)
  }
  
  # --- 3. SAMPLING AREA ------------------------------------------------------
  
  if (inherits(site, "Spatial")) site <- sf::st_as_sf(site)
  if (is.na(sf::st_crs(site))) {
    msg("  Site has no CRS; assuming EPSG:4326.")
    sf::st_crs(site) <- 4326
  }
  site_proj <- sf::st_union(sf::st_transform(sf::st_make_valid(site), crs_sf))
  
  msg(paste0("Eroding study area by ", round(plot_radius, 1), " m ..."))
  site_buffered <- sf::st_buffer(site_proj, dist = -plot_radius)
  if (length(site_buffered) == 0L || all(sf::st_is_empty(site_buffered)))
    stop("Negative buffer of ", round(plot_radius, 1), " m leaves no area. Reduce 'plot_radius'.")
  
  original_area <- as.numeric(sum(sf::st_area(site_proj)))
  buffered_area <- as.numeric(sum(sf::st_area(site_buffered)))
  area_loss_pct <- (1 - buffered_area / original_area) * 100
  if (area_loss_pct > 15)
    warning("Boundary constraint removes ", round(area_loss_pct, 1),
            "% of the study area. Consider reducing 'plot_radius'.")
  
  # Exclusion areas
  excl_proj <- NULL
  has_exclusions <- !is.null(excl_areas)
  if (has_exclusions) {
    msg("  Processing exclusion areas ...")
    if (inherits(excl_areas, "Spatial")) excl_areas <- sf::st_as_sf(excl_areas)
    if (is.na(sf::st_crs(excl_areas))) {
      msg("  Exclusion areas have no CRS; assuming EPSG:4326.")
      sf::st_crs(excl_areas) <- 4326
    }
    excl_proj <- sf::st_union(sf::st_transform(sf::st_make_valid(excl_areas), crs_sf))
    site_sampling <- sf::st_make_valid(sf::st_difference(site_buffered, excl_proj))
    if (length(site_sampling) == 0L || all(sf::st_is_empty(site_sampling)))
      stop("Exclusion areas leave no area for sampling.")
  } else {
    site_sampling <- site_buffered
  }
  
  sampling_area <- as.numeric(sum(sf::st_area(site_sampling)))
  msg(paste0("  Sampling area: ", round(sampling_area / 1e6, 2), " km2 (",
             round(area_loss_pct, 1), "% lost to boundary constraint)."))
  
  # Check plot vs cell size
  cell_area <- prod(terra::res(ndvi))
  if (areaplot > cell_area)
    warning("Plot area (", round(areaplot), " m2) > NDVI cell area (",
            round(cell_area), " m2). Each plot is represented by its centre cell only.")
  msg(paste0("  Plot area: ", round(areaplot, 1), " m2; NDVI cell: ",
             round(cell_area, 1), " m2."))
  
  # Check packing density
  packing_ratio <- (nplot * areaplot) / sampling_area
  if (packing_ratio > 0.3)
    warning("Plots occupy ~", round(packing_ratio * 100),
            "% of the sampling area. Valid configurations will be rare.")
  
  # Crop and mask rasters to site
  site_vect <- terra::vect(sf::st_as_sf(site_proj))
  ndvi <- terra::mask(terra::crop(ndvi, site_vect), site_vect)
  ignorance <- terra::mask(terra::crop(ignorance, site_vect), site_vect)
  
  # --- 4. GENERATE AND SCORE CONFIGURATIONS -----------------------------------
  
  msg(paste0("Generating ", perm, " configurations (block size ", block_size, ") ..."))
  
  min_sep <- 2 * plot_radius
  
  ndvi_between_var <- rep(NA_real_, perm)
  mean_ignorance   <- rep(NA_real_, perm)
  mean_nn_dist     <- rep(NA_real_, perm)
  
  reject_overlap <- logical(perm)
  reject_npoints <- logical(perm)
  reject_na      <- logical(perm)
  configs        <- vector("list", perm)
  
  # Robust point drawing with retry
  draw_points <- function(n) {
    got <- list()
    have <- 0L
    tries <- 0L
    while (have < n && tries < 25L) {
      s <- sf::st_sample(site_sampling, size = n - have, type = "random")
      s <- suppressWarnings(sf::st_cast(s, "POINT"))
      if (length(s) > 0L) {
        got[[length(got) + 1L]] <- s
        have <- have + length(s)
      }
      tries <- tries + 1L
    }
    if (have < n)
      stop("st_sample() could not draw ", n, " points after 25 attempts.")
    do.call(c, got)[seq_len(n)]
  }
  
  block_starts <- seq(1L, perm, by = block_size)
  pb <- if (verbose) utils::txtProgressBar(min = 0, max = perm, style = 3) else NULL
  
  for (bs in block_starts) {
    
    be <- min(bs + block_size - 1L, perm)
    n_conf <- be - bs + 1L
    n_pts <- n_conf * nplot
    
    pts <- draw_points(n_pts)
    pts_sf <- sf::st_as_sf(pts)
    
    coords_all <- sf::st_coordinates(pts)
    ndvi_all <- terra::extract(ndvi, pts_sf, ID = FALSE)[, 1]
    igno_all <- terra::extract(ignorance, pts_sf, ID = FALSE)[, 1]
    
    if (has_exclusions) {
      excl_hit <- apply(sf::st_intersects(pts_sf, excl_proj, sparse = FALSE), 1, any)
    } else {
      excl_hit <- rep(FALSE, n_pts)
    }
    
    idx <- matrix(seq_len(n_pts), nrow = nplot)
    
    for (j in seq_len(n_conf)) {
      
      i <- bs + j - 1L
      rows <- idx[, j]
      
      if (any(excl_hit[rows])) { reject_npoints[i] <- TRUE; next }
      
      coords <- coords_all[rows, , drop = FALSE]
      dmat <- as.matrix(stats::dist(coords))
      diag(dmat) <- Inf
      
      # Overlap check: centres closer than 2 * plot_radius
      if (min(dmat) < min_sep) { reject_overlap[i] <- TRUE; next }
      
      ndvi_vals <- ndvi_all[rows]
      igno_vals <- igno_all[rows]
      if (anyNA(ndvi_vals) || anyNA(igno_vals)) { reject_na[i] <- TRUE; next }
      
      # Mean NND: rewards even spacing across the study area
      nn <- apply(dmat, 1, min)
      
      ndvi_between_var[i] <- stats::var(ndvi_vals)
      mean_ignorance[i]   <- mean(igno_vals)
      mean_nn_dist[i]     <- mean(nn)
      
      configs[[i]] <- data.frame(
        config_id = i,
        plot_id   = seq_len(nplot),
        x = coords[, 1], y = coords[, 2],
        ndvi = ndvi_vals, ignorance = igno_vals,
        nn_dist = nn,
        stringsAsFactors = FALSE
      )
    }
    
    if (!is.null(pb)) utils::setTxtProgressBar(pb, be)
  }
  if (!is.null(pb)) close(pb)
  
  valid <- !(reject_overlap | reject_npoints | reject_na)
  
  msg(paste0("  Valid: ", sum(valid), "/", perm,
             " (overlap: ", sum(reject_overlap),
             ", exclusion: ", sum(reject_npoints),
             ", NA: ", sum(reject_na), ")"))
  
  if (!any(valid))
    stop("No valid configuration in ", perm, " attempts. ",
         "Reduce 'nplot'/'plot_radius' or increase 'perm'.")
  if (mean(valid) < 0.10)
    warning("Only ", round(mean(valid) * 100, 1),
            "% valid. Increase 'perm' for a better solution.")
  
  # --- 5. NORMALIZE, WEIGHT, SELECT -------------------------------------------
  
  msg("Scoring ...")
  
  normalize <- function(x) {
    ok <- !is.na(x)
    if (!any(ok)) return(x)
    rng <- max(x[ok]) - min(x[ok])
    if (rng == 0) return(ifelse(ok, 0.5, NA_real_))
    (x - min(x[ok])) / rng
  }
  
  scores <- data.frame(
    config_id        = seq_len(perm),
    ndvi_between_var = ndvi_between_var,
    mean_ignorance   = mean_ignorance,
    mean_nn_dist     = mean_nn_dist,
    valid = valid,
    reject_overlap = reject_overlap,
    reject_npoints = reject_npoints,
    reject_na      = reject_na,
    stringsAsFactors = FALSE
  )
  
  # Zero-weight objectives contribute nothing
  scores$ndvi_norm <- if (ndvi.weight == 0) 0 else normalize(scores$ndvi_between_var)
  scores$igno_norm <- if (igno.weight == 0) 0 else normalize(scores$mean_ignorance)
  scores$dist_norm <- if (dist.weight == 0) 0 else normalize(scores$mean_nn_dist)
  
  scores$final_score <- (weights[["ndvi"]] * scores$ndvi_norm +
                           weights[["igno"]] * scores$igno_norm +
                           weights[["dist"]] * scores$dist_norm) / sum(weights)
  
  valid_scores <- scores[scores$valid, ]
  best <- valid_scores[which.max(valid_scores$final_score), ]
  
  best_solution <- configs[[best$config_id]]
  best_solution_sf <- sf::st_as_sf(best_solution, coords = c("x", "y"), crs = crs_sf)
  best_buffers_sf  <- sf::st_buffer(best_solution_sf, dist = plot_radius)
  
  msg(paste0("Best configuration: #", best$config_id,
             " (score ", round(best$final_score, 3), ")."))
  
  # --- 6. FIELD SHEET AND STATISTICS ------------------------------------------
  
  wgs <- sf::st_coordinates(sf::st_transform(best_solution_sf, 4326))
  field_sheet <- data.frame(
    Plot_ID         = sprintf("PLOT_%03d", best_solution$plot_id),
    Longitude_WGS84 = round(wgs[, 1], 6),
    Latitude_WGS84  = round(wgs[, 2], 6),
    NDVI            = round(best_solution$ndvi, 3),
    Ignorance       = round(best_solution$ignorance, 2),
    Notes           = "",
    stringsAsFactors = FALSE
  )
  
  end_time <- Sys.time()
  elapsed  <- as.numeric(difftime(end_time, start_time, units = "secs"))
  epsg_used <- if (is.na(crs_sf$epsg)) "custom" else as.character(crs_sf$epsg)
  
  statistics <- data.frame(
    Statistic = c("Started", "Finished", "Elapsed time (s)", "Working CRS (EPSG)",
                  "Seed", "Exclusion areas", "Number of plots", "Plot radius (m)",
                  "Plot area (m2)", "NDVI cell size (m)", "Sampling area (km2)",
                  "Configurations tested", "Block size", "Valid configurations",
                  "Rejected - overlap", "Rejected - in exclusion area",
                  "Rejected - missing values",
                  "NDVI weight", "Ignorance weight", "Distance weight",
                  "Best configuration", "Best final score",
                  "Best NDVI between-plot variance", "Best mean ignorance",
                  "Best mean nearest-neighbour distance (m)"),
    Value = c(format(start_time), format(end_time), round(elapsed, 2),
              epsg_used,
              if (is.null(seed)) "not set" else as.character(seed),
              if (has_exclusions) "yes" else "no",
              nplot, round(plot_radius, 2), round(areaplot, 1),
              round(sqrt(cell_area), 1), round(sampling_area / 1e6, 2),
              perm, block_size, sum(valid),
              sum(reject_overlap), sum(reject_npoints), sum(reject_na),
              weights[["ndvi"]], weights[["igno"]], weights[["dist"]],
              best$config_id, round(best$final_score, 4),
              signif(best$ndvi_between_var, 4), round(best$mean_ignorance, 3),
              round(best$mean_nn_dist, 1)),
    stringsAsFactors = FALSE
  )
  
  # --- 7. PLOTS ---------------------------------------------------------------
  
  excl_plot <- NULL
  if (has_exclusions) {
    excl_plot <- suppressWarnings(sf::st_intersection(excl_proj, site_proj))
    if (length(excl_plot) == 0L || all(sf::st_is_empty(excl_plot))) excl_plot <- NULL
  }
  
  add_boundaries <- function(p) {
    p <- p + ggplot2::geom_sf(data = site_proj, fill = NA, colour = "black",
                              linewidth = 0.8, inherit.aes = FALSE)
    if (!is.null(excl_plot))
      p <- p + ggplot2::geom_sf(data = excl_plot, fill = "grey60", colour = "black",
                                linewidth = 0.4, linetype = "dotted",
                                alpha = 0.3, inherit.aes = FALSE)
    p
  }
  
  map_layer <- function(r, palette, direction, legend, title) {
    p <- ggplot2::ggplot() +
      tidyterra::geom_spatraster(data = r) +
      ggplot2::geom_sf(data = best_buffers_sf, fill = NA, colour = "white", linewidth = 1.2) +
      ggplot2::geom_sf(data = best_solution_sf, colour = "white", size = 2.2,
                       shape = 21, fill = "black", stroke = 1) +
      ggplot2::scale_fill_distiller(palette = palette, name = legend,
                                    na.value = "transparent", direction = direction) +
      ggplot2::ggtitle(title) + ggplot2::theme_minimal()
    add_boundaries(p)
  }
  
  plot_ndvi      <- map_layer(ndvi, "YlGn", 1, "NDVI", "Best configuration on NDVI")
  plot_ignorance <- map_layer(ignorance, "Spectral", -1, "MRFI", "Best configuration on MRFI")
  
  plot_scores <- ggplot2::ggplot(
    valid_scores,
    ggplot2::aes(x = .data$ndvi_norm, y = .data$igno_norm,
                 size = .data$dist_norm, colour = .data$final_score)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_point(data = best, colour = "black", size = 5, shape = 18) +
    ggplot2::scale_color_distiller(palette = "Spectral", name = "Final score") +
    ggplot2::ggtitle("Objective space") +
    ggplot2::xlab("NDVI variance (normalized)") +
    ggplot2::ylab("Mean ignorance (normalized)") +
    ggplot2::labs(size = "Dispersion") +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold"))
  
  contributions <- data.frame(
    Objective = factor(c("Dispersion", "Ignorance", "NDVI variance"),
                       levels = c("Dispersion", "Ignorance", "NDVI variance")),
    Value = c(best$dist_norm, best$igno_norm, best$ndvi_norm) *
      c(weights[["dist"]], weights[["igno"]], weights[["ndvi"]]) / sum(weights)
  )
  
  plot_contributions <- ggplot2::ggplot(
    contributions, ggplot2::aes(x = .data$Objective, y = .data$Value)) +
    ggplot2::geom_col(fill = "#4C9F70", alpha = 0.85, width = 0.65) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", .data$Value)),
                       hjust = -0.15, size = 3.5) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(limits = c(0, max(contributions$Value) * 1.25)) +
    ggplot2::labs(
      title = sprintf("Weighted contributions (final score: %.3f)", best$final_score),
      x = NULL, y = "Contribution") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank())
  
  plots <- list(ndvi = plot_ndvi, ignorance = plot_ignorance,
                scores = plot_scores, contributions = plot_contributions)
  
  # --- 8. OPTIONAL FILE OUTPUT ------------------------------------------------
  
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    utils::write.csv(best_solution,
                     file.path(output_dir, paste0(output_prefix, "_best-solution.csv")),
                     row.names = FALSE)
    utils::write.csv(scores,
                     file.path(output_dir, paste0(output_prefix, "_all-scores.csv")),
                     row.names = FALSE)
    utils::write.csv(field_sheet,
                     file.path(output_dir, paste0(output_prefix, "_field-sheet.csv")),
                     row.names = FALSE)
    
    bbox <- sf::st_bbox(site_proj)
    wide <- (bbox[["xmax"]] - bbox[["xmin"]]) / (bbox[["ymax"]] - bbox[["ymin"]]) > 1.3
    map_ncol <- if (wide) 1 else 2
    
    grDevices::pdf(file.path(output_dir, paste0(output_prefix, "_plots.pdf")),
                   width = 11.69, height = 8.27, onefile = TRUE)
    on.exit(grDevices::dev.off(), add = TRUE)
    
    gridExtra::grid.arrange(
      plot_ndvi, plot_ignorance, ncol = map_ncol,
      top = grid::textGrob("Best sampling configuration",
                           gp = grid::gpar(fontsize = 14, fontface = "bold")))
    gridExtra::grid.arrange(
      plot_contributions, plot_scores, ncol = 1, heights = c(0.45, 1),
      top = grid::textGrob("Optimization diagnostics",
                           gp = grid::gpar(fontsize = 14, fontface = "bold")))
    grid::grid.draw(gridExtra::grid.arrange(
      gridExtra::tableGrob(statistics, rows = NULL),
      top = grid::textGrob("Summary statistics",
                           gp = grid::gpar(fontsize = 14, fontface = "bold"))))
    
    msg(paste0("Files written to ", output_dir, "."))
  }
  
  # --- 9. RETURN --------------------------------------------------------------
  
  list(
    best_solution    = best_solution,
    best_solution_sf = best_solution_sf,
    best_scores      = best,
    scores           = scores,
    field_sheet      = field_sheet,
    plots            = plots,
    statistics       = statistics
  )
}