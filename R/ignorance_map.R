#' @title Map of Relative Floristic Ignorance
#'
#' @description Computes a Map of Relative Floristic Ignorance (MRFI) using
#' modern `sf` and `terra` packages for high performance. Aligned with MRFI 
#' philosophy: shows all areas where data exists, letting IRFI values indicate 
#' sampling quality rather than excluding poorly-sampled areas.
#' 
#' The function automatically generates both quantile and continuous visualizations,
#' providing comprehensive analysis in a standardized 4-page PDF output.
#'
#' @param data_flor A data.frame with 5 columns: 'Taxon' (species identity), 
#' 'Long' (longitude), 'Lat' (latitude), 'uncertainty' (spatial uncertainty 
#' radius in meters), and 'year' (year of occurrence record).
#' @param site An `sf` object (or 'SpatialPolygonsDataFrame') representing the 
#' study area. If no CRS is set, assumes EPSG:4326 (WGS84).
#' @param year_study Numeric year of analysis (e.g., 2025). Defaults to current 
#' system year if not specified.
#' @param excl_areas Optional `sf` object (or 'SpatialPolygonsDataFrame') to 
#' delimit unsuitable areas (e.g., marine surfaces for terrestrial flora) to be 
#' excluded from calculations. If no CRS is set, assumes EPSG:4326.
#' @param CRS.new Numeric EPSG code for projected CRS used in calculations (must 
#' be in meters). Default = 3035 (ETRS89-LAEA Europe).
#' @param tau Numeric. Percentual value of taxa loss over 100 years (0 <= tau < 100). 
#' Represents the temporal decay of floristic records.
#' @param cellsize Numeric. Resolution of output raster in meters.
#' @param verbose Logical. If TRUE (default), prints progress messages during execution.
#' @param check_overlap Logical. If TRUE (default), checks and plots spatial overlap 
#' between occurrence points and study area.
#' @param output_dir Character. Directory path for output files. Defaults to working directory.
#' @param output_prefix Character. Prefix for output filenames. Default = "MRFI".
#' @param site_buffer Logical. If TRUE, expands the study area boundary to analyze 
#' a larger region beyond the original site. If FALSE (default), analyzes only the 
#' original boundary. Note: Edge effects from external records are already handled 
#' by the spatial uncertainty mechanism - points outside the site automatically 
#' contribute to cells within the study area if their uncertainty buffers reach in. 
#' This parameter is for expanding the analysis region itself, not for fixing edge 
#' effects. Buffer width controlled by buffer_width parameter.
#' @param buffer_width Numeric. Width of buffer in meters when site_buffer=TRUE. 
#' If NULL (default), uses cellsize as buffer distance. If specified, overrides 
#' default. Must be > 0. Buffer will NOT extend into exclusion areas.
#' @param mask_method Character. Method for final raster masking. One of:
#'   \itemize{
#'     \item{"touches"}{ Include any cell intersecting site boundary (default, MRFI-aligned)}
#'     \item{"none"}{ Include entire raster extent, no masking (shows all data)}
#'   }
#' @param use_coverage_weighting Logical. If TRUE (default), weights ignorance 
#' contribution by fraction of cell covered by each buffer (matches original algorithm, 
#' slower but more accurate). If FALSE, uses faster binary approach where any buffer 
#' touch contributes full ignorance value (faster but may underestimate ignorance).
#'
#' @return A list with 4 objects:
#' \itemize{
#'   \item{MRFI}{ A `terra SpatRaster` of the Map of Relative Floristic Ignorance}
#'   \item{RICH}{ A `terra SpatRaster` of species richness computed without uncertainties}
#'   \item{Uncertainties}{ A data.frame of uncertainty and year values for all records 
#'   used in computation}
#'   \item{Statistics}{ A data.frame summarizing settings and results}
#' }
#' 
#' @details 
#' PDF Output Structure (4 pages, A4 landscape):
#' \itemize{
#'   \item{Page 1}{ Quantile comparison (MRFI + Richness) - Primary analysis figure}
#'   \item{Page 2}{ Continuous comparison (MRFI + Richness) - Raw data view}
#'   \item{Page 3}{ Data diagnostics (temporal + spatial uncertainty histograms)}
#'   \item{Page 4}{ Summary statistics table}
#' }
#'
#' @importFrom sf st_as_sf st_make_valid st_crs st_transform st_buffer st_intersects st_intersection st_union st_bbox st_coordinates st_geometry st_difference st_sf
#' @importFrom terra rast values mask crop vect rasterize extract ext global writeRaster app
#' @importFrom utils txtProgressBar setTxtProgressBar write.csv
#' @importFrom ggplot2 ggplot aes geom_tile geom_sf coord_equal theme_classic labs theme xlab ylab scale_fill_gradientn guide_legend ggtitle geom_histogram coord_cartesian scale_y_continuous after_stat element_text
#' @importFrom grid unit grid.draw textGrob gpar
#' @importFrom gridExtra grid.arrange tableGrob
#' @importFrom grDevices pdf dev.off
#' @importFrom stats quantile
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales rescale
#' 
#' @export
#' 
#' @examples \dontrun{
#' # Load example data (requires 'ignobioR' package data)
#' data(floratus)
#' data(park)
#' data(unsuitablezone)
#'
#' # Example 1: Basic usage (recommended - edge effects already handled)
#' mrfi <- ignorance_map(
#'   data_flor = floratus,
#'   site = park,
#'   excl_areas = unsuitablezone,
#'   tau = 20,
#'   cellsize = 2000
#' )
#'
#' # Example 2: With buffer to EXPAND analysis area beyond original site
#' # (defaults to cellsize width expansion)
#' mrfi_buffered <- ignorance_map(
#'   data_flor = floratus,
#'   site = park,
#'   tau = 20,
#'   cellsize = 2000,
#'   site_buffer = TRUE  # Analyzes original site + 2000m buffer zone
#' )
#'
#' # Example 3: Custom buffer to analyze larger region (5km beyond site)
#' mrfi_5km <- ignorance_map(
#'   data_flor = floratus,
#'   site = park,
#'   excl_areas = unsuitablezone,
#'   tau = 20,
#'   cellsize = 2000,
#'   site_buffer = TRUE,
#'   buffer_width = 5000  # Analyzes original site + 5km buffer zone
#' )
#'
#' # Example 4: Fast mode without coverage weighting
#' mrfi_fast <- ignorance_map(
#'   data_flor = floratus,
#'   site = park,
#'   tau = 20,
#'   cellsize = 2000,
#'   use_coverage_weighting = FALSE
#' )
#'
#' # View results
#' terra::plot(mrfi$MRFI)
#' terra::plot(mrfi$RICH)
#' print(mrfi$Statistics)
#' }

ignorance_map <- function(data_flor, site, year_study = NULL, excl_areas = NULL,
                          CRS.new = 3035, tau, cellsize, verbose = TRUE,
                          check_overlap = TRUE, output_dir = file.path(getwd(), "output"),
                          output_prefix = "MRFI", site_buffer = FALSE,
                          buffer_width = NULL,
                          mask_method = "touches",
                          use_coverage_weighting = TRUE) {
  
  # ============================================================================
  # HELPER FUNCTIONS
  # ============================================================================
  
  msg <- function(...) if (verbose) message(...)
  start_time <- Sys.time()
  
  # ============================================================================
  # SECTION 1: SETTINGS AND INPUT VALIDATION
  # ============================================================================
  
  msg("Checking settings and inputs...")
  
  if (is.null(year_study)) year_study <- as.numeric(format(Sys.Date(), "%Y"))
  
  if (year_study < 1800 || year_study > 2100) {
    warning("year_study seems unusual. Double-check the value.")
  }
  
  if (any(is.na(data_flor$year))) {
    stop("Missing values (NA) found in 'year' column. Please remove or impute missing years.")
  }
  
  if (max(data_flor$year, na.rm = TRUE) > year_study) {
    warning("Some occurrence dates are more recent than year_study.")
  }
  
  if (tau < 0 || tau >= 100) stop("0 <= tau < 100 is required.")
  if (cellsize <= 0) stop("cellsize must be positive.")
  if (CRS.new <= 0 || !is.numeric(CRS.new)) stop("CRS.new must be a valid positive EPSG code.")
  if (!is.logical(site_buffer)) stop("site_buffer must be TRUE or FALSE.")
  
  if (site_buffer) {
    if (is.null(buffer_width)) {
      buffer_distance <- cellsize
      msg(paste0("Site buffer enabled: expanding analysis area by ", cellsize, "m (cellsize)..."))
    } else {
      if (!is.numeric(buffer_width) || buffer_width <= 0) {
        stop("buffer_width must be a positive numeric value (in meters).")
      }
      buffer_distance <- buffer_width
      msg(paste0("Site buffer enabled: expanding analysis area by ", buffer_distance, "m..."))
    }
    use_buffer <- TRUE
  } else {
    if (!is.null(buffer_width)) {
      warning("buffer_width specified but site_buffer=FALSE. Buffer width will be ignored.")
    }
    buffer_distance <- 0
    use_buffer <- FALSE
    msg("Analyzing original site boundary only (edge effects handled via uncertainty inclusion)...")
  }
  
  valid_methods <- c("touches", "none")
  if (!mask_method %in% valid_methods) {
    stop("mask_method must be one of: ", paste(valid_methods, collapse = ", "), ".")
  }
  
  if (!is.logical(use_coverage_weighting)) {
    stop("use_coverage_weighting must be TRUE or FALSE.")
  }
  
  if (use_coverage_weighting) {
    msg("Coverage weighting: ENABLED (accurate, matches original algorithm, slower)...")
  } else {
    msg("Coverage weighting: DISABLED (faster, binary touch method)...")
  }
  
  req_cols <- c("Taxon", "Long", "Lat", "uncertainty", "year")
  if (!all(req_cols %in% names(data_flor))) {
    stop("data_flor must contain columns: ", paste(req_cols, collapse = ", "), ".")
  }
  
  # Check for missing values in critical columns
  if (any(is.na(data_flor$uncertainty))) {
    stop("Missing values (NA) found in 'uncertainty' column. Please remove or impute missing values.")
  }
  
  if (any(is.na(data_flor$Lat)) || any(is.na(data_flor$Long))) {
    stop("Missing values (NA) found in 'Lat' or 'Long' columns. Please remove records with missing coordinates.")
  }
  
  total_initial_records <- nrow(data_flor)
  
  if (any(2 * data_flor$uncertainty < (cellsize / 20))) {
    warning("Some records have uncertainty very small relative to cellsize. ",
            "They may not be well-captured by the raster grid.")
  }
  
  msg("Inputs validated.")
  
  # ============================================================================
  # SECTION 2: COORDINATE SYSTEM SETUP AND REPROJECTION
  # ============================================================================
  
  msg(paste("Reprojecting inputs to EPSG:", CRS.new, "..."))
  
  crs_sf <- sf::st_crs(CRS.new)
  crs_terra <- paste0("EPSG:", CRS.new)
  
  if (inherits(site, "Spatial")) site <- sf::st_as_sf(site)
  
  if (is.na(sf::st_crs(site))) {
    msg("Input 'site' has no CRS. Assuming EPSG:4326.")
    sf::st_crs(site) <- 4326
  }
  
  site_proj_original <- sf::st_transform(sf::st_make_valid(site), crs_sf)
  
  excl_proj <- NULL
  has_exclusions <- !is.null(excl_areas)
  
  if (has_exclusions) {
    msg("Processing exclusion areas...")
    
    if (inherits(excl_areas, "Spatial")) excl_areas <- sf::st_as_sf(excl_areas)
    
    if (is.na(sf::st_crs(excl_areas))) {
      msg("Input 'excl_areas' has no CRS. Assuming EPSG:4326.")
      sf::st_crs(excl_areas) <- 4326
    }
    
    excl_proj <- sf::st_union(sf::st_transform(sf::st_make_valid(excl_areas), crs_sf))
  }
  
  # ============================================================================
  # SECTION 3: SITE BUFFER LOGIC (FOR EXPANDING ANALYSIS AREA)
  # ============================================================================
  
  if (use_buffer) {
    msg(paste0("Expanding study area boundary for analysis (buffer = ", buffer_distance, "m)..."))
    
    site_proj_buffered <- sf::st_buffer(site_proj_original, dist = buffer_distance)
    
    if (has_exclusions) {
      msg("  Removing exclusion areas from expanded site...")
      site_proj_buffered <- sf::st_difference(site_proj_buffered, excl_proj)
    }
    
    site_proj_buffered <- sf::st_make_valid(site_proj_buffered)
    site_proj_processing <- site_proj_buffered
    site_proj_mask <- site_proj_buffered
    
  } else {
    msg("Using original site boundary (no expansion)...")
    site_proj_processing <- site_proj_original
    site_proj_mask <- site_proj_original
    site_proj_buffered <- NULL
  }
  
  site_vect_processing <- terra::vect(site_proj_processing)
  site_vect_mask <- terra::vect(site_proj_mask)
  
  # ============================================================================
  # SECTION 4: FLORISTIC DATA PREPARATION
  # ============================================================================
  
  msg("Standardizing input: data_flor...")
  
  if (inherits(data_flor, "sf")) {
    msg("data_flor is already an sf object.")
    pts_sf <- data_flor
  } else if (inherits(data_flor, "Spatial")) {
    msg("Converting sp object to sf...")
    pts_sf <- sf::st_as_sf(data_flor)
  } else if (is.data.frame(data_flor)) {
    msg("Converting data.frame to sf (using Long/Lat)...")
    required_cols <- c("Long", "Lat")
    if (!all(required_cols %in% names(data_flor))) {
      stop("If 'data_flor' is a data.frame, it must include 'Long' and 'Lat' columns.")
    }
    pts_sf <- sf::st_as_sf(data_flor, coords = c("Long", "Lat"), crs = 4326, remove = FALSE)
  } else {
    stop("Unsupported data_flor type: must be data.frame, sf, or Spatial object.")
  }
  
  pts_proj <- sf::st_transform(pts_sf, crs_sf)
  pts_in_site <- sf::st_intersection(pts_proj, site_proj_original)
  
  # ============================================================================
  # SECTION 5: OPTIONAL OVERLAP CHECK AND VISUALIZATION
  # ============================================================================
  
  if (check_overlap) {
    overlap_check <- sf::st_intersects(pts_proj, site_proj_processing, sparse = FALSE)
    n_overlap <- sum(overlap_check)
    msg(paste("Number of occurrence points overlapping the processing area:", n_overlap))
    
    if (n_overlap == 0) {
      warning("No points overlap the study area. Check projections or coordinates.")
    } else {
      plot(sf::st_geometry(site_proj_original), col = NA, border = "black", 
           main = "Points vs Site Boundaries", lwd = 2)
      
      if (use_buffer) {
        plot(sf::st_geometry(site_proj_buffered), col = NA, border = "black", 
             add = TRUE, lty = 2, lwd = 1)
      }
      
      if (has_exclusions) {
        plot(sf::st_geometry(excl_proj), col = "lightgray", border = "black", 
             add = TRUE, lty = 3, lwd = 0.5)
      }
      
      points(sf::st_coordinates(pts_proj), col = "red", pch = 20, cex = 0.5)
      
      legend_items <- c("Original site")
      legend_cols <- c("black")
      legend_lty <- c(1)
      legend_lwd <- c(2)
      
      if (use_buffer) {
        legend_items <- c(legend_items, "Expanded site (for analysis)")
        legend_cols <- c(legend_cols, "black")
        legend_lty <- c(legend_lty, 2)
        legend_lwd <- c(legend_lwd, 1)
      }
      
      if (has_exclusions) {
        legend_items <- c(legend_items, "Excluded areas")
        legend_cols <- c(legend_cols, "black")
        legend_lty <- c(legend_lty, 3)
        legend_lwd <- c(legend_lwd, 0.5)
      }
      
      legend("topright", legend = legend_items, col = legend_cols, 
             lty = legend_lty, lwd = legend_lwd, cex = 0.8)
    }
  }
  
  # ============================================================================
  # SECTION 6: OPTIMIZED POINT FILTERING AND BUFFER CREATION
  # ============================================================================
  
  msg("Filtering records and creating buffers (optimized)...")
  msg("  Including external points whose uncertainty buffers reach into study area...")
  
  max_uncertainty <- max(pts_proj$uncertainty)
  site_search_area <- sf::st_buffer(site_proj_processing, dist = max_uncertainty)
  
  potential_pts_idx <- sf::st_intersects(pts_proj, site_search_area, sparse = FALSE)[, 1]
  pts_filtered <- pts_proj[potential_pts_idx, ]
  
  msg(paste("Pre-filtered to", nrow(pts_filtered), "potentially relevant points..."))
  
  buffers_filtered_sf <- sf::st_buffer(pts_filtered, dist = pts_filtered$uncertainty)
  intersects_idx <- sf::st_intersects(buffers_filtered_sf, site_proj_processing, sparse = FALSE)[, 1]
  
  pts_computed <- pts_filtered[intersects_idx, ]
  buffers_computed_sf <- buffers_filtered_sf[intersects_idx, ]
  
  if (nrow(pts_computed) == 0) stop("No occurrence buffers intersect the study area.")
  
  if (has_exclusions) {
    msg("Removing exclusion areas from occurrence buffers...")
    buffers_computed_sf <- sf::st_make_valid(
      sf::st_difference(buffers_computed_sf, excl_proj)
    )
  }
  
  msg(paste("Retained", nrow(pts_computed), "records for computation."))
  
  # ============================================================================
  # SECTION 7: RASTER TEMPLATE CREATION
  # ============================================================================
  
  msg("Creating template raster with guaranteed complete coverage...")
  
  site_bbox <- sf::st_bbox(site_proj_processing)
  bbox_expanded <- site_bbox
  bbox_expanded["xmin"] <- site_bbox["xmin"] - cellsize
  bbox_expanded["xmax"] <- site_bbox["xmax"] + cellsize
  bbox_expanded["ymin"] <- site_bbox["ymin"] - cellsize
  bbox_expanded["ymax"] <- site_bbox["ymax"] + cellsize
  
  r_template <- terra::rast(
    extent = terra::ext(bbox_expanded),
    resolution = cellsize,
    crs = crs_terra
  )
  terra::values(r_template) <- NA
  
  msg(paste0("  Template extent expanded by ", cellsize, "m on all sides for complete coverage."))
  
  # Check if cellsize is appropriate for study area
  site_area <- as.numeric(sf::st_area(site_proj_processing))
  cell_area <- cellsize^2
  n_cells_approx <- site_area / cell_area
  
  if (n_cells_approx < 100) {
    site_dims <- c(
      site_bbox["xmax"] - site_bbox["xmin"],
      site_bbox["ymax"] - site_bbox["ymin"]
    )
    min_dim <- min(site_dims)
    
    warning(paste0("Cell size (", cellsize, "m) is very large for this study area. ",
                   "Expected ~", round(n_cells_approx), " cells covering ", 
                   round(site_area / 1e6, 2), " km^2. ",
                   "Consider reducing cellsize to improve spatial resolution. ",
                   "Suggested maximum: ", round(min_dim / 10), "m (10% of smallest dimension)."))
  }
  
  # ============================================================================
  # SECTION 8: SPECIES RICHNESS MAP CALCULATION
  # ============================================================================
  
  msg("Calculating species richness map (RICH)...")
  
  pts_vect <- terra::vect(pts_computed)
  r_rich <- terra::rasterize(
    pts_vect,
    r_template,
    field = "Taxon",
    fun = function(x) length(unique(x))
  )
  r_rich[is.na(r_rich)] <- 0
  
  # ============================================================================
  # SECTION 9: SPATIO-TEMPORAL SCORE CALCULATION
  # ============================================================================
  
  msg("Calculating spatio-temporal scores...")
  
  v_buff <- terra::vect(buffers_computed_sf)
  cell_data <- terra::extract(r_template, v_buff, cells = TRUE, ID = TRUE)
  counts_per_id <- table(cell_data$ID)
  
  spatial_count <- rep(1, nrow(pts_computed))
  spatial_count[as.integer(names(counts_per_id))] <- as.numeric(counts_per_id)
  
  pts_computed$spatial_score <- 1 / spatial_count
  pts_computed$time_score <- (1 - (tau / 100))^((year_study - pts_computed$year) / 100)
  pts_computed$st_ignorance <- pts_computed$spatial_score * pts_computed$time_score
  
  # ============================================================================
  # SECTION 10: PER-TAXON IGNORANCE RASTERIZATION
  # ============================================================================
  
  taxa_list <- unique(pts_computed$Taxon)
  
  if (use_coverage_weighting) {
    msg("Drafting Map of Relative Floristic Ignorance (processing by taxon)...")
    msg("  Using coverage-weighted rasterization (accurate, slower)...")
    
    pb <- utils::txtProgressBar(min = 0, max = length(taxa_list), style = 3)
    
    raster_list <- suppressWarnings(
      lapply(seq_along(taxa_list), function(i) {
        tname <- taxa_list[i]
        taxon_idx <- pts_computed$Taxon == tname
        pts_taxon <- pts_computed[taxon_idx, ]
        bufs_taxon_sf <- buffers_computed_sf[taxon_idx, ]
        
        tax_r <- r_template
        terra::values(tax_r) <- 0
        
        for (j in 1:nrow(pts_taxon)) {
          single_buf <- sf::st_sf(geometry = sf::st_geometry(bufs_taxon_sf[j, ]))
          coverage_r <- terra::rasterize(terra::vect(single_buf), r_template, cover = TRUE)
          coverage_r[is.na(coverage_r)] <- 0
          weighted_r <- coverage_r * pts_taxon$st_ignorance[j]
          tax_r <- max(tax_r, weighted_r, na.rm = TRUE)
        }
        
        utils::setTxtProgressBar(pb, i)
        return(tax_r)
      })
    )
    close(pb)
    
  } else {
    msg("Drafting Map of Relative Floristic Ignorance (processing by taxon)...")
    msg("  Using binary touch method (fast, no coverage weighting)...")
    
    pb <- utils::txtProgressBar(min = 0, max = length(taxa_list), style = 3)
    
    raster_list <- suppressWarnings(
      lapply(seq_along(taxa_list), function(i) {
        tname <- taxa_list[i]
        taxon_idx <- pts_computed$Taxon == tname
        pts_taxon <- pts_computed[taxon_idx, ]
        bufs_taxon_sf <- buffers_computed_sf[taxon_idx, ]
        
        bufs_taxon_clean <- sf::st_sf(
          st_ignorance = pts_taxon$st_ignorance,
          geometry = sf::st_geometry(bufs_taxon_sf)
        )
        
        tax_r <- terra::rasterize(
          terra::vect(bufs_taxon_clean),
          r_template,
          field = "st_ignorance",
          fun = "max",
          touches = TRUE
        )
        tax_r[is.na(tax_r)] <- 0
        
        utils::setTxtProgressBar(pb, i)
        return(tax_r)
      })
    )
    close(pb)
  }
  
  # ============================================================================
  # SECTION 11: MRFI RASTER FINALIZATION
  # ============================================================================
  
  msg("Finalizing rasters...")
  
  raster_list_valid <- raster_list[sapply(raster_list, function(x) inherits(x, "SpatRaster"))]
  
  if (length(raster_list_valid) > 0) {
    raster_stack <- terra::rast(raster_list_valid)
    raster_sum <- terra::app(raster_stack, fun = sum, na.rm = TRUE)
  } else {
    raster_sum <- r_template
    terra::values(raster_sum) <- 0
  }
  
  rmax <- max(terra::values(raster_sum), na.rm = TRUE)
  if (!is.finite(rmax)) rmax <- 0
  mrfi_r <- rmax - raster_sum
  
  # ============================================================================
  # SECTION 12: MASKING APPLICATION
  # ============================================================================
  
  msg(paste0("Applying mask with method: '", mask_method, "'..."))
  
  if (mask_method == "touches") {
    site_mask_r <- terra::rasterize(site_vect_mask, r_template, field = 1, touches = TRUE)
    site_mask_r[site_mask_r == 0] <- NA
    mrfi_final <- terra::mask(terra::crop(mrfi_r, site_mask_r), site_mask_r)
    rich_final <- terra::mask(terra::crop(r_rich, site_mask_r), site_mask_r)
  } else if (mask_method == "none") {
    msg("  No masking applied - showing entire raster extent.")
    mrfi_final <- mrfi_r
    rich_final <- r_rich
  }
  
  # ============================================================================
  # SECTION 13: EXCLUSION AREA REMOVAL
  # ============================================================================
  
  if (has_exclusions) {
    msg("Removing cells that are >=95% within excluded areas...")
    excl_mask_r <- terra::rasterize(terra::vect(excl_proj), r_template, cover = TRUE)
    
    # Identify cells where exclusion coverage is >=95%
    cells_to_remove <- !is.na(excl_mask_r) & excl_mask_r >= 0.95
    
    # Remove only those cells
    mrfi_final[cells_to_remove] <- NA
    rich_final[cells_to_remove] <- NA
    
    n_removed <- sum(terra::values(cells_to_remove), na.rm = TRUE)
    msg(paste0("  Removed ", n_removed, " cells (>=95% in exclusion areas)"))
  }
  
  # ============================================================================
  # SECTION 14: STATISTICS COMPILATION
  # ============================================================================
  
  msg("Compiling statistics...")
  end_time <- Sys.time()
  
  total_used_records <- nrow(pts_computed)
  excluded_records_count <- total_initial_records - total_used_records
  
  statistics_df <- data.frame(
    Statistic = c(
      "Started", "Finished", "Elapsed time (secs)",
      "Site buffer applied (m)", "Processing & masking boundary",
      "Final masking method", "Coverage weighting", "Exclusion areas applied",
      "CRS (EPSG code)", "Cell size (m)", "100 years % loss ratio (tau)",
      "Total initial records", "Total occurrences within original site",
      "Total occurrences computed (buffers)", "Records excluded from analysis",
      "Occ. uncertainty (median, m)", "Occ. dates (median, year)"
    ),
    Value = c(
      as.character(start_time), as.character(end_time),
      round(as.numeric(end_time - start_time, units = "secs")),
      if (use_buffer) round(buffer_distance, 1) else "None",
      if (use_buffer) paste0("Original site + ", round(buffer_distance, 1), "m expansion") else "Original site only",
      mask_method,
      if (use_coverage_weighting) "Enabled (accurate)" else "Disabled (fast)",
      if (has_exclusions) "Yes" else "No",
      as.character(CRS.new), cellsize, tau,
      total_initial_records, nrow(pts_in_site), nrow(pts_computed),
      excluded_records_count,
      round(median(pts_computed$uncertainty, na.rm = TRUE)),
      round(median(pts_computed$year, na.rm = TRUE))
    )
  )
  
  # ============================================================================
  # SECTION 15: PLOT GENERATION
  # ============================================================================
  
  msg("Generating plots (quantile + continuous versions)...")
  
  excl_plot <- NULL
  if (has_exclusions) {
    excl_plot <- sf::st_intersection(excl_proj, site_proj_processing)
    if (length(excl_plot) == 0 || all(sf::st_is_empty(excl_plot))) {
      excl_plot <- NULL
    }
  }
  
  mrfi_df <- as.data.frame(mrfi_final, xy = TRUE)
  colnames(mrfi_df) <- c("x", "y", "value")
  
  rich_df <- as.data.frame(rich_final, xy = TRUE)
  colnames(rich_df) <- c("x", "y", "value")
  
  mrfi_max_val <- terra::global(mrfi_final, "max", na.rm = TRUE)$max
  rich_max_val <- terra::global(rich_final, "max", na.rm = TRUE)$max
  
  # Calculate quantile breaks
  msg("  Calculating quantile breaks (8 quantiles)...")
  
  mrfi_values <- terra::values(mrfi_final, na.rm = TRUE)
  
  if (length(unique(mrfi_values)) >= 8) {
    mrfi_breaks_quant <- stats::quantile(mrfi_values, probs = seq(0, 1, length.out = 9), na.rm = TRUE)
    mrfi_breaks_quant <- unique(mrfi_breaks_quant)
  } else {
    n_unique <- length(unique(mrfi_values))
    mrfi_breaks_quant <- stats::quantile(mrfi_values, probs = seq(0, 1, length.out = n_unique + 1), na.rm = TRUE)
    mrfi_breaks_quant <- unique(mrfi_breaks_quant)
  }
  
  rich_values_nonzero <- terra::values(rich_final, na.rm = TRUE)
  rich_values_nonzero <- rich_values_nonzero[rich_values_nonzero > 0]
  
  if (length(rich_values_nonzero) > 0) {
    rich_breaks_quant <- stats::quantile(rich_values_nonzero, probs = seq(0, 1, length.out = 8), na.rm = TRUE)
    rich_breaks_quant <- unique(rich_breaks_quant)
    if (length(rich_breaks_quant) < 8) {
      rich_breaks_quant <- seq(min(rich_values_nonzero), max(rich_values_nonzero), length.out = 8)
    }
    rich_breaks_quant <- c(0, rich_breaks_quant)
  } else {
    rich_breaks_quant <- seq(0, max(terra::values(rich_final, na.rm = TRUE)), length.out = 9)
  }
  
  mrfi_breaks_cont <- seq(0, mrfi_max_val, length.out = 9)
  rich_breaks_cont <- seq(0, rich_max_val, length.out = 9)
  
  # Helper function to add boundaries
  add_boundaries <- function(p) {
    p <- p + ggplot2::geom_sf(data = site_proj_original, fill = NA, color = "black", 
                              linewidth = 1, inherit.aes = FALSE)
    if (use_buffer && !is.null(site_proj_buffered)) {
      p <- p + ggplot2::geom_sf(data = site_proj_buffered, fill = NA, color = "black", 
                                linewidth = 0.8, linetype = "dashed", inherit.aes = FALSE)
    }
    if (has_exclusions && !is.null(excl_plot)) {
      p <- p + ggplot2::geom_sf(data = excl_plot, fill = "gray", color = "black", 
                                linewidth = 0.5, linetype = "dotted", inherit.aes = FALSE, alpha = 0.3)
    }
    return(p)
  }
  
  # MRFI - Quantile
  p1_quant <- ggplot2::ggplot(mrfi_df) +
    ggplot2::coord_equal() +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "right", legend.direction = 'vertical',
                   legend.key.width = grid::unit(0.6, "cm"),
                   plot.title = ggplot2::element_text(size = 11, face = "bold")) +
    ggplot2::geom_tile(ggplot2::aes(x = .data$x, y = .data$y, fill = .data$value),
                       alpha = 0.8, colour = "black", linewidth = 0.1) +
    ggplot2::scale_fill_gradientn(
      colors = rev(RColorBrewer::brewer.pal(9, "Spectral")),
      values = scales::rescale(mrfi_breaks_quant),
      breaks = mrfi_breaks_quant, labels = round(mrfi_breaks_quant, 1),
      limits = c(min(mrfi_breaks_quant), max(mrfi_breaks_quant)),
      na.value = "transparent",
      guide = ggplot2::guide_legend(title = "IRFI", keyheight = grid::unit(1.2, "lines"),
                                    keywidth = grid::unit(1.2, "lines"),
                                    label.theme = ggplot2::element_text(size = 9),
                                    title.theme = ggplot2::element_text(size = 10, face = "bold"))
    ) +
    ggplot2::ggtitle("MRFI - Quantile Scale (Octiles)") +
    ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude")
  p1_quant <- add_boundaries(p1_quant)
  
  # Richness - Quantile
  p2_quant <- ggplot2::ggplot(rich_df) +
    ggplot2::coord_equal() +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "right", legend.direction = 'vertical',
                   legend.key.width = grid::unit(0.6, "cm"),
                   plot.title = ggplot2::element_text(size = 11, face = "bold")) +
    ggplot2::geom_tile(ggplot2::aes(x = .data$x, y = .data$y, fill = .data$value),
                       alpha = 0.8, colour = "black", linewidth = 0.1) +
    ggplot2::scale_fill_gradientn(
      colors = RColorBrewer::brewer.pal(9, "Spectral"),
      values = scales::rescale(rich_breaks_quant),
      breaks = rich_breaks_quant, labels = round(rich_breaks_quant, 1),
      limits = c(min(rich_breaks_quant), max(rich_breaks_quant)),
      na.value = "transparent",
      guide = ggplot2::guide_legend(title = "N taxa", keyheight = grid::unit(1.2, "lines"),
                                    keywidth = grid::unit(1.2, "lines"),
                                    label.theme = ggplot2::element_text(size = 9),
                                    title.theme = ggplot2::element_text(size = 10, face = "bold"))
    ) +
    ggplot2::ggtitle("Species Richness - Quantile Scale (Octiles)") +
    ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude")
  p2_quant <- add_boundaries(p2_quant)
  
  # MRFI - Continuous
  p1_cont <- ggplot2::ggplot(mrfi_df) +
    ggplot2::coord_equal() +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "right", legend.direction = 'vertical',
                   legend.key.width = grid::unit(0.6, "cm"),
                   plot.title = ggplot2::element_text(size = 11, face = "bold")) +
    ggplot2::geom_tile(ggplot2::aes(x = .data$x, y = .data$y, fill = .data$value),
                       alpha = 0.8, colour = "black", linewidth = 0.1) +
    ggplot2::scale_fill_gradientn(
      colors = rev(RColorBrewer::brewer.pal(9, "Spectral")),
      limits = c(0, mrfi_max_val),
      breaks = mrfi_breaks_cont, labels = round(mrfi_breaks_cont, 1),
      na.value = "transparent",
      guide = ggplot2::guide_legend(title = "IRFI", keyheight = grid::unit(1.2, "lines"),
                                    keywidth = grid::unit(1.2, "lines"),
                                    label.theme = ggplot2::element_text(size = 9),
                                    title.theme = ggplot2::element_text(size = 10, face = "bold"))
    ) +
    ggplot2::ggtitle("MRFI - Continuous Scale") +
    ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude")
  p1_cont <- add_boundaries(p1_cont)
  
  # Richness - Continuous
  p2_cont <- ggplot2::ggplot(rich_df) +
    ggplot2::coord_equal() +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "right", legend.direction = 'vertical',
                   legend.key.width = grid::unit(0.6, "cm"),
                   plot.title = ggplot2::element_text(size = 11, face = "bold")) +
    ggplot2::geom_tile(ggplot2::aes(x = .data$x, y = .data$y, fill = .data$value),
                       alpha = 0.8, colour = "black", linewidth = 0.1) +
    ggplot2::scale_fill_gradientn(
      colors = RColorBrewer::brewer.pal(9, "Spectral"),
      limits = c(0, rich_max_val),
      breaks = rich_breaks_cont, labels = round(rich_breaks_cont, 1),
      na.value = "transparent",
      guide = ggplot2::guide_legend(title = "N taxa", keyheight = grid::unit(1.2, "lines"),
                                    keywidth = grid::unit(1.2, "lines"),
                                    label.theme = ggplot2::element_text(size = 9),
                                    title.theme = ggplot2::element_text(size = 10, face = "bold"))
    ) +
    ggplot2::ggtitle("Species Richness - Continuous Scale") +
    ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude")
  p2_cont <- add_boundaries(p2_cont)
  
  # Histograms
  
  # 1. Temporal distribution histogram
  p3 <- ggplot2::ggplot(pts_computed) +
    ggplot2::aes(x = .data$year) +
    ggplot2::geom_histogram(alpha = 0.6, fill = "#FF6666",
                            binwidth = diff(range(pts_computed$year)) / 30,
                            color = "black", linewidth = 0.3) +
    ggplot2::coord_cartesian(xlim = c(min(pts_computed$year), year_study)) +
    ggplot2::ggtitle("Temporal Distribution - Occurrence Dates") +
    ggplot2::xlab("Year of occurrence") + 
    ggplot2::ylab("Number of records") +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11, face = "bold"))
  
  # 3. Spatial uncertainty distribution histogram
  p4 <- ggplot2::ggplot(pts_computed) +
    ggplot2::aes(x = .data$uncertainty) +
    ggplot2::geom_histogram(alpha = 0.6, fill = "#66B2FF",
                            binwidth = diff(range(pts_computed$uncertainty)) / 30,
                            color = "black", linewidth = 0.3) +
    ggplot2::coord_cartesian(xlim = c(min(pts_computed$uncertainty), 
                                      max(pts_computed$uncertainty))) +
    ggplot2::ggtitle("Spatial Distribution - Occurrence Uncertainty") +
    ggplot2::xlab("Uncertainty (m)") + 
    ggplot2::ylab("Number of records") +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11, face = "bold"))
  
  # ============================================================================
  # SECTION 16: FILE OUTPUT
  # ============================================================================
  
  msg("Saving output files...")
  
  pdf_path <- file.path(output_dir, paste0(output_prefix, "_output.pdf"))
  tif_path <- file.path(output_dir, paste0(output_prefix, "_map.tif"))
  csv_path <- file.path(output_dir, paste0(output_prefix, "_taxa.csv"))
  
  # Detect study area aspect ratio for optimal layout
  site_bbox <- sf::st_bbox(site_proj_mask)
  bbox_width <- site_bbox$xmax - site_bbox$xmin
  bbox_height <- site_bbox$ymax - site_bbox$ymin
  aspect_ratio <- bbox_width / bbox_height
  
  # Determine layout: wide areas (E-W) use vertical stacking, tall/square use horizontal
  if (aspect_ratio > 1.5) {
    layout_ncol <- 1
    layout_type <- "vertical"
    msg("  Study area is wide (E-W oriented): using vertical layout for better map visibility")
  } else {
    layout_ncol <- 2
    layout_type <- "horizontal"
    msg("  Study area is tall/square: using side-by-side layout")
  }
  
  grDevices::pdf(pdf_path, width = 11.69, height = 8.27, onefile = TRUE)
  
  gridExtra::grid.arrange(p1_quant, p2_quant, ncol = layout_ncol,
                          top = grid::textGrob("Page 1: Quantile Comparison (Primary Analysis)", 
                                               gp = grid::gpar(fontsize = 14, fontface = "bold")))
  
  gridExtra::grid.arrange(p1_cont, p2_cont, ncol = layout_ncol,
                          top = grid::textGrob("Page 2: Continuous Scale (Raw Data View)", 
                                               gp = grid::gpar(fontsize = 14, fontface = "bold")))
  
  # Page 3 with diagnostic histograms
  gridExtra::grid.arrange(
    p3, p4,
    ncol = 2,
    top = grid::textGrob("Page 3: Data Diagnostics", 
                         gp = grid::gpar(fontsize = 14, fontface = "bold"))
  )
  
  grid::grid.draw(
    gridExtra::grid.arrange(
      top = grid::textGrob("Page 4: Summary Statistics", 
                           gp = grid::gpar(fontsize = 14, fontface = "bold")),
      gridExtra::tableGrob(statistics_df)
    )
  )
  
  grDevices::dev.off()
  
  terra::writeRaster(mrfi_final, filename = tif_path, overwrite = TRUE)
  utils::write.csv(taxa_list, row.names = FALSE, csv_path)
  
  msg(paste0("Done! Files saved to: ", output_dir))
  msg(paste0("PDF layout: ", layout_type, " (aspect ratio: ", round(aspect_ratio, 2), ")"))
  msg("PDF structure: Page 1 (Quantile), Page 2 (Continuous), Page 3 (Diagnostics), Page 4 (Statistics)")
  
  # ============================================================================
  # SECTION 17: CONSOLE DISPLAY
  # ============================================================================
  
  print(p1_quant)
  print(p2_quant)
  
  # ============================================================================
  # SECTION 18: RETURN RESULTS
  # ============================================================================
  
  return(list(
    MRFI = mrfi_final,
    RICH = rich_final,
    Uncertainties = data.frame(
      uncertainty = pts_computed$uncertainty,
      year = pts_computed$year,
      Taxon = pts_computed$Taxon
    ),
    Statistics = statistics_df
  ))
}