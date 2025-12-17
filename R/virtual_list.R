#' @title Virtual Floristic List
#'
#' @description
#' Computes a Virtual Floristic List (VFL): taxa potentially occurring within
#' a study site, with probabilities based on spatial uncertainty and temporal
#' decay using inclusion-exclusion principle.
#'
#' @param data_flor Data frame with columns: 'Taxon', 'Long', 'Lat',
#'   'uncertainty' (radius in meters), 'year'.
#' @param site sf polygon or SpatialPolygonsDataFrame of study area.
#' @param year_study Numeric year of analysis (default = current year).
#' @param excl_areas Optional sf polygon(s) of unsuitable areas to exclude.
#' @param CRS.new Numeric EPSG code for projected CRS (default = 3035).
#' @param tau Percent taxa loss per 100 years (0 ≤ tau < 100).
#' @param upperlimit Maximum number of records per taxon used in probability
#'   calculation (default = 20). Prevents computational explosion for well-sampled
#'   taxa. Higher values increase accuracy but dramatically slow computation:
#'   \itemize{
#'     \item{10}{ Very fast, good accuracy - for exploratory analysis}
#'     \item{20}{ Fast, very good accuracy - recommended default}
#'     \item{30}{ Slow, excellent accuracy - for publication-quality results}
#'   }
#' @param min_probability Minimum probability threshold (%) for inclusion in VFL
#'   (default = 0). Set to 5-10 to filter unlikely taxa.
#' @param verbose Logical; print progress messages (default = TRUE).
#' @param check_overlap Logical; plot spatial overlap diagnostic (default = TRUE).
#' @param output_dir Directory for output files (default = working directory).
#' @param output_prefix Filename prefix (default = "VFL").
#'
#' @return List with:
#' \itemize{
#'   \item{\code{VFL}}{Data frame: Taxon, probability, records, max, min.}
#'   \item{\code{Statistics}}{Metadata table.}
#'   \item{\code{Plots}}{Named list of 2 ggplot objects (histograms).}
#'   \item{\code{spatial_data}}{sf objects for further analysis.}
#' }
#'
#' @details
#' The function uses the inclusion-exclusion principle to aggregate probabilities
#' from multiple records of the same taxon. For a taxon with n records, the
#' computation requires 2^n combinations. The \code{upperlimit} parameter prevents
#' exponential explosion: with n > upperlimit, only the upperlimit records with
#' highest probabilities are used.
#'
#' PDF Output Structure (2 pages, A4 landscape):
#' \itemize{
#'   \item{Page 1}{ Distribution analysis (probability + temporal histograms)}
#'   \item{Page 2}{ Summary statistics table}
#' }
#'
#' @importFrom sf st_as_sf st_make_valid st_crs st_transform st_buffer st_intersection st_union st_area st_geometry st_difference st_coordinates st_intersects
#' @importFrom ggplot2 ggplot aes geom_histogram geom_point theme_classic xlab ylab ggtitle geom_abline element_text
#' @importFrom utils write.csv txtProgressBar setTxtProgressBar combn
#' @importFrom grDevices pdf dev.off rgb
#' @importFrom grid grid.draw textGrob gpar
#' @importFrom gridExtra grid.arrange tableGrob
#' @export
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(floratus)
#' data(park)
#'
#' # Basic usage (creates VFL_output.pdf and VFL_probabilities.csv)
#' vfl <- virtual_list(
#'   data_flor = floratus,
#'   site = park,
#'   tau = 30
#' )
#'
#' # Custom naming (creates Park_output.pdf and Park_probabilities.csv)
#' vfl <- virtual_list(
#'   data_flor = floratus,
#'   site = park,
#'   tau = 30,
#'   output_prefix = "Park"
#' )
#'
#' # Filter unlikely taxa (>5% probability only)
#' vfl_filtered <- virtual_list(
#'   data_flor = floratus,
#'   site = park,
#'   tau = 30,
#'   min_probability = 5
#' )
#'
#' # High accuracy mode
#' vfl_accurate <- virtual_list(
#'   data_flor = floratus,
#'   site = park,
#'   tau = 30,
#'   upperlimit = 30
#' )
#' }

virtual_list <- function(data_flor, site, year_study = NULL, excl_areas = NULL, 
                         CRS.new = 3035, tau, upperlimit = 20, 
                         min_probability = 0, verbose = TRUE,
                         check_overlap = TRUE,
                         output_dir = output_dir = getwd(),
                         output_prefix = "VFL") {
  
  # ============================================================================
  # HELPER FUNCTIONS
  # ============================================================================
  
  msg <- function(...) if (verbose) message(...)
  start_time <- Sys.time()
  
  # Compute taxon probability using inclusion-exclusion principle
  compute_taxon_prob <- function(probs, max_records) {
    n_records <- length(probs)
    
    if (n_records == 0) {
      return(c(prob = 0, n_records = 0, max_prob = 0, min_prob = 0))
    }
    
    # Limit number of records to prevent combinatorial explosion
    # 2^n combinations: 2^20 = 1M (manageable), 2^30 = 1B (slow), 2^40 = 1T (crash)
    if (n_records > max_records) {
      probs <- sort(probs, decreasing = TRUE)[1:max_records]
      n_records <- max_records
    }
    
    if (n_records == 1) {
      return(c(prob = probs[1], n_records = 1, max_prob = probs[1], min_prob = probs[1]))
    }
    
    # Inclusion-exclusion principle: P(A₁∪...∪Aₙ) = ΣP(Aᵢ) - ΣP(Aᵢ∩Aⱼ) + ...
    total_prob <- 0
    
    for (k in 1:n_records) {
      combs <- utils::combn(n_records, k)
      sign_factor <- (-1)^(k + 1)
      k_sum <- sum(apply(combs, 2, function(idx) prod(probs[idx])))
      total_prob <- total_prob + sign_factor * k_sum
    }
    
    return(c(
      prob = total_prob,
      n_records = length(probs),
      max_prob = max(probs),
      min_prob = min(probs)
    ))
  }
  
  # ============================================================================
  # SECTION 1: INPUT VALIDATION
  # ============================================================================
  
  msg("Validating inputs...")
  
  if (is.null(year_study)) year_study <- as.numeric(format(Sys.Date(), "%Y"))
  
  if (tau < 0 || tau >= 100) stop("tau must satisfy: 0 <= tau < 100")
  if (upperlimit < 1) stop("upperlimit must be at least 1")
  if (upperlimit > 35) {
    warning("Large upperlimit (>35) may cause very slow computation or memory issues. ",
            "Recommended range: 10-30. Consider using 20-25 for most analyses.")
  }
  
  if (!is.numeric(min_probability) || min_probability < 0 || min_probability >= 100) {
    stop("min_probability must be between 0 and 100 (percentage)")
  }
  
  req_cols <- c("Taxon", "Long", "Lat", "uncertainty", "year")
  missing_cols <- setdiff(req_cols, names(data_flor))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if (max(data_flor$year, na.rm = TRUE) > year_study) {
    warning("Some occurrence dates are more recent than year_study")
  }
  
  total_initial_records <- nrow(data_flor)
  total_initial_taxa <- length(unique(data_flor$Taxon))
  
  msg("Inputs validated.")
  
  # ============================================================================
  # SECTION 2: CRS SETUP AND REPROJECTION
  # ============================================================================
  
  msg(paste("Reprojecting to EPSG:", CRS.new))
  
  crs_sf <- sf::st_crs(CRS.new)
  
  if (inherits(site, "Spatial")) site <- sf::st_as_sf(site)
  if (is.na(sf::st_crs(site))) {
    msg("Site has no CRS; assuming EPSG:4326")
    sf::st_crs(site) <- 4326
  }
  site_proj <- sf::st_transform(sf::st_make_valid(site), crs_sf)
  
  excl_proj <- NULL
  has_exclusions <- !is.null(excl_areas)
  
  if (has_exclusions) {
    if (inherits(excl_areas, "Spatial")) excl_areas <- sf::st_as_sf(excl_areas)
    if (is.na(sf::st_crs(excl_areas))) {
      msg("Exclusion areas have no CRS; assuming EPSG:4326")
      sf::st_crs(excl_areas) <- 4326
    }
    excl_proj <- sf::st_union(sf::st_transform(sf::st_make_valid(excl_areas), crs_sf))
    msg("Exclusion areas will be removed from buffers")
  }
  
  if (inherits(data_flor, "sf")) {
    pts_sf <- data_flor
  } else if (inherits(data_flor, "Spatial")) {
    pts_sf <- sf::st_as_sf(data_flor)
  } else {
    pts_sf <- sf::st_as_sf(data_flor, coords = c("Long", "Lat"), crs = 4326, remove = FALSE)
  }
  
  pts_proj <- sf::st_transform(pts_sf, crs_sf)
  
  # ============================================================================
  # SECTION 3: CREATE BUFFERS
  # ============================================================================
  
  msg("Creating uncertainty buffers...")
  
  pts_proj$record_id <- seq_len(nrow(pts_proj))
  buffers_sf <- sf::st_buffer(pts_proj, dist = pts_proj$uncertainty)
  buffers_sf$area_buffer <- as.numeric(sf::st_area(buffers_sf))
  
  if (has_exclusions) {
    msg("Removing exclusion areas from buffers...")
    buffers_sf <- sf::st_make_valid(sf::st_difference(buffers_sf, excl_proj))
    buffers_sf$area_buffer <- as.numeric(sf::st_area(buffers_sf))
  }
  
  # ============================================================================
  # SECTION 4: OPTIONAL OVERLAP CHECK
  # ============================================================================
  
  if (check_overlap) {
    msg("Checking spatial overlap between occurrences and study area...")
    
    intersects_test <- sf::st_intersects(buffers_sf, site_proj, sparse = FALSE)[, 1]
    n_intersecting <- sum(intersects_test)
    
    msg(paste("  ", n_intersecting, "of", nrow(buffers_sf), "buffers intersect the study area"))
    
    if (n_intersecting == 0) {
      warning("No buffers intersect the study area! Check CRS and coordinates.")
    }
    
    plot(sf::st_geometry(site_proj), col = NA, border = "black", 
         main = "Spatial Overlap Check", lwd = 2)
    
    if (has_exclusions) {
      plot(sf::st_geometry(excl_proj), col = "lightgray", border = "gray50",
           add = TRUE, lty = 3)
    }
    
    points(sf::st_coordinates(pts_proj), col = "red", pch = 20, cex = 0.5)
    
    legend_items <- "Study site"
    legend_cols <- "black"
    legend_pch <- NA
    legend_lty <- 1
    
    if (has_exclusions) {
      legend_items <- c(legend_items, "Exclusion areas")
      legend_cols <- c(legend_cols, "gray50")
      legend_pch <- c(NA, NA)
      legend_lty <- c(1, 3)
    }
    
    legend_items <- c(legend_items, "Occurrence points")
    legend_cols <- c(legend_cols, "red")
    legend_pch <- c(legend_pch, 20)
    legend_lty <- c(legend_lty, NA)
    
    legend("topright", legend = legend_items, col = legend_cols,
           pch = legend_pch, lty = legend_lty, cex = 0.8)
  }
  
  # ============================================================================
  # SECTION 5: INTERSECTION WITH STUDY AREA
  # ============================================================================
  
  msg("Computing buffer intersections with study area...")
  
  intersects_idx <- sf::st_intersects(buffers_sf, site_proj, sparse = FALSE)[, 1]
  
  if (sum(intersects_idx) == 0) {
    stop("No occurrence buffers intersect the study area. Check CRS and coordinates.")
  }
  
  buffers_intersecting <- buffers_sf[intersects_idx, ]
  msg(paste("Retained", nrow(buffers_intersecting), "buffers that intersect site"))
  
  inter_sf <- sf::st_intersection(buffers_intersecting, site_proj)
  inter_sf$area_intersection <- as.numeric(sf::st_area(inter_sf))
  
  # ============================================================================
  # SECTION 6: CALCULATE PROBABILITIES
  # ============================================================================
  
  msg("Calculating spatial and temporal probabilities...")
  
  inter_sf$p_spatial <- inter_sf$area_intersection / inter_sf$area_buffer
  inter_sf$p_temporal <- (1 - (tau / 100))^((year_study - inter_sf$year) / 100)
  inter_sf$p_st <- inter_sf$p_spatial * inter_sf$p_temporal
  
  if (any(inter_sf$p_st > 1.0, na.rm = TRUE)) {
    warning("Some probabilities exceed 1.0; capping at 1.0")
    inter_sf$p_st <- pmin(inter_sf$p_st, 1.0)
  }
  if (any(inter_sf$p_st < 0, na.rm = TRUE)) {
    warning("Some probabilities are negative; setting to 0")
    inter_sf$p_st <- pmax(inter_sf$p_st, 0)
  }
  
  # ============================================================================
  # SECTION 7: AGGREGATE BY TAXON (INCLUSION-EXCLUSION PRINCIPLE)
  # ============================================================================
  
  msg("Aggregating probabilities by taxon using inclusion-exclusion...")
  
  taxa_list <- unique(inter_sf$Taxon)
  n_taxa <- length(taxa_list)
  
  pb <- utils::txtProgressBar(min = 0, max = n_taxa, style = 3)
  
  results <- lapply(seq_along(taxa_list), function(i) {
    taxon <- taxa_list[i]
    taxon_probs <- inter_sf$p_st[inter_sf$Taxon == taxon]
    result <- compute_taxon_prob(taxon_probs, upperlimit)
    utils::setTxtProgressBar(pb, i)
    result
  })
  
  close(pb)
  
  # ============================================================================
  # SECTION 8: BUILD VFL DATA FRAME
  # ============================================================================
  
  msg("Building Virtual Floristic List...")
  
  VFL <- data.frame(
    Taxon = taxa_list,
    Estimated_Spatiotemporal_probability = sapply(results, `[`, 1),
    Number_of_records = sapply(results, `[`, 2),
    Max_probability = sapply(results, `[`, 3),
    Min_probability = sapply(results, `[`, 4),
    stringsAsFactors = FALSE
  )
  
  VFL <- VFL[VFL$Estimated_Spatiotemporal_probability > (min_probability / 100), ]
  
  if (nrow(VFL) == 0) {
    warning(paste("No taxa meet the minimum probability threshold of", min_probability, "%"))
  }
  
  VFL$Estimated_Spatiotemporal_probability <- 
    round(VFL$Estimated_Spatiotemporal_probability * 100, 1)
  VFL$Max_probability <- round(VFL$Max_probability * 100, 1)
  VFL$Min_probability <- round(VFL$Min_probability * 100, 1)
  
  VFL <- VFL[order(-VFL$Estimated_Spatiotemporal_probability, VFL$Taxon), ]
  rownames(VFL) <- NULL
  
  # ============================================================================
  # SECTION 9: COMPILE STATISTICS
  # ============================================================================
  
  msg("Compiling statistics...")
  end_time <- Sys.time()
  
  statistics_df <- data.frame(
    Statistic = c(
      "Started", "Finished", "Elapsed time (secs)",
      "CRS (EPSG code)", "Year of study", "Tau (%)",
      "Upperlimit", "Min probability filter (%)",
      "Exclusion areas applied",
      "Total initial records", "Total initial taxa",
      "Records intersecting site", "Taxa in VFL (after filtering)",
      "Occ. uncertainty (median, m)", "Occ. dates (median, year)"
    ),
    Value = c(
      as.character(start_time), as.character(end_time),
      round(as.numeric(difftime(end_time, start_time, units = "secs")), 1),
      as.character(CRS.new), year_study, tau,
      upperlimit, min_probability,
      if (has_exclusions) "Yes" else "No",
      total_initial_records, total_initial_taxa,
      nrow(inter_sf), nrow(VFL),
      round(median(inter_sf$uncertainty, na.rm = TRUE)),
      round(median(inter_sf$year, na.rm = TRUE))
    ),
    stringsAsFactors = FALSE
  )
  
  # ============================================================================
  # SECTION 10: GENERATE DIAGNOSTIC PLOTS
  # ============================================================================
  
  msg("Generating diagnostic plots...")
  
  # 1. Probability distribution histogram
  p_prob_hist <- ggplot2::ggplot(VFL, ggplot2::aes(x = Estimated_Spatiotemporal_probability)) +
    ggplot2::geom_histogram(fill = "#66CC66", alpha = 0.7, bins = 30, 
                            color = "black", linewidth = 0.3) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold")) +
    ggplot2::xlab("Spatio-temporal probability (%)") +
    ggplot2::ylab("Number of taxa") +
    ggplot2::ggtitle("Probability Distribution")
  
  # 2. Temporal distribution histogram
  p_temporal <- ggplot2::ggplot(as.data.frame(inter_sf), ggplot2::aes(x = year)) +
    ggplot2::geom_histogram(fill = "#FF6666", alpha = 0.7, bins = 30, 
                            color = "black", linewidth = 0.3) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold")) +
    ggplot2::xlab("Year of occurrence") +
    ggplot2::ylab("Number of records") +
    ggplot2::ggtitle("Temporal Distribution")
  
  # 3. Spatial uncertainty distribution histogram
  p_spatial <- ggplot2::ggplot(as.data.frame(inter_sf), ggplot2::aes(x = uncertainty)) +
    ggplot2::geom_histogram(fill = "#66B2FF", alpha = 0.7, bins = 30,
                            color = "black", linewidth = 0.3) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold")) +
    ggplot2::xlab("Uncertainty (m)") +
    ggplot2::ylab("Number of records") +
    ggplot2::ggtitle("Spatial Uncertainty Distribution")
  
  plots_list <- list(
    prob_histogram = p_prob_hist,
    temporal_distribution = p_temporal,
    spatial_distribution = p_spatial
  )
  
  # ============================================================================
  # SECTION 11: SAVE OUTPUTS
  # ============================================================================
  
  msg("Saving output files...")
  
  # Define output paths FIRST (before using them!)
  pdf_path <- file.path(output_dir, paste0(output_prefix, "_output.pdf"))
  csv_path <- file.path(output_dir, paste0(output_prefix, "_probabilities.csv"))
  
  # Save CSV
  utils::write.csv(VFL, csv_path, row.names = FALSE)
  
  # Create 2-page PDF with 3 histograms
  grDevices::pdf(pdf_path, width = 11.69, height = 8.27, onefile = TRUE)
  
  # Page 1: Distribution Analysis (3 plots: probability full width, temporal + spatial below)
  gridExtra::grid.arrange(
    p_prob_hist,
    gridExtra::arrangeGrob(p_temporal, p_spatial, ncol = 2),
    ncol = 1,
    heights = c(1, 1),
    top = grid::textGrob("Page 1: Distribution Analysis", 
                         gp = grid::gpar(fontsize = 14, fontface = "bold"))
  )
  
  # Page 2: Summary Statistics
  grid::grid.draw(
    gridExtra::grid.arrange(
      top = grid::textGrob("Page 2: Summary Statistics", 
                           gp = grid::gpar(fontsize = 14, fontface = "bold")),
      gridExtra::tableGrob(statistics_df)
    )
  )
  
  grDevices::dev.off()
  
  msg(paste0("Done! Files saved to: ", output_dir))
  msg("PDF structure: Page 1 (Distribution Analysis with 3 histograms), Page 2 (Statistics)")
  
  # ============================================================================
  # SECTION 12: CONSOLE DISPLAY
  # ============================================================================
  
  if (verbose) {
    print(p_prob_hist)
  }
  
  # ============================================================================
  # SECTION 13: RETURN RESULTS
  # ============================================================================
  
  return(list(
    VFL = VFL,
    Statistics = statistics_df,
    Plots = plots_list,
    spatial_data = list(
      site = site_proj,
      buffers = buffers_intersecting,
      intersections = inter_sf
    )
  ))
}