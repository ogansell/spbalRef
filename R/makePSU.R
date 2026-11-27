#' Create Two-Stage Spatial Sampling Design with Primary Sampling Units
#'
#' @description
#' Creates a two-stage spatial sampling design by dividing a study area into
#' approximately equal-area primary sampling units (PSUs) using balanced clustering
#' and Voronoi tessellation, then selecting secondary sample points within each PSU
#' using spatially balanced sampling (GRTS). The function can optimize PSU boundaries
#' to minimize area variation while maintaining spatial balance, or use a fast
#' single-pass approach for model-based inference.
#'
#' @param sf_poly An \code{sf} object representing the sample frame (study area polygon).
#'   Can be a single polygon or multiple polygons that will be unioned.
#' @param psu_area Numeric. Target area for each PSU in hectares.
#' @param cluster_population Integer. Number of random points to generate for the
#'   clustering (default = 200000 for optimization, 50000 for quick_mode). Higher
#'   values may improve results but increase computation time.
#' @param optimization_start Integer. Starting number of points for the optimization
#'   loop (default = 10000). Must be less than \code{cluster_population - 10000}.
#'   Ignored if \code{quick_mode = TRUE}.
#' @param optimization_step Integer. Step size between iterations in the
#'   optimization loop (default = 3000). Ignored if \code{quick_mode = TRUE}.
#' @param optimization_metric Character. Metric to optimize PSU areas. One of
#'   \code{"cv"} (coefficient of variation, default), \code{"range"}, or \code{"sd"}
#'   (standard deviation). Ignored if \code{quick_mode = TRUE}.
#' @param n_samples Integer. Number of secondary sample points to select within
#'   each PSU.
#' @param min_radius Numeric or NULL. Minimum distance (in meters) between
#'   secondary sample points, passed to \code{\link[spbal]{BAS}} as
#'   \code{minRadius}. Default \code{NULL} (no minimum distance constraint).
#' @param island Character. Which New Zealand island master sample to use.
#'   One of \code{"North"} or \code{"South"}.
#' @param min.dist Numeric. Minimum distance between PSU centroids for GRTS
#'   sampling (default = 0).
#' @param base_n Integer. Number of base PSUs to select. If \code{NULL} (default),
#'   all PSUs are selected.
#' @param over_n Integer. Number of oversample PSUs to select. If \code{NULL}
#'   (default), no oversample is created.
#' @param legacy_sites An optional \code{sf} object containing legacy monitoring
#'   sites to prioritize in PSU selection.
#' @param legacy_id_col Character. Name of the column in \code{legacy_sites}
#'   containing site identifiers (default = "PlotID").
#' @param grts_seltype Character. GRTS selection type passed to
#'   \code{\link[spsurvey]{grts}}. Default is \code{"equal"}.
#' @param grts_caty_var Character. Optional stratification variable name for
#'   categorical GRTS sampling.
#' @param grts_aux_var Character. Optional auxiliary variable name for GRTS sampling.
#' @param seed Integer. Random seed for reproducibility (default = 25945).
#' @param plot Logical. Whether to create diagnostic plots (default = TRUE).
#' @param verbose Logical. Whether to print progress messages (default = TRUE).
#' @param quick_mode Logical. If \code{TRUE}, uses single-pass balanced clustering
#'   without optimization (default = FALSE). Recommended for model-based inference
#'   where PSU size uniformity is less critical. Reduces runtime by ~80-90\%.
#'
#' @return
#' An object of class \code{"psu"} containing:
#' \describe{
#'   \item{primary_sample_units}{An \code{sf} object with PSU polygons and sample order}
#'   \item{secondary_sample_units}{An \code{sf} object with sample points within PSUs}
#'   \item{optimization_summary}{A data frame summarizing optimization metrics (NULL in quick_mode)}
#'   \item{optimization_data}{Full optimization data (NULL in quick_mode)}
#'   \item{optimal_n_points}{The number of points used}
#'   \item{pielou_evenness}{Pielou evenness index for spatial balance}
#'   \item{psu_area_stats}{List of PSU area statistics (mean, sd, cv, range, etc.)}
#'   \item{optimization_metric}{The metric used for optimization (NA in quick_mode)}
#'   \item{settings}{List of function settings used}
#' }
#'
#' @details
#' The function implements a multi-step process:
#'
#' \strong{Standard mode (quick_mode = FALSE):}
#' \enumerate{
#'   \item Generates random points across the study area
#'   \item Iterates through different point densities, using balanced clustering
#'         at each iteration to create cluster centroids
#'   \item Creates Voronoi polygons around centroids and clips to study area
#'   \item Calculates area statistics and spatial balance metrics
#'   \item Selects the iteration that minimizes the chosen optimization metric
#'   \item Applies GRTS to PSU centroids to determine sampling order
#'   \item Generates secondary sample points within selected PSUs using BAS
#' }
#'
#' \strong{Quick mode (quick_mode = TRUE):}
#' \enumerate{
#'   \item Generates random points (fewer than standard mode)
#'   \item Single-pass balanced clustering to create PSU centroids
#'   \item Creates Voronoi polygons and clips to study area
#'   \item Applies GRTS to PSU centroids
#'   \item Generates secondary sample points using BAS
#' }
#'
#' The balanced clustering algorithm (\code{\link[anticlust]{balanced_clustering}})
#' tends to produce more equal-sized PSUs than standard k-means clustering.
#'
#' If \code{legacy_sites} are provided, PSUs containing these sites will be
#' prioritized in the GRTS selection using the legacy variable option.
#'
#' @section When to Use Quick Mode:
#' Quick mode is recommended when:
#' \itemize{
#'   \item Using model-based inference (GLMMs, GAMs, occupancy models)
#'   \item PSU size uniformity is not critical for your analysis
#'   \item You want faster results (1-2 minutes vs 10+ minutes)
#'   \item Spatial balance is more important than equal PSU areas
#' }
#'
#' Quick mode typically produces CV of 0.15-0.25 vs 0.10-0.15 for optimized mode.
#' For most ecological applications using model-based inference, this difference
#' is negligible. In practice, quick mode often produces results equal to or
#' better than the optimized mode, especially for irregular study areas.
#'
#' @section Optimization Metrics:
#' \describe{
#'   \item{cv}{Coefficient of variation - often most robust for comparing
#'             relative variation across different PSU sizes}
#'   \item{range}{Difference between largest and smallest PSU}
#'   \item{sd}{Standard deviation of PSU areas}
#' }
#'
#' @section Spatial Balance:
#' Spatial balance is measured using Pielou's evenness index via
#' \code{\link[spsurvey]{sp_balance}}. Higher values indicate better spatial
#' distribution of PSUs. Note that elongated or irregular study areas will
#' naturally have lower Pielou values, which is acceptable as long as PSUs
#' provide good spatial coverage.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(spbalRef)
#'
#' # Load sample frame
#' sample_frame <- st_read("study_area.shp")
#'
#' # Quick mode - fast, for model-based inference (recommended)
#' results_quick <- makePSU(
#'   sf_poly = sample_frame,
#'   psu_area = 300,
#'   n_samples = 25,
#'   island = "North",
#'   quick_mode = TRUE
#' )
#'
#' # Standard mode - optimized PSU sizes
#' results_optimized <- makePSU(
#'   sf_poly = sample_frame,
#'   psu_area = 300,
#'   n_samples = 25,
#'   island = "North",
#'   quick_mode = FALSE
#' )
#'
#' # Compare results
#' print(results_quick)
#' print(results_optimized)
#'
#' # With legacy sites
#' legacy <- st_read("existing_sites.shp")
#'
#' results <- makePSU(
#'   sf_poly = sample_frame,
#'   psu_area = 500,
#'   n_samples = 20,
#'   island = "South",
#'   legacy_sites = legacy,
#'   base_n = 10,
#'   over_n = 5,
#'   quick_mode = TRUE
#' )
#' }
#'
#' @seealso
#' \code{\link[spsurvey]{grts}} for GRTS sampling details
#' \code{\link[anticlust]{balanced_clustering}} for clustering algorithm
#' \code{\link[spbal]{BAS}} for secondary sample generation
#' \code{\link[spsurvey]{sp_balance}} for spatial balance metrics
#'
#' @export
#' @importFrom sf st_union st_as_sf st_transform st_crs st_sample st_area st_cast st_bbox st_intersection st_centroid st_coordinates st_geometry st_set_geometry st_join st_intersects
#' @importFrom dplyr bind_rows group_by summarise mutate arrange filter left_join rename select slice_head tally first "%>%"
#' @importFrom tibble as_tibble
#' @importFrom tidyr drop_na
#' @importFrom stringr str_replace_all
#' @importFrom spsurvey grts sp_rbind sp_balance
#' @importFrom anticlust balanced_clustering
#' @importFrom dismo voronoi
#' @importFrom rlang sym
makePSU <- function(sf_poly,
                    psu_area,
                    cluster_population = NULL,
                    optimization_start = 10000,
                    optimization_step = 3000,
                    optimization_metric = c("cv", "range", "sd"),
                    n_samples,
                    min_radius = NULL,
                    island = c("North", "South"),
                    min.dist = 0,
                    base_n = NULL,
                    over_n = NULL,
                    legacy_sites = NULL,
                    legacy_id_col = "PlotID",
                    grts_seltype = "equal",
                    grts_caty_var = NULL,
                    grts_aux_var = NULL,
                    seed = 25945,
                    plot = TRUE,
                    verbose = TRUE,
                    quick_mode = FALSE) {

  # Input validation
  island <- match.arg(island)
  optimization_metric <- match.arg(optimization_metric)

  if(psu_area <= 0) {
    stop("psu_area must be positive")
  }

  # Set default cluster_population based on mode
  if(is.null(cluster_population)) {
    cluster_population <- if(quick_mode) 50000 else 200000
  }

  if(!quick_mode && cluster_population < 20000) {
    warning("cluster_population may be too small for optimization loop")
  }

  if(!quick_mode && optimization_start >= cluster_population - 10000) {
    stop("optimization_start must be less than cluster_population - 10000")
  }

  # Set seed for replicability
  set.seed(seed)

  if(verbose) {
    message("Preparing sample frame...")
    if(quick_mode) {
      message("Running in QUICK MODE - single-pass clustering")
    }
  }

  # Save original polygon
  sf_poly_original <- sf_poly

  # Combine the polygon
  sf_poly <- st_union(sf_poly)
  sf_poly <- st_as_sf(as_Spatial(sf_poly))

  # Calculate number of PSUs needed
  n_areas <- round(as.numeric(sum(st_area(sf_poly)/10000)/psu_area))

  if(n_areas < 2) {
    stop("Sample frame too small for specified psu_area. Increase psu_area or sample frame size.")
  }

  # Cache bounding box reference
  bb_ref <- if(island == "North") {
    spbalRef::north_nz
  } else {
    spbalRef::south_nz
  }

  # CHECK IF POLYGON IS IN THE SPECIFIED ISLAND BOUNDING BOX
  # Transform polygon to bbox CRS for comparison
  sf_poly_check <- st_transform(sf_poly, st_crs(bb_ref$bb))
  poly_bbox <- st_bbox(sf_poly_check)

  # bb_ref$bb is a polygon, not a bbox - extract its bbox
  master_bbox <- st_bbox(bb_ref$bb)

  # Convert both to simple numeric vectors
  poly_bbox_vec <- c(
    xmin = as.numeric(poly_bbox["xmin"]),
    ymin = as.numeric(poly_bbox["ymin"]),
    xmax = as.numeric(poly_bbox["xmax"]),
    ymax = as.numeric(poly_bbox["ymax"])
  )

  master_bbox_vec <- c(
    xmin = as.numeric(master_bbox["xmin"]),
    ymin = as.numeric(master_bbox["ymin"]),
    xmax = as.numeric(master_bbox["xmax"]),
    ymax = as.numeric(master_bbox["ymax"])
  )

  # Check if polygon overlaps with master sample bounding box
  poly_overlaps <- !(poly_bbox_vec["xmin"] > master_bbox_vec["xmax"] ||
                       poly_bbox_vec["xmax"] < master_bbox_vec["xmin"] ||
                       poly_bbox_vec["ymin"] > master_bbox_vec["ymax"] ||
                       poly_bbox_vec["ymax"] < master_bbox_vec["ymin"])

  if(!poly_overlaps) {
    # Try to detect which island it actually is
    other_island <- if(island == "North") "South" else "North"
    other_bb_ref <- if(island == "North") spbalRef::south_nz else spbalRef::north_nz
    sf_poly_check_other <- st_transform(sf_poly, st_crs(other_bb_ref$bb))
    poly_bbox_other <- st_bbox(sf_poly_check_other)
    other_bbox <- st_bbox(other_bb_ref$bb)  # Extract bbox from polygon

    # Convert to simple numeric vectors
    poly_bbox_other_vec <- c(
      xmin = as.numeric(poly_bbox_other["xmin"]),
      ymin = as.numeric(poly_bbox_other["ymin"]),
      xmax = as.numeric(poly_bbox_other["xmax"]),
      ymax = as.numeric(poly_bbox_other["ymax"])
    )

    other_bbox_vec <- c(
      xmin = as.numeric(other_bbox["xmin"]),
      ymin = as.numeric(other_bbox["ymin"]),
      xmax = as.numeric(other_bbox["xmax"]),
      ymax = as.numeric(other_bbox["ymax"])
    )

    other_overlaps <- !(poly_bbox_other_vec["xmin"] > other_bbox_vec["xmax"] ||
                          poly_bbox_other_vec["xmax"] < other_bbox_vec["xmin"] ||
                          poly_bbox_other_vec["ymin"] > other_bbox_vec["ymax"] ||
                          poly_bbox_other_vec["ymax"] < other_bbox_vec["ymin"])

    if(other_overlaps) {
      stop(sprintf(
        "Sample frame does not overlap with %s Island bounding box, but DOES overlap with %s Island.\nPlease set island = '%s'.",
        island, other_island, other_island
      ))
    } else {
      stop(sprintf(
        "Sample frame does not overlap with %s Island bounding box.\n  Your polygon extent: xmin=%.0f, xmax=%.0f, ymin=%.0f, ymax=%.0f\n  %s Island bbox: xmin=%.0f, xmax=%.0f, ymin=%.0f, ymax=%.0f\n  Check your polygon's coordinate system.",
        island,
        poly_bbox_vec["xmin"], poly_bbox_vec["xmax"], poly_bbox_vec["ymin"], poly_bbox_vec["ymax"],
        island,
        master_bbox_vec["xmin"], master_bbox_vec["xmax"], master_bbox_vec["ymin"], master_bbox_vec["ymax"]
      ))
    }
  }

  if(verbose) {
    message(sprintf("Creating %d PSUs of ~%.1f ha each", n_areas, psu_area))
    message(sprintf("Using %s Island master sample", island))
  }

  # Transform to appropriate CRS
  sf_poly <- st_transform(sf_poly, st_crs(bb_ref$bb))

  if(verbose) message("Generating random points for clustering...")

  # Generate random points for clustering
  points_rnd <- st_sample(sf_poly, size = cluster_population, type = "random")

  # QUICK MODE: Single-pass clustering
  if(quick_mode) {

    if(verbose) message("Performing balanced clustering (single pass)...")

    # Make matrix of coordinates FIRST
    points <- do.call(rbind, st_geometry(points_rnd)) %>%
      as_tibble() %>%
      setNames(c("lon", "lat")) %>%
      tidyr::drop_na()

    # Calculate how many points we can use (must be divisible by n_areas)
    n_points_final <- floor(nrow(points) / n_areas) * n_areas

    if(n_points_final < n_areas) {
      stop(sprintf(
        "Not enough valid points (%d) for %d clusters. Increase cluster_population.",
        nrow(points), n_areas
      ))
    }

    # Subset to exact number needed
    points <- points[1:n_points_final, ]

    if(verbose) {
      message(sprintf("  Using %d points (%d per cluster)",
                      n_points_final, n_points_final / n_areas))
    }

    # Balanced clustering
    bc_clust <- anticlust::balanced_clustering(
      x = points,
      K = n_areas,
      method = "centroid"
    )

    bc_clust <- cbind(points, bc_clust)
    bc_clust <- st_as_sf(bc_clust, coords = c("lon", "lat"))
    bc_clust <- bc_clust %>%
      group_by(bc_clust) %>%
      summarise(geometry = st_union(geometry)) %>%
      st_centroid()

    bc_clust <- do.call(rbind, st_geometry(bc_clust)) %>%
      as_tibble() %>%
      setNames(c("lon", "lat")) %>%
      as.matrix()

    # Create voronoi polygons
    if(verbose) message("Creating PSU boundaries...")
    sf_poly_cast <- st_cast(sf_poly, "POLYGON")
    voronoi_polys <- dismo::voronoi(bc_clust, ext = st_bbox(sf_poly_cast))
    voronoi_polys <- st_as_sf(voronoi_polys)
    voronoi_polys <- st_set_crs(voronoi_polys, st_crs(sf_poly))

    # Clip to sample frame
    voronoi_sf <- st_as_sf(voronoi_polys)
    equal_areas <- st_intersection(voronoi_sf,
                                   st_transform(sf_poly, st_crs(voronoi_sf)))

    # Calculate areas
    equal_areas$area <- st_area(equal_areas)/10000
    equal_areas <- equal_areas %>%
      group_by(id) %>%
      summarise()
    equal_areas$area <- st_area(equal_areas)/10000

    # Make centroids for PSU selection
    equal_area_centres <- st_centroid(equal_areas)
    equal_area_centres <- st_cast(equal_area_centres, "POINT")
    coords <- st_coordinates(equal_area_centres)
    equal_area_centres$X <- coords[, 1]
    equal_area_centres$Y <- coords[, 2]

    # Check spatial balance
    pielou <- sp_balance(
      object = equal_area_centres,
      sframe = sf_poly,
      stratum_var = NULL,
      ip = NULL,
      metrics = "pielou",
      extents = FALSE
    )

    psu_size <- n_points_final

    # Set optimization outputs to NULL for quick mode
    dat.sum1 <- NULL
    dat.sum <- NULL

    if(verbose) {
      message(sprintf("PSUs created with n_points = %d", psu_size))
      message(sprintf("PSU area CV = %.4f (quick mode - not optimized)",
                      sd(as.numeric(equal_areas$area)) / mean(as.numeric(equal_areas$area))))
    }

  } else {
    # STANDARD MODE: Optimization loop

    # Pre-allocate storage with lists for better performance
    optimization_seq <- seq(optimization_start,
                            length(points_rnd) - 10000,
                            by = optimization_step)

    n_iterations <- length(optimization_seq)
    dat_list <- vector("list", n_iterations)
    PSU_list <- vector("list", n_iterations)
    PSU_centres_list <- vector("list", n_iterations)

    if(verbose) {
      message(sprintf("Running optimization loop (%d iterations)...", n_iterations))
    }

    # Iterate clustering at different point densities to find optimal PSU size
    for(idx in seq_along(optimization_seq)) {
      i <- optimization_seq[idx]

      if(verbose && idx %% 5 == 0) {
        message(sprintf("  Testing iteration %d/%d (n_points = %d)",
                        idx, n_iterations, i))
      }

      # Find number that divides evenly
      n_points <- ceiling(i/n_areas) * n_areas
      points_rnd1 <- points_rnd[1:n_points]

      # Make matrix of coordinates
      points <- do.call(rbind, st_geometry(points_rnd1)) %>%
        as_tibble() %>%
        setNames(c("lon", "lat")) %>%
        tidyr::drop_na()

      # Ensure exact divisibility after dropping NAs
      n_points_adjusted <- floor(nrow(points) / n_areas) * n_areas
      points <- points[1:n_points_adjusted, ]

      # Balanced clustering - produces more even-sized polygons
      bc_clust <- anticlust::balanced_clustering(
        x = points,
        K = n_areas,
        method = "centroid"
      )

      bc_clust <- cbind(points, bc_clust)
      bc_clust <- st_as_sf(bc_clust, coords = c("lon", "lat"))
      bc_clust <- bc_clust %>%
        group_by(bc_clust) %>%
        summarise(geometry = st_union(geometry)) %>%
        st_centroid()

      bc_clust <- do.call(rbind, st_geometry(bc_clust)) %>%
        as_tibble() %>%
        setNames(c("lon", "lat")) %>%
        as.matrix()

      # Create voronoi polygons
      sf_poly_cast <- st_cast(sf_poly, "POLYGON")
      voronoi_polys <- dismo::voronoi(bc_clust, ext = st_bbox(sf_poly_cast))
      voronoi_polys <- st_as_sf(voronoi_polys)
      voronoi_polys <- st_set_crs(voronoi_polys, st_crs(sf_poly))

      # Clip to sample frame
      voronoi_sf <- st_as_sf(voronoi_polys)
      equal_areas <- st_intersection(voronoi_sf,
                                     st_transform(sf_poly, st_crs(voronoi_sf)))

      # Calculate areas
      equal_areas$area <- st_area(equal_areas)/10000
      equal_areas <- equal_areas %>%
        group_by(id) %>%
        summarise()
      equal_areas$area <- st_area(equal_areas)/10000
      equal_areas$iteration <- n_points

      # Make centroids for PSU selection
      equal_area_centres <- st_centroid(equal_areas)
      equal_area_centres <- st_cast(equal_area_centres, "POINT")
      coords <- st_coordinates(equal_area_centres)
      equal_area_centres$X <- coords[, 1]
      equal_area_centres$Y <- coords[, 2]
      equal_area_centres$iteration <- n_points

      # Check spatial balance
      pielou <- sp_balance(
        object = equal_area_centres,
        sframe = sf_poly,
        stratum_var = NULL,
        ip = NULL,
        metrics = "pielou",
        extents = FALSE
      )

      # Store summary statistics
      dat <- st_set_geometry(equal_areas, NULL)
      dat$n_points <- n_points
      dat$var <- var(dat$area)
      dat$mean <- mean(dat$area)
      dat$sd <- sd(dat$area)
      dat$cv <- sd(dat$area) / mean(dat$area)
      dat$pielou <- pielou$value

      dat_list[[idx]] <- dat
      PSU_list[[idx]] <- equal_areas
      PSU_centres_list[[idx]] <- equal_area_centres
    }

    # Combine results
    if(verbose) message("Combining optimization results...")
    dat.sum <- bind_rows(dat_list)
    PSU_store <- bind_rows(PSU_list)
    PSU_centres_store <- bind_rows(PSU_centres_list)

    # Clear large objects to free memory
    rm(dat_list, PSU_list, PSU_centres_list, points_rnd)
    gc()

    # Find optimal iteration based on selected metric
    dat.sum1 <- dat.sum %>%
      group_by(n_points) %>%
      summarise(
        mean = mean(area),
        sd = sd(area),
        cv = mean(cv),
        median = median(area),
        range_min = min(area),
        range_max = max(area),
        iqr = IQR(area),
        pielou = first(pielou)
      ) %>%
      mutate(range = as.numeric(range_max) - as.numeric(range_min))

    # Select optimal based on metric
    dat.sum1 <- dat.sum1 %>%
      arrange(!!sym(optimization_metric))

    # Get optimal PSU size
    psu_size <- dat.sum1[1, ]$n_points
    pielou <- dat.sum1[1, ]$pielou

    if(verbose) {
      message(sprintf("Optimal iteration: n_points = %d (%s = %.4f)",
                      psu_size,
                      optimization_metric,
                      dat.sum1[1, ][[optimization_metric]]))
    }

    # Extract optimal iteration
    equal_areas <- PSU_store %>% filter(iteration == psu_size)
    equal_area_centres <- PSU_centres_store %>% filter(iteration == psu_size)
  }

  # Handle legacy sites if provided (common to both modes)
  if(!is.null(legacy_sites)) {
    if(verbose) message("Processing legacy sites...")

    # Check if legacy sites are in frame
    legacy_sites <- st_transform(legacy_sites, st_crs(sf_poly))
    legacy_in_frame <- st_intersects(legacy_sites, sf_poly, sparse = FALSE)
    if(sum(legacy_in_frame) < nrow(legacy_sites)) {
      warning(sprintf("%d legacy site(s) outside sample frame and will be ignored",
                      nrow(legacy_sites) - sum(legacy_in_frame)))
      legacy_sites <- legacy_sites[legacy_in_frame[,1], ]
    }

    legacy_filter <- st_join(
      st_transform(legacy_sites, st_crs(equal_areas)),
      equal_areas[, "id"]
    ) %>%
      tidyr::drop_na(id)

    if(nrow(legacy_filter) > 0) {
      equal_area_centres <- equal_area_centres %>%
        left_join(
          st_set_geometry(legacy_filter, NULL)[, c(legacy_id_col, "id")],
          by = "id"
        )
    } else {
      warning("No legacy sites intersected with PSUs")
    }
  }

  # Sample PSUs using GRTS (common to both modes)
  if(is.null(base_n)) {
    base_n <- n_areas
  }

  if(is.null(over_n) || over_n == 0) {
    over_n <- 0
    use_over <- FALSE
  } else {
    use_over <- TRUE
  }

  if(verbose) message("Selecting PSUs using GRTS...")

  set.seed(seed + 1)

  # Call GRTS - only include n_over if we actually want an oversample
  if(use_over) {
    equal_area_sample <- grts(
      sframe = equal_area_centres,
      n_base = base_n,
      n_over = over_n,
      legacy_var = if(!is.null(legacy_sites) && legacy_id_col %in% names(equal_area_centres)) {
        legacy_id_col
      } else {
        NULL
      },
      seltype = grts_seltype,
      caty_var = grts_caty_var,
      aux_var = grts_aux_var,
      mindis = min.dist,
      stratum_var = NULL,
      maxtry = 500
    )
    equal_area_sample <- sp_rbind(equal_area_sample)
  } else {
    # Don't include n_over at all when we don't want oversample
    equal_area_sample <- grts(
      sframe = equal_area_centres,
      n_base = base_n,
      legacy_var = if(!is.null(legacy_sites) && legacy_id_col %in% names(equal_area_centres)) {
        legacy_id_col
      } else {
        NULL
      },
      seltype = grts_seltype,
      caty_var = grts_caty_var,
      aux_var = grts_aux_var,
      mindis = min.dist,
      stratum_var = NULL,
      maxtry = 500
    )
    # Just use the base sites when there's no oversample
    equal_area_sample <- equal_area_sample$sites_base
  }

  # Join sample info back to PSUs
  equal_areas <- equal_areas %>%
    st_join(equal_area_sample[, c("siteID", "siteuse")]) %>%
    rename(PSU_sample_order = "siteID")

  if(verbose) message("Generating secondary sample points...")

  # Generate secondary sample points using BAS
  points_rnd <- spbal::BAS(
    shapefile = sf_poly,
    n = 10000,
    boundingbox = bb_ref$bb,
    seeds = bb_ref$seed,
    verbose = FALSE,
    minRadius = min_radius
  )

  # Add unique IDs
  if (!is.null(min_radius)) {
    bas_result <- points_rnd$minRadius
    bas_result$uniqueID <- paste0(bb_ref$mastersample_name, bas_result$SiteID)
    points_rnd <- bas_result
  } else {
    points_rnd <- spbalRef::uniqueID(bb_ref, points_rnd$sample)
    points_rnd$uniqueID <- paste0(points_rnd$sample, points_rnd$SiteID)
  }

  # Only join against PSUs that were actually selected (Base/Over). Joining
  # against every PSU (including unselected ones, which have
  # PSU_sample_order = NA) lumped BAS points that fell in non-selected PSUs
  # into a bogus PSU_sample_order = NA group that then got sampled from too.
  selected_psus <- equal_areas %>%
    filter(!is.na(PSU_sample_order))

  # Select secondary sample points within each PSU
  sample_points <- points_rnd %>%
    st_join(dplyr::select(selected_psus, PSU_sample_order, siteuse)) %>%
    filter(!is.na(PSU_sample_order)) %>%
    group_by(PSU_sample_order) %>%
    mutate(
      Transect_sample_order = rank(SiteID),
      PSU_sample_order = as.numeric(
        str_replace_all(PSU_sample_order, "[-[:alpha:]]", "")
      )
    ) %>%
    arrange(PSU_sample_order) %>%
    slice_head(n = n_samples) %>%
    dplyr::select(PSU_sample_order, Transect_sample_order,
                  siteuse, uniqueID)

  # Check if all selected PSUs have enough sample points. Compare against the
  # full set of selected PSU IDs -- not just the groups present in
  # sample_points -- so a PSU that received ZERO BAS points (which never
  # forms a group_by() group at all) is still detected instead of silently
  # vanishing from secondary_sample_units.
  selected_ids <- as.numeric(
    str_replace_all(selected_psus$PSU_sample_order, "[-[:alpha:]]", "")
  )
  points_per_psu <- tibble::tibble(PSU_sample_order = selected_ids) %>%
    left_join(
      sf::st_drop_geometry(sample_points) %>%
        group_by(PSU_sample_order) %>%
        tally(),
      by = "PSU_sample_order"
    ) %>%
    mutate(n = ifelse(is.na(n), 0L, n))

  n_empty <- sum(points_per_psu$n == 0)
  n_short <- sum(points_per_psu$n > 0 & points_per_psu$n < n_samples)

  if(n_empty > 0) {
    warning(sprintf(
      "%d PSU(s) received ZERO secondary sample points and are missing entirely from secondary_sample_units (PSU_sample_order: %s)",
      n_empty,
      paste(sort(points_per_psu$PSU_sample_order[points_per_psu$n == 0]), collapse = ", ")
    ))
  }
  if(n_short > 0) {
    warning(sprintf(
      "%d PSU(s) have fewer than %d sample points (PSU_sample_order: %s)",
      n_short, n_samples,
      paste(sort(points_per_psu$PSU_sample_order[points_per_psu$n > 0 & points_per_psu$n < n_samples]), collapse = ", ")
    ))
  }

  # Calculate final statistics
  psu_area_stats <- list(
    target_area = psu_area,
    mean = mean(as.numeric(equal_areas$area)),
    sd = sd(as.numeric(equal_areas$area)),
    cv = sd(as.numeric(equal_areas$area)) / mean(as.numeric(equal_areas$area)),
    median = median(as.numeric(equal_areas$area)),
    min = min(as.numeric(equal_areas$area)),
    max = max(as.numeric(equal_areas$area)),
    range = max(as.numeric(equal_areas$area)) - min(as.numeric(equal_areas$area))
  )

  # Create output list
  samp_objects <- list(
    primary_sample_units = equal_areas,
    secondary_sample_units = sample_points,
    optimization_summary = dat.sum1,
    optimization_data = dat.sum,
    optimal_n_points = psu_size,
    pielou_evenness = if(quick_mode) pielou$value else dat.sum1[1, ]$pielou,
    psu_area_stats = psu_area_stats,
    optimization_metric = if(quick_mode) NA_character_ else optimization_metric,
    settings = list(
      psu_area = psu_area,
      cluster_population = cluster_population,
      n_samples = n_samples,
      min_radius = min_radius,
      island = island,
      seed = seed,
      base_n = base_n,
      over_n = over_n,
      quick_mode = quick_mode
    )
  )

  # Create plots
  if(plot) {
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

    # 1. Histogram of PSU areas
    hist(
      as.numeric(equal_areas$area),
      main = sprintf("PSU Areas (n = %d)%s", n_areas,
                     if(quick_mode) " - Quick Mode" else ""),
      xlab = "Area (hectares)",
      col = "lightblue",
      border = "white"
    )
    abline(v = psu_area, col = "red", lwd = 2, lty = 2)
    legend("topright", legend = sprintf("Target: %.1f ha", psu_area),
           col = "red", lty = 2, lwd = 2, bty = "n")

    # 2. Map of PSUs with sample order
    plot(st_geometry(equal_areas),
         col = ifelse(equal_areas$siteuse == "Base", "lightgreen",
                      ifelse(equal_areas$siteuse == "Over", "lightyellow", "lightgray")),
         border = "gray40",
         main = "Primary Sample Units")
    legend("topright",
           legend = c("Base", "Over", "Not selected"),
           fill = c("lightgreen", "lightyellow", "lightgray"),
           border = "gray40",
           bty = "n")

    par(mfrow = c(1, 1))
  }

  if(verbose) {
    message("\nSummary:")
    message(sprintf("  Mode: %s", if(quick_mode) "Quick (single-pass)" else "Optimized"))
    message(sprintf("  PSUs created: %d", n_areas))
    message(sprintf("  PSUs selected (base): %d", base_n))
    if(!is.null(over_n) && over_n > 0) {
      message(sprintf("  PSUs selected (oversample): %d", over_n))
    }
    message(sprintf("  Mean PSU area: %.1f ha (SD = %.1f, CV = %.3f)",
                    psu_area_stats$mean, psu_area_stats$sd, psu_area_stats$cv))
    message(sprintf("  PSU area range: %.1f - %.1f ha",
                    psu_area_stats$min, psu_area_stats$max))
    message(sprintf("  Pielou evenness: %.4f", samp_objects$pielou_evenness))
  }

  class(samp_objects) <- c("psu", "list")
  return(samp_objects)
}

#' Print Method for psu Objects
#'
#' @description
#' Prints a summary of a two-stage spatial sampling design created by
#' \code{\link{makePSU}}.
#'
#' @param x An object of class \code{"psu"} created by
#'   \code{\link{makePSU}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return
#' Invisibly returns the input object \code{x}.
#'
#' @export
print.psu <- function(x, ...) {
  cat("Two-stage spatial sample design\n")
  cat("================================\n\n")
  cat(sprintf("Mode: %s\n",
              if(isTRUE(x$settings$quick_mode)) "Quick (single-pass)" else "Optimized"))
  cat(sprintf("Primary sampling units: %d created, %d selected\n",
              nrow(x$primary_sample_units),
              sum(!is.na(x$primary_sample_units$PSU_sample_order))))
  cat(sprintf("Secondary sampling units: %d\n",
              nrow(x$secondary_sample_units)))
  cat(sprintf("\nPSU area statistics:\n"))
  cat(sprintf("  Target: %.1f ha\n", x$psu_area_stats$target_area))
  cat(sprintf("  Mean: %.1f ha (SD = %.1f)\n",
              x$psu_area_stats$mean, x$psu_area_stats$sd))
  cat(sprintf("  CV: %.3f\n", x$psu_area_stats$cv))
  cat(sprintf("  Range: %.1f - %.1f ha\n",
              x$psu_area_stats$min, x$psu_area_stats$max))
  cat(sprintf("\nSpatial balance (Pielou): %.4f\n", x$pielou_evenness))
  if(!is.na(x$optimization_metric)) {
    cat(sprintf("Optimized using: %s\n", x$optimization_metric))
  }
  invisible(x)
}
