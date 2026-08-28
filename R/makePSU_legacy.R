#' Create Two-Stage Spatial Sampling Design with Legacy-Seeded PSUs (experimental)
#'
#' @description
#' Experimental variant of \code{\link{makePSU}} that guarantees each legacy
#' monitoring site anchors its own primary sampling unit (PSU). Rather than
#' clustering PSU centroids purely at random and only checking afterwards
#' which PSU a legacy site happened to fall into, this function reserves one
#' cluster per legacy site up front: the points nearest to each legacy site
#' are assigned to that site's cluster first, and the balanced clustering
#' algorithm is only run on the remaining points to fill the rest of the
#' PSUs. The legacy PSU's centroid is the centroid of its assigned points
#' (not the literal legacy coordinate, since legacy coordinates are often
#' themselves a mean of historical transect points).
#'
#' Legacy PSUs are always included in the base sample; GRTS fills the
#' remaining base/oversample slots from the non-legacy PSUs.
#'
#' Only single-pass ("quick mode") clustering is implemented for now -- no
#' optimization loop over PSU size. This is a standalone function for
#' testing the legacy-seeding approach before it is folded into
#' \code{\link{makePSU}}.
#'
#' @inheritParams makePSU
#' @param legacy_sites An \code{sf} object containing legacy monitoring sites
#'   to seed PSUs around. Required (unlike \code{makePSU}, where it is
#'   optional).
#'
#' @return An object of class \code{"psu"}, structured the same as the output
#'   of \code{\link{makePSU}}. \code{primary_sample_units} additionally has a
#'   \code{is_legacy} logical column and the \code{legacy_id_col} column
#'   identifying which legacy site(s) seeded each PSU.
#'
#' @export
#' @importFrom sf st_union st_as_sf st_transform st_crs st_sample st_area
#'   st_cast st_bbox st_intersection st_centroid st_coordinates st_geometry
#'   st_set_geometry st_join st_intersects st_distance
#' @importFrom dplyr bind_rows group_by summarise mutate arrange filter
#'   left_join rename select slice_head tally first n_distinct
#' @importFrom tibble as_tibble
#' @importFrom tidyr drop_na
#' @importFrom stringr str_replace_all
#' @importFrom spsurvey grts sp_rbind sp_balance
#' @importFrom anticlust balanced_clustering
#' @importFrom dismo voronoi
makePSU_legacy <- function(sf_poly,
                           psu_area,
                           cluster_population = 50000,
                           n_samples,
                           min_radius = NULL,
                           island = c("North", "South"),
                           min.dist = 0,
                           base_n = NULL,
                           over_n = NULL,
                           legacy_sites,
                           legacy_id_col = "PlotID",
                           grts_seltype = "equal",
                           grts_caty_var = NULL,
                           grts_aux_var = NULL,
                           seed = 25945,
                           plot = TRUE,
                           verbose = TRUE) {

  # Input validation
  island <- match.arg(island)

  if(psu_area <= 0) {
    stop("psu_area must be positive")
  }

  if(is.null(legacy_sites) || nrow(legacy_sites) == 0) {
    stop("legacy_sites is required for makePSU_legacy() (use makePSU() if you have none)")
  }

  if(!(legacy_id_col %in% names(legacy_sites))) {
    stop(sprintf("legacy_id_col '%s' not found in legacy_sites", legacy_id_col))
  }

  # Set seed for replicability
  set.seed(seed)

  if(verbose) message("Preparing sample frame...")

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

  if(verbose) {
    message(sprintf("Creating %d PSUs of ~%.1f ha each", n_areas, psu_area))
    message(sprintf("Using %s Island master sample", island))
  }

  # Transform to appropriate CRS
  sf_poly <- st_transform(sf_poly, st_crs(bb_ref$bb))
  legacy_sites <- st_transform(legacy_sites, st_crs(sf_poly))

  # Drop legacy sites outside the frame
  legacy_in_frame <- st_intersects(legacy_sites, sf_poly, sparse = FALSE)[, 1]
  if(sum(legacy_in_frame) < nrow(legacy_sites)) {
    warning(sprintf("%d legacy site(s) outside sample frame and will be ignored",
                    nrow(legacy_sites) - sum(legacy_in_frame)))
    legacy_sites <- legacy_sites[legacy_in_frame, ]
  }

  # Collapse legacy sites down to one anchor per site, but keep note of
  # duplicates so multiple legacy sites that fall very close together don't
  # each try to claim their own PSU quota independently.
  n_legacy <- nrow(legacy_sites)

  if(n_legacy == 0) {
    stop("No legacy sites fall within the sample frame")
  }

  if(n_legacy >= n_areas) {
    stop(sprintf(
      "Number of legacy sites (%d) must be less than the number of PSUs (%d). Decrease psu_area or reduce legacy_sites.",
      n_legacy, n_areas
    ))
  }

  if(verbose) {
    message(sprintf("Seeding %d PSU(s) around legacy sites (%d PSU(s) total)",
                    n_legacy, n_areas))
  }

  if(verbose) message("Generating random points for clustering...")

  # Generate random points for clustering
  points_rnd <- st_sample(sf_poly, size = cluster_population, type = "random")

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

  points <- points[1:n_points_final, ]
  quota <- n_points_final / n_areas

  if(verbose) {
    message(sprintf("  Using %d points (%d per cluster)", n_points_final, quota))
  }

  # --- Legacy-seeded clustering -------------------------------------------
  # Reserve clusters 1:n_legacy for legacy sites: assign each legacy site the
  # nearest `quota` still-unassigned points. Remaining points are balanced-
  # clustered as usual to fill clusters (n_legacy + 1):n_areas.

  legacy_coords <- st_coordinates(legacy_sites)
  point_coords <- as.matrix(points[, c("lon", "lat")])

  cluster_assignment <- rep(NA_integer_, nrow(points))
  unassigned <- rep(TRUE, nrow(points))

  if(verbose) message("Assigning points to legacy-seeded clusters...")

  for(i in seq_len(n_legacy)) {
    idx_pool <- which(unassigned)

    d <- sqrt(
      (point_coords[idx_pool, 1] - legacy_coords[i, 1])^2 +
      (point_coords[idx_pool, 2] - legacy_coords[i, 2])^2
    )

    take_n <- min(quota, length(idx_pool))
    nearest_idx <- idx_pool[order(d)[seq_len(take_n)]]

    cluster_assignment[nearest_idx] <- i
    unassigned[nearest_idx] <- FALSE
  }

  # Balanced-cluster whatever points remain into the non-legacy PSUs
  remaining_idx <- which(unassigned)
  n_remaining_clusters <- n_areas - n_legacy

  if(verbose) message("Performing balanced clustering on remaining points...")

  bc_remaining <- anticlust::balanced_clustering(
    x = points[remaining_idx, ],
    K = n_remaining_clusters,
    method = "centroid"
  )

  cluster_assignment[remaining_idx] <- bc_remaining + n_legacy

  bc_clust <- cluster_assignment

  # --- Build centroids and Voronoi PSUs (same as makePSU quick_mode) -----

  bc_clust_sf <- cbind(points, bc_clust) %>%
    st_as_sf(coords = c("lon", "lat"))
  bc_clust_sf <- bc_clust_sf %>%
    group_by(bc_clust) %>%
    summarize(geometry = st_union(geometry)) %>%
    st_centroid()

  bc_clust_mat <- do.call(rbind, st_geometry(bc_clust_sf)) %>%
    as_tibble() %>%
    setNames(c("lon", "lat")) %>%
    as.matrix()

  if(verbose) message("Creating PSU boundaries...")
  sf_poly_cast <- st_cast(sf_poly, "POLYGON")
  voronoi_polys <- dismo::voronoi(bc_clust_mat, ext = st_bbox(sf_poly_cast))
  voronoi_polys <- st_as_sf(voronoi_polys)
  voronoi_polys <- st_set_crs(voronoi_polys, st_crs(sf_poly))

  # Clip to sample frame
  voronoi_sf <- st_as_sf(voronoi_polys)
  equal_areas <- st_intersection(voronoi_sf,
                                 st_transform(sf_poly, st_crs(voronoi_sf)))

  equal_areas$area <- st_area(equal_areas)/10000
  equal_areas <- equal_areas %>%
    group_by(id) %>%
    summarize()
  equal_areas$area <- st_area(equal_areas)/10000

  # Make centroids for PSU selection
  equal_area_centres <- st_centroid(equal_areas)
  equal_area_centres <- st_cast(equal_area_centres, "POINT")
  coords <- st_coordinates(equal_area_centres)
  equal_area_centres$X <- coords[, 1]
  equal_area_centres$Y <- coords[, 2]

  # Tag legacy PSUs: attach the legacy_id_col of whichever legacy site(s)
  # actually fall within each PSU polygon (should normally be the site(s)
  # that seeded it, since it was clustered around their nearest points).
  legacy_join <- st_join(
    legacy_sites[, legacy_id_col],
    equal_areas[, "id"]
  ) %>%
    st_set_geometry(NULL) %>%
    tidyr::drop_na(id) %>%
    group_by(id) %>%
    summarise(
      !!legacy_id_col := paste(unique(.data[[legacy_id_col]]), collapse = "; "),
      .groups = "drop"
    )

  equal_area_centres <- equal_area_centres %>%
    left_join(legacy_join, by = "id") %>%
    mutate(is_legacy = !is.na(.data[[legacy_id_col]]))

  n_tagged <- sum(equal_area_centres$is_legacy)
  if(n_tagged < n_legacy) {
    warning(sprintf(
      "%d legacy site(s) did not fall inside their seeded PSU polygon (only %d of %d PSUs tagged). Voronoi boundaries may not follow cluster assignment exactly.",
      n_legacy - n_tagged, n_tagged, n_legacy
    ))
  }

  # Check spatial balance
  pielou <- sp_balance(
    object = equal_area_centres,
    sframe = sf_poly,
    stratum_var = NULL,
    ip = NULL,
    metrics = "pielou",
    extents = FALSE
  )

  if(verbose) {
    message(sprintf("PSU area CV = %.4f (legacy-seeded, single pass)",
                    sd(as.numeric(equal_areas$area)) / mean(as.numeric(equal_areas$area))))
  }

  # --- Sample PSUs: legacy PSUs always included, GRTS fills the rest -----

  if(is.null(base_n)) {
    base_n <- n_areas
  }

  if(base_n < n_tagged) {
    stop(sprintf(
      "base_n (%d) is smaller than the number of legacy-tagged PSUs (%d); legacy PSUs must always be included",
      base_n, n_tagged
    ))
  }

  if(is.null(over_n) || over_n == 0) {
    over_n <- 0
    use_over <- FALSE
  } else {
    use_over <- TRUE
  }

  if(verbose) message("Selecting PSUs (legacy PSUs guaranteed, remainder via GRTS)...")

  set.seed(seed + 1)

  legacy_centres <- equal_area_centres %>% filter(is_legacy)
  candidate_centres <- equal_area_centres %>% filter(!is_legacy)

  grts_base_n <- base_n - n_tagged

  siteID_all <- character(nrow(equal_area_centres))
  siteuse_all <- character(nrow(equal_area_centres))

  if(grts_base_n > 0 || use_over) {
    if(use_over) {
      grts_sample <- grts(
        sframe = candidate_centres,
        n_base = grts_base_n,
        n_over = over_n,
        seltype = grts_seltype,
        caty_var = grts_caty_var,
        aux_var = grts_aux_var,
        mindis = min.dist,
        stratum_var = NULL,
        maxtry = 500
      )
      grts_sample <- sp_rbind(grts_sample)
    } else {
      grts_sample <- grts(
        sframe = candidate_centres,
        n_base = grts_base_n,
        seltype = grts_seltype,
        caty_var = grts_caty_var,
        aux_var = grts_aux_var,
        mindis = min.dist,
        stratum_var = NULL,
        maxtry = 500
      )
      grts_sample <- grts_sample$sites_base
    }
  } else {
    grts_sample <- NULL
  }

  # Legacy PSUs get sample order 1:n_tagged (always Base), GRTS-selected
  # PSUs are numbered after them.
  legacy_centres$PSU_sample_order <- paste0("Legacy-", seq_len(nrow(legacy_centres)))
  legacy_centres$siteuse <- "Base"

  if(!is.null(grts_sample) && nrow(grts_sample) > 0) {
    # spsurvey::grts() does not preserve the input frame's "id" attribute --
    # it returns its own siteID numbering -- so match selected sites back to
    # candidate_centres spatially (as makePSU() does), not by an id join.
    grts_join <- candidate_centres %>%
      st_join(dplyr::select(st_transform(grts_sample, st_crs(candidate_centres)),
                             siteuse)) %>%
      filter(!is.na(siteuse)) %>%
      arrange(siteuse) %>%
      mutate(PSU_sample_order = paste0(
        ifelse(siteuse == "Base", "Base-", "Over-"),
        row_number()
      ))
  } else {
    grts_join <- candidate_centres[0, ]
    grts_join$PSU_sample_order <- character(0)
    grts_join$siteuse <- character(0)
  }

  equal_area_sample <- bind_rows(
    st_set_geometry(legacy_centres, NULL)[, c("id", "PSU_sample_order", "siteuse")],
    st_set_geometry(grts_join, NULL)[, c("id", "PSU_sample_order", "siteuse")]
  )

  # Join sample info back to PSUs
  equal_areas <- equal_areas %>%
    left_join(equal_area_sample, by = "id")

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

  if(!is.null(min_radius)) {
    bas_result <- points_rnd$minRadius
    bas_result$uniqueID <- paste0(bb_ref$mastersample_name, bas_result$SiteID)
    points_rnd <- bas_result
  } else {
    points_rnd <- spbalRef::uniqueID(bb_ref, points_rnd$sample)
    points_rnd$uniqueID <- paste0(points_rnd$sample, points_rnd$SiteID)
  }

  selected_psus <- equal_areas %>%
    filter(!is.na(PSU_sample_order))

  sample_points <- points_rnd %>%
    st_join(dplyr::select(selected_psus, PSU_sample_order, siteuse)) %>%
    filter(!is.na(PSU_sample_order)) %>%
    group_by(PSU_sample_order) %>%
    mutate(Transect_sample_order = rank(SiteID)) %>%
    arrange(PSU_sample_order) %>%
    slice_head(n = n_samples) %>%
    dplyr::select(PSU_sample_order, Transect_sample_order,
                  siteuse, uniqueID)

  points_per_psu <- selected_psus %>%
    st_set_geometry(NULL) %>%
    dplyr::select(PSU_sample_order) %>%
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

  samp_objects <- list(
    primary_sample_units = equal_areas,
    secondary_sample_units = sample_points,
    optimization_summary = NULL,
    optimization_data = NULL,
    optimal_n_points = n_points_final,
    pielou_evenness = pielou$value,
    psu_area_stats = psu_area_stats,
    optimization_metric = NA_character_,
    settings = list(
      psu_area = psu_area,
      cluster_population = cluster_population,
      n_samples = n_samples,
      min_radius = min_radius,
      island = island,
      seed = seed,
      base_n = base_n,
      over_n = over_n,
      quick_mode = TRUE,
      n_legacy = n_legacy
    )
  )

  if(plot) {
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

    hist(
      as.numeric(equal_areas$area),
      main = sprintf("PSU Areas (n = %d) - Legacy-Seeded", n_areas),
      xlab = "Area (hectares)",
      col = "lightblue",
      border = "white"
    )
    abline(v = psu_area, col = "red", lwd = 2, lty = 2)
    legend("topright", legend = sprintf("Target: %.1f ha", psu_area),
           col = "red", lty = 2, lwd = 2, bty = "n")

    plot_col <- ifelse(is.na(equal_areas$siteuse), "lightgray",
                 ifelse(equal_areas$siteuse == "Base" &
                          grepl("^Legacy-", equal_areas$PSU_sample_order), "orange",
                 ifelse(equal_areas$siteuse == "Base", "lightgreen", "lightyellow")))

    plot(st_geometry(equal_areas),
         col = plot_col,
         border = "gray40",
         main = "Primary Sample Units")
    legend("topright",
           legend = c("Legacy (Base)", "Base", "Over", "Not selected"),
           fill = c("orange", "lightgreen", "lightyellow", "lightgray"),
           border = "gray40",
           bty = "n")

    par(mfrow = c(1, 1))
  }

  if(verbose) {
    message("\nSummary:")
    message("  Mode: Legacy-seeded (single-pass)")
    message(sprintf("  PSUs created: %d (%d legacy-seeded)", n_areas, n_legacy))
    message(sprintf("  PSUs selected (base): %d", base_n))
    if(!is.null(over_n) && over_n > 0) {
      message(sprintf("  PSUs selected (oversample): %d", over_n))
    }
    message(sprintf("  Mean PSU area: %.1f ha (SD = %.1f, CV = %.3f)",
                    psu_area_stats$mean, psu_area_stats$sd, psu_area_stats$cv))
    message(sprintf("  PSU area range: %.1f - %.1f ha",
                    psu_area_stats$min, psu_area_stats$max))
    message(sprintf("  Pielou evenness: %.4f", pielou$value))
  }

  class(samp_objects) <- c("psu", "list")
  return(samp_objects)
}
