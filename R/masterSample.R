#' Draw a Master Sample from an NZ Island Master Sample
#'
#' @description
#' Convenience wrapper around \code{\link[spbal]{BAS}} for drawing a
#' spatially balanced sample of points from the New Zealand master sample,
#' with the mastersample \code{uniqueID} (e.g. \code{"NZ_North12345"})
#' already appended so the result can be tracked against the shared national
#' master sample. Handles picking the correct bounding box/seed for the
#' chosen island and both the plain and \code{minRadius} variants of
#' \code{\link[spbal]{BAS}}'s output.
#'
#' @param sf_poly An \code{sf} object representing the area to sample within.
#'   Can be a single polygon or multiple polygons that will be unioned.
#' @param n Integer. Number of sample points to draw.
#' @param island Character. Which New Zealand island master sample to use.
#'   One of \code{"North"} or \code{"South"}.
#' @param min_radius Numeric or NULL. Minimum distance (in metres) between
#'   sample points, passed to \code{\link[spbal]{BAS}} as \code{minRadius}.
#'   Default \code{NULL} (no minimum distance constraint).
#' @param verbose Logical. Passed through to \code{\link[spbal]{BAS}}
#'   (default \code{FALSE}).
#'
#' @return An \code{sf} point data frame with columns \code{SiteID},
#'   \code{spbalSeqID} (or the minRadius equivalent), \code{uniqueID}, and
#'   \code{geometry}, in the CRS of the island master sample.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(spbalRef)
#'
#' sample_frame <- st_read("study_area.shp")
#'
#' pts <- masterSample(sample_frame, n = 50, island = "North")
#'
#' # With a minimum spacing between points
#' pts_spaced <- masterSample(sample_frame, n = 50, island = "North",
#'                            min_radius = 100)
#' }
#'
#' @export
#' @importFrom sf st_union st_as_sf st_transform st_crs st_bbox
#' @importFrom spbal BAS
masterSample <- function(sf_poly,
                         n,
                         island = c("North", "South"),
                         min_radius = NULL,
                         verbose = FALSE) {

  island <- match.arg(island)

  if(n <= 0) {
    stop("n must be positive")
  }

  bb_ref <- if(island == "North") {
    spbalRef::north_nz
  } else {
    spbalRef::south_nz
  }

  # Combine and transform to the master sample's CRS
  sf_poly <- st_union(sf_poly)
  sf_poly <- st_as_sf(sf::as_Spatial(sf_poly))
  sf_poly <- st_transform(sf_poly, st_crs(bb_ref$bb))

  # Check the polygon overlaps the chosen island's bounding box
  poly_bbox <- st_bbox(sf_poly)
  master_bbox <- st_bbox(bb_ref$bb)

  poly_overlaps <- !(poly_bbox["xmin"] > master_bbox["xmax"] ||
                       poly_bbox["xmax"] < master_bbox["xmin"] ||
                       poly_bbox["ymin"] > master_bbox["ymax"] ||
                       poly_bbox["ymax"] < master_bbox["ymin"])

  if(!poly_overlaps) {
    other_island <- if(island == "North") "South" else "North"
    stop(sprintf(
      "Sample frame does not overlap with %s Island bounding box. Check your polygon's coordinate system, or try island = '%s'.",
      island, other_island
    ))
  }

  bas_result <- spbal::BAS(
    shapefile = sf_poly,
    n = n,
    boundingbox = bb_ref$bb,
    seeds = bb_ref$seed,
    verbose = verbose,
    minRadius = min_radius
  )

  if(!is.null(min_radius)) {
    out <- bas_result$minRadius
    out$uniqueID <- paste0(bb_ref$mastersample_name, out$SiteID)
  } else {
    bas_result <- spbalRef::uniqueID(bb_ref, bas_result)
    out <- bas_result$sample
  }

  out
}
