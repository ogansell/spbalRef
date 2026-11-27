#' Generate Random Transects from Sample Points
#'
#' @description
#' Generates one or more randomly-bearinged transect lines per input point,
#' using the point as either the start or the midpoint of each transect.
#' Bearings are drawn uniformly at random and corrected for magnetic
#' declination (via \code{\link[wmm]{GetMagneticFieldWMM}}) so the returned
#' bearings are magnetic, not true, bearings.
#'
#' @param transect.point An \code{sf} (or \code{sp}) point object giving the
#'   locations to generate transects from. Coerced to \code{sf} if not
#'   already.
#' @param line.length Numeric. Length of each transect, in the units of
#'   \code{transect.point}'s CRS.
#' @param n.bearings Integer. Number of random-bearing transects to generate
#'   per point.
#' @param transect.make Character. Whether \code{transect.point} is the
#'   \code{"start"} of each transect or its \code{"middle"}. Partial matching
#'   is allowed; defaults to \code{"start"}.
#'
#' @return An \code{sf} object with one row per point per bearing, giving the
#'   transect geometry plus \code{bearing}, \code{bearing.number},
#'   \code{StartPointX}/\code{StartPointY} and \code{EndPointX}/\code{EndPointY}
#'   columns.
#'
#' @details
#' The installed \pkg{wmm} coefficient model only covers a roughly five-year
#' epoch. If the current date falls outside that epoch, magnetic declination
#' is instead calculated using the last date the installed model supports,
#' with a warning, rather than failing outright.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(spbalRef)
#'
#' pts <- st_as_sf(data.frame(lon = 172.5, lat = -41.2),
#'                  coords = c("lon", "lat"), crs = 4326)
#' pts <- st_transform(pts, 2193)
#'
#' transects <- makeTransects(pts, line.length = 100, n.bearings = 2,
#'                            transect.make = "start")
#' }
#'
#' @export
#' @importFrom sf st_as_sf st_coordinates st_transform st_touches st_geometry st_multilinestring st_sfc st_sf st_set_crs st_join st_buffer st_crs st_intersects st_contains as_Spatial
#' @importFrom lwgeom st_startpoint st_endpoint
#' @importFrom geosphere bearing
#' @importFrom igraph components graph_from_adj_list
#' @importFrom spatialEco sampleTransect
makeTransects <- function(transect.point,line.length,n.bearings,transect.make = c("start","middle")){

  transect.make <- match.arg(transect.make)

  # Updating to work for sf only. Start here...
  if (class(transect.point)[1] != "sf")
  {
    transect.point <-  st_as_sf(transect.point)
  }


  #Set up empty object
  dat <-NULL


  #Generate transects looping over the number of transects required
  for (x in 1:n.bearings) {

  name1 <-paste0("bearing",x)



  wgs84 <-"+proj=longlat +ellps=WGS84 +datum=WGS84"

  #Convert to wgs84 for calculating magnetic declination
  transect.point$EastingWGS <-st_coordinates(st_transform(transect.point,wgs84))[,1]
  transect.point$NorthingWGS <-st_coordinates(st_transform(transect.point,wgs84))[,2]

  #Calculate magnetic declination
  #wmm's installed coefficient model only covers a ~5yr epoch (currently through
  #2024-12-31) - fall back to that boundary date if Sys.time() is out of range,
  #rather than erroring, since the slow secular-variation drift makes a slightly
  #stale declination acceptable here.
  wmm.time <-tryCatch({
    wmm::GetMagneticFieldWMM(lon = 0, lat = 0, height = 0, time = Sys.time())
    Sys.time()
  }, error = function(e) {
    fallback <-as.POSIXct("2024-12-31")
    warning("wmm package does not support the current date; using ", fallback, " instead.")
    fallback
  })

  transect.point$magdec <-round(mapply(
    FUN = function(lon, lat) {
      wmm::GetMagneticFieldWMM(
        lon = lon,
        lat = lat,
        height = 0,
        time = wmm.time
      )$d
    },
    lon = transect.point$EastingWGS,
    lat = transect.point$NorthingWGS
  ),digits = 0)

  #Generate a bearing for the point
  azimuth <-sample(0:359,nrow(transect.point),replace = TRUE)

  #Calculate the bearing adjusted for magnetic declination

  azimuth_magdec <-(azimuth + transect.point$magdec) %% 360

  backbearing <-  ifelse(azimuth_magdec > 180,azimuth_magdec - 180,
                         azimuth_magdec + 180)

  transect.point$bearing.number <- name1

  #If the random point is the middle of the transect then::
  if(transect.make == "middle"){

  transect1 <-sampleTransect(
    x = st_as_sf(transect.point),
    distance = rep(line.length/2, nrow(transect.point)),
    azimuth = azimuth_magdec)

  transect2 <-sampleTransect(
    x = st_as_sf(transect.point),
    distance = rep(line.length/2, nrow(transect.point)),
    azimuth = backbearing)

  #Make 1 lines object to merge the 2 lines together
  lines <-rbind(transect1$transects, transect2$transects)

  #But if the transect is generated from the start point then:
  } else {
    transect1 <-sampleTransect(
      x = st_as_sf(transect.point),
      distance = rep(line.length, nrow(transect.point)),
      azimuth = azimuth_magdec)
  lines <-transect1$transects
  }

  #Strip the object back to only the geometry
  lines <-st_geometry(lines)
  #Identify where lines touch
  com <- components(graph_from_adj_list(st_touches(lines)))
  #Merge lines using tapply
  mlines <- do.call(st_sfc,
                    tapply(lines,
                           com$membership,
                           function(i){
                             st_multilinestring(lines[i])}))
  #Convert to sf
  mlines <-st_sf(mlines)
  #Set the crs
  mlines <-st_set_crs(mlines,st_crs(transect.point))
  #Join with point based data
  transect <-st_join(mlines, st_buffer(transect.point,20),
                     join = st_intersects)
  #Get start and end of each line point
  if(transect.make == "middle"){
  st.points <- st_startpoint(transect)}else{
  st.points <-transect.point[,"geometry"]}
  st.points <- st_join(st_as_sf(as_Spatial(st.points)),st_buffer(transect,10), join = st_contains)
  end.points <-st_endpoint(transect)
  end.points <- st_join(st_as_sf(as_Spatial(end.points)),st_buffer(transect,10), join = st_contains)
  #Get bearing for each line from the start point
  transect$bearing <-bearing(p1 = as_Spatial(st_transform(st.points,wgs84)),p2 = as_Spatial(st_transform(end.points,wgs84)))
  transect$bearing <-as.integer(ifelse( transect$bearing < 0, transect$bearing + 360, transect$bearing))
  transect$backbearing <-NULL
  transect$StartPointX <-st_coordinates(st.points)[,1]
  transect$StartPointY <-st_coordinates(st.points)[,2]
  transect$EndPointX <-st_coordinates(end.points)[,1]
  transect$EndPointY <-st_coordinates(end.points)[,2]
  transect$azimuth_magdec <-NULL
  transect$magdec <-NULL
  dat <-rbind(transect,dat)

}
  return(dat)
}
