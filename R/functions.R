#Define functions for use in spbalRef package

#' uniqueID
#'
#'Append a mastersample uniqueID to a sample draw from spbal
#'
#' @param mastersample The existing mastersample to use
#' @param sample The sample drawn from spbal
#'
#' @return An sf point file from spbal with the uniqueID of an existing mastersample appended to it
#' @export
#'
#' @examples
#' #Append the uniqueID to a sample draw from the mastersample for the South Island, New Zealand
#' uniqueID(south_nz,BASsample)
#'
uniqueID <- function(mastersample,sample){
  sample$sample$uniqueID <-paste0(mastersample$mastersample_name,
                                  sample$sample$SiteID)
  return(sample)
}
