#' Launch the Interactive Sampling Shiny App
#'
#' @description
#' Launches a Shiny app for drawing or uploading a sample frame, generating a
#' master sample (\code{\link{masterSample}}) or a two-stage PSU sample
#' (\code{\link{makePSU}}, optionally seeded with legacy sites), and
#' extending any of the resulting points into random transects
#' (\code{\link{makeTransects}}). Results can be downloaded as zipped
#' GeoJSON.
#'
#' @param ... Additional arguments passed to \code{\link[shiny]{runApp}}
#'   (e.g. \code{launch.browser}, \code{port}).
#'
#' @return Does not return; runs the Shiny app until stopped.
#'
#' @examples
#' \dontrun{
#' runSampleApp()
#' }
#'
#' @export
runSampleApp <- function(...) {
  required <- c("shiny", "leaflet", "leaflet.extras", "zip")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(sprintf(
      "The following packages are required to run the app but are not installed: %s",
      paste(missing, collapse = ", ")
    ))
  }

  app_dir <- system.file("shiny_app", package = "spbalRef")
  if (app_dir == "") {
    stop("Could not find the shiny_app directory. Try re-installing spbalRef.")
  }

  shiny::runApp(app_dir, ...)
}
