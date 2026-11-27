## Shiny app for spbalRef: draw/upload a sample frame, generate a master
## sample or two-stage (PSU) sample, optionally seed PSUs with legacy sites,
## and optionally extend any point output into random transects.
##
## Launch with spbalRef::runSampleApp()

library(shiny)
library(leaflet)
library(leaflet.extras)
library(sf)

options(shiny.maxRequestSize = 50 * 1024^2)

# ---- helpers ----------------------------------------------------------

#' Read an uploaded spatial file (geojson/json/kml/gpkg, or a zipped shapefile)
read_uploaded_sf <- function(datapath, name) {
  ext <- tolower(tools::file_ext(name))

  if (ext == "zip") {
    tmp_dir <- tempfile()
    dir.create(tmp_dir)
    utils::unzip(datapath, exdir = tmp_dir)
    shp <- list.files(tmp_dir, pattern = "\\.shp$", full.names = TRUE, recursive = TRUE)
    if (length(shp) == 0) stop("No .shp file found inside the uploaded zip.")
    return(sf::st_read(shp[1], quiet = TRUE))
  }

  if (ext %in% c("geojson", "json", "kml", "gpkg")) {
    return(sf::st_read(datapath, quiet = TRUE))
  }

  stop(sprintf("Unsupported file type '.%s'. Use GeoJSON, KML, GPKG, or a zipped Shapefile.", ext))
}

#' Read an uploaded CSV of legacy points (needs lon/lat columns)
read_uploaded_csv_points <- function(datapath, id_col) {
  df <- utils::read.csv(datapath, stringsAsFactors = FALSE)

  lon_col <- intersect(c("lon", "Lon", "LON", "longitude", "x", "X"), names(df))[1]
  lat_col <- intersect(c("lat", "Lat", "LAT", "latitude", "y", "Y"), names(df))[1]

  if (is.na(lon_col) || is.na(lat_col)) {
    stop("CSV must contain longitude/latitude columns (e.g. 'lon'/'lat' or 'x'/'y').")
  }
  if (!(id_col %in% names(df))) {
    stop(sprintf("CSV does not contain the legacy ID column '%s'.", id_col))
  }

  sf::st_as_sf(df, coords = c(lon_col, lat_col), crs = 4326)
}

#' Convert a leaflet.extras draw-event feature (list parsed from GeoJSON) into
#' an sf polygon in EPSG:4326
draw_feature_to_sf <- function(feature) {
  coords <- feature$geometry$coordinates[[1]]
  ring <- do.call(rbind, lapply(coords, function(pt) c(pt[[1]], pt[[2]])))
  poly <- sf::st_polygon(list(ring))
  sf::st_sf(geometry = sf::st_sfc(poly, crs = 4326))
}

#' Zip a named list of sf objects up as GeoJSON files and return the zip path
zip_sf_results <- function(sf_list) {
  tmp_dir <- tempfile()
  dir.create(tmp_dir)
  files <- character(0)

  for (nm in names(sf_list)) {
    obj <- sf_list[[nm]]
    if (is.null(obj) || nrow(obj) == 0) next
    out_path <- file.path(tmp_dir, paste0(nm, ".geojson"))
    sf::st_write(obj, out_path, driver = "GeoJSON", quiet = TRUE, delete_dsn = TRUE)
    files <- c(files, out_path)
  }

  zip_path <- tempfile(fileext = ".zip")
  if (length(files) > 0) {
    zip::zip(zip_path, files = files, mode = "cherry-pick")
  }
  zip_path
}

# ---- UI ----------------------------------------------------------------

ui <- fluidPage(
  titlePanel("spbalRef: Sample Frame, Master Sample, PSUs & Transects"),

  sidebarLayout(
    sidebarPanel(
      width = 4,

      h4("1. Sample frame"),
      helpText("Draw a polygon on the map (use the polygon tool) or upload one below."),
      fileInput("frame_file", "Or upload a sample frame",
                accept = c(".geojson", ".json", ".kml", ".gpkg", ".zip")),
      actionButton("clear_frame", "Clear frame", icon = icon("eraser")),
      tags$hr(),

      radioButtons("island", "NZ Island (must contain the sample frame)",
                   choices = c("North" = "North", "South" = "South"),
                   selected = "North", inline = TRUE),
      tags$hr(),

      h4("2. Sample type"),
      radioButtons("mode", NULL,
                   choices = c("Master sample points only" = "master",
                               "Two-stage sample (PSUs)" = "psu"),
                   selected = "master"),

      conditionalPanel(
        condition = "input.mode == 'master'",
        numericInput("master_n", "Number of sample points", value = 50, min = 1),
        numericInput("master_min_radius", "Minimum spacing between points (m, optional)",
                     value = NA, min = 0)
      ),

      conditionalPanel(
        condition = "input.mode == 'psu'",
        numericInput("psu_area", "Target PSU area (ha)", value = 300, min = 0.1),
        numericInput("n_samples", "Secondary sample points per PSU", value = 25, min = 1),
        numericInput("psu_min_radius", "Minimum spacing between secondary points (m, optional)",
                     value = NA, min = 0),
        checkboxInput("quick_mode", "Quick mode (fast, single-pass clustering)", value = TRUE),
        numericInput("base_n", "PSUs to select as base sample (blank = all created)",
                     value = NA, min = 1),
        numericInput("over_n", "Oversample PSUs (blank = none)", value = NA, min = 0),
        numericInput("seed", "Random seed", value = 25945),
        tags$hr(),

        h4("Legacy sites (optional)"),
        helpText("PSUs containing an uploaded legacy site will be prioritised ",
                 "for base-sample selection by makePSU()."),
        fileInput("legacy_file", "Upload legacy sites (CSV with lon/lat, or GeoJSON/KML/GPKG/zipped Shapefile)",
                  accept = c(".csv", ".geojson", ".json", ".kml", ".gpkg", ".zip")),
        textInput("legacy_id_col", "Legacy ID column name", value = "PlotID")
      ),
      tags$hr(),

      h4("3. Transects (optional extension)"),
      checkboxInput("make_transects", "Generate transects from the sample points", value = FALSE),
      conditionalPanel(
        condition = "input.make_transects == true",
        numericInput("transect_length", "Transect length (m)", value = 100, min = 1),
        numericInput("transect_n_bearings", "Transects per point", value = 1, min = 1),
        radioButtons("transect_make", "Point is the transect's...",
                     choices = c("Start" = "start", "Middle" = "middle"),
                     selected = "start", inline = TRUE)
      ),
      tags$hr(),

      actionButton("generate", "Generate", icon = icon("play"), class = "btn-primary"),
      tags$hr(),
      downloadButton("download_results", "Download results (GeoJSON, zipped)")
    ),

    mainPanel(
      width = 8,
      leafletOutput("map", height = 700),
      verbatimTextOutput("status")
    )
  )
)

# ---- Server --------------------------------------------------------------

server <- function(input, output, session) {

  frame_sf <- reactiveVal(NULL)
  legacy_sf <- reactiveVal(NULL)
  results <- reactiveValues(primary = NULL, secondary = NULL,
                            master_points = NULL, transects = NULL)
  status_text <- reactiveVal("Draw or upload a sample frame, then click Generate.")

  output$map <- renderLeaflet({
    leaflet() %>%
      addTiles() %>%
      setView(lng = 172.5, lat = -41.5, zoom = 5) %>%
      addDrawToolbar(
        targetGroup = "drawn",
        polylineOptions = FALSE,
        circleOptions = FALSE,
        rectangleOptions = FALSE,
        markerOptions = FALSE,
        circleMarkerOptions = FALSE,
        editOptions = editToolbarOptions(edit = TRUE, remove = TRUE)
      )
  })

  # --- sample frame: drawn ---
  observeEvent(input$map_draw_new_feature, {
    feature <- input$map_draw_new_feature
    tryCatch({
      frame_sf(draw_feature_to_sf(feature))
      status_text("Sample frame set from drawn polygon.")
    }, error = function(e) {
      status_text(paste("Could not read drawn polygon:", conditionMessage(e)))
    })
  })

  observeEvent(input$map_draw_deleted_features, {
    frame_sf(NULL)
    status_text("Sample frame cleared.")
  })

  # --- sample frame: uploaded ---
  observeEvent(input$frame_file, {
    req(input$frame_file)
    tryCatch({
      f <- read_uploaded_sf(input$frame_file$datapath, input$frame_file$name)
      f <- sf::st_zm(f, drop = TRUE)
      frame_sf(f)
      status_text(sprintf("Sample frame uploaded (%d feature(s)).", nrow(f)))

      leafletProxy("map") %>%
        clearGroup("frame") %>%
        addPolygons(data = sf::st_transform(f, 4326), group = "frame",
                    color = "blue", weight = 2, fillOpacity = 0.1)
    }, error = function(e) {
      status_text(paste("Could not read uploaded sample frame:", conditionMessage(e)))
    })
  })

  observeEvent(input$clear_frame, {
    frame_sf(NULL)
    leafletProxy("map") %>% clearGroup("frame") %>% clearGroup("drawn")
    status_text("Sample frame cleared.")
  })

  # draw the frame whenever it's set (drawn polygons already render themselves
  # via the draw toolbar, but re-draw uploaded frames here)
  observeEvent(frame_sf(), {
    f <- frame_sf()
    if (is.null(f)) return(invisible())
  })

  # --- legacy sites: uploaded ---
  observeEvent(input$legacy_file, {
    req(input$legacy_file)
    tryCatch({
      ext <- tolower(tools::file_ext(input$legacy_file$name))
      pts <- if (ext == "csv") {
        read_uploaded_csv_points(input$legacy_file$datapath, input$legacy_id_col)
      } else {
        read_uploaded_sf(input$legacy_file$datapath, input$legacy_file$name)
      }
      legacy_sf(pts)
      status_text(sprintf("Legacy sites uploaded (%d site(s)).", nrow(pts)))

      leafletProxy("map") %>%
        clearGroup("legacy") %>%
        addCircleMarkers(data = sf::st_transform(pts, 4326), group = "legacy",
                          radius = 5, color = "darkorange", fillOpacity = 0.9)
    }, error = function(e) {
      status_text(paste("Could not read legacy sites:", conditionMessage(e)))
    })
  })

  # --- Generate ---
  observeEvent(input$generate, {

    f <- frame_sf()
    if (is.null(f)) {
      status_text("Please draw or upload a sample frame first.")
      return(invisible())
    }

    results$primary <- NULL
    results$secondary <- NULL
    results$master_points <- NULL
    results$transects <- NULL

    withProgress(message = "Generating sample...", value = 0.1, {

      out <- tryCatch({

        if (input$mode == "master") {

          incProgress(0.3, detail = "Drawing master sample")
          min_radius <- if (is.na(input$master_min_radius)) NULL else input$master_min_radius

          pts <- spbalRef::masterSample(
            sf_poly = f,
            n = input$master_n,
            island = input$island,
            min_radius = min_radius,
            verbose = FALSE
          )
          results$master_points <- pts
          list(points_for_transects = pts)

        } else {

          incProgress(0.2, detail = "Building PSUs")
          min_radius <- if (is.na(input$psu_min_radius)) NULL else input$psu_min_radius
          base_n <- if (is.na(input$base_n)) NULL else input$base_n
          over_n <- if (is.na(input$over_n)) NULL else input$over_n

          psu_result <- spbalRef::makePSU(
            sf_poly = f,
            psu_area = input$psu_area,
            n_samples = input$n_samples,
            min_radius = min_radius,
            island = input$island,
            base_n = base_n,
            over_n = over_n,
            legacy_sites = legacy_sf(),
            legacy_id_col = input$legacy_id_col,
            seed = input$seed,
            quick_mode = input$quick_mode,
            plot = FALSE,
            verbose = FALSE
          )

          results$primary <- psu_result$primary_sample_units
          results$secondary <- psu_result$secondary_sample_units
          list(points_for_transects = psu_result$secondary_sample_units)
        }

      }, error = function(e) {
        status_text(paste("Error:", conditionMessage(e)))
        NULL
      })

      if (is.null(out)) return(invisible())

      if (isTRUE(input$make_transects) && !is.null(out$points_for_transects) &&
          nrow(out$points_for_transects) > 0) {
        incProgress(0.3, detail = "Generating transects")
        tryCatch({
          results$transects <- spbalRef::makeTransects(
            transect.point = out$points_for_transects,
            line.length = input$transect_length,
            n.bearings = input$transect_n_bearings,
            transect.make = input$transect_make
          )
        }, error = function(e) {
          status_text(paste("Sample generated, but transects failed:", conditionMessage(e)))
        })
      }

      incProgress(0.2, detail = "Drawing map")
    })

    if (is.null(results$primary) && is.null(results$master_points)) return(invisible())

    status_text(sprintf(
      "Done. %s%s%s",
      if (!is.null(results$primary)) sprintf("%d PSUs, %d secondary points. ", nrow(results$primary), nrow(results$secondary)) else "",
      if (!is.null(results$master_points)) sprintf("%d master sample points. ", nrow(results$master_points)) else "",
      if (!is.null(results$transects)) sprintf("%d transects.", nrow(results$transects)) else ""
    ))

    proxy <- leafletProxy("map") %>%
      clearGroup("primary") %>% clearGroup("secondary") %>%
      clearGroup("master") %>% clearGroup("transects")

    if (!is.null(results$primary)) {
      pal_col <- ifelse(is.na(results$primary$siteuse), "lightgray",
                  ifelse(results$primary$siteuse == "Legacy", "orange",
                  ifelse(results$primary$siteuse == "Base", "lightgreen", "lightyellow")))
      proxy <- proxy %>%
        addPolygons(data = sf::st_transform(results$primary, 4326), group = "primary",
                    color = "gray40", weight = 1, fillColor = pal_col, fillOpacity = 0.5,
                    label = ~as.character(PSU_sample_order))
    }
    if (!is.null(results$secondary) && nrow(results$secondary) > 0) {
      proxy <- proxy %>%
        addCircleMarkers(data = sf::st_transform(results$secondary, 4326), group = "secondary",
                          radius = 3, color = "blue", fillOpacity = 0.8)
    }
    if (!is.null(results$master_points) && nrow(results$master_points) > 0) {
      proxy <- proxy %>%
        addCircleMarkers(data = sf::st_transform(results$master_points, 4326), group = "master",
                          radius = 4, color = "blue", fillOpacity = 0.8)
    }
    if (!is.null(results$transects) && nrow(results$transects) > 0) {
      proxy <- proxy %>%
        addPolylines(data = sf::st_transform(results$transects, 4326), group = "transects",
                      color = "red", weight = 2)
    }
  })

  output$status <- renderText(status_text())

  output$download_results <- downloadHandler(
    filename = function() paste0("spbalRef_results_", Sys.Date(), ".zip"),
    content = function(file) {
      sf_list <- list(
        sample_frame = frame_sf(),
        legacy_sites = legacy_sf(),
        primary_sample_units = results$primary,
        secondary_sample_units = results$secondary,
        master_sample_points = results$master_points,
        transects = results$transects
      )
      zip_path <- zip_sf_results(sf_list)
      file.copy(zip_path, file)
    },
    contentType = "application/zip"
  )
}

shinyApp(ui, server)
