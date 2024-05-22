library(sf)


#Define CRS projections

nztm <-"+proj=tmerc +lat_0=0 +lon_0=173 +k=0.9996 +x_0=1600000 +y_0=10000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
wgs84 <-"+proj=longlat +ellps=WGS84 +datum=WGS84"
epsg3005 <- "+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"


# Define bounding boxes

#Add coordinates for a mastersample bounding box

south_nz_bb <-data.frame(x = c(1089354,1721164 ),
                      y = c(4747979,5516919))

north_nz_bb <- data.frame(x = c(1510593, 2092000),
                       y = c(5390569,6223164))

auckland_islands_nz_bb <- data.frame(x = c(3490761,3525662),
                 y = c(4355830,4405196))

bc_ca_bb <- data.frame(x = c(85148,1280999),
                      y = c(33745,1351981))

bb.df <- c("xmin" = 85148, "ymin" = 33745, "xmax" = 1280999, "ymax" = 1351981)


#Convert bounding box to an sf object
south_nz_bb <-south_nz_bb %>%
  st_as_sf(coords = c("x","y"), crs = nztm) %>%
  st_bbox() %>% st_as_sfc()

north_nz_bb <-north_nz_bb %>%
  st_as_sf(coords = c("x","y"), crs = nztm) %>%
  st_bbox() %>% st_as_sfc()

auckland_islands_nz_bb <-auckland_islands_nz_bb %>%
  st_as_sf(coords = c("x","y"), crs = nztm) %>%
  st_bbox() %>% st_as_sfc()

bc_ca_bb <-bc_ca_bb %>%
  st_as_sf(coords = c("x","y"), crs = epsg3005) %>%
  st_bbox() %>% st_as_sfc()

# Add mastersample seeds

south_nz_seed <- c(4887260, 18041662)
north_nz_seed <- c(5137598, 8906854)
auckland_islands_nz_seed <- c(1925006, 6242772)
bc_ca_seed <- c(37916, 85846)

# Make a list of the objects and save for use

south_nz <-list(bb = south_nz_bb, seed = south_nz_seed,
                mastersample_name = "NZ_South")
usethis::use_data(south_nz,overwrite = TRUE)

north_nz <-list(bb = north_nz_bb, seed = north_nz_seed,
                mastersample_name = "NZ_North")
usethis::use_data(north_nz,overwrite = TRUE)

auckland_islands_nz <-list(bb = auckland_islands_nz_bb,
                           seed = auckland_islands_nz_seed,
                           mastersample_name = "NZ_AucklandIslands")
usethis::use_data(auckland_islands_nz,overwrite = TRUE)

bc_ca <-list(bb = bc_ca_bb,
             seed = bc_ca_seed,
             mastersample_name = "BC_Canada")
usethis::use_data(bc_ca,overwrite = TRUE)



