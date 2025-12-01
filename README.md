# spbalRef

A companion package to spbal for storing seeds and bounding boxes for existing master samples plus additional convenience functions for 2 stage sampling and transect generation.

## Two-Stage Spatial Sampling with makePSU

The `makePSU()` function creates optimized two-stage spatial sampling designs by dividing a study area into primary sampling units (PSUs) and selecting spatially-balanced sample points within each PSU.

### Quick Example

```r
library(spbalRef)
library(sf)

# Load your study area
study_area <- st_read("conservation_area.shp")

# Create PSUs and sample points
results <- makePSU(
  sf_poly = study_area,
  psu_area = 300,        # Target 300 hectare PSUs
  n_samples = 25,        # 25 sample points per PSU
  island = "North"       # Use North Island master sample
)

# View summary
print(results)

# Access components
psus <- results$primary_sample_units
sample_points <- results$secondary_sample_units

# Export for field work
st_write(psus, "PSUs.shp")
st_write(sample_points, "field_points.gpx", driver = "GPX")
```

### Key Features

- **Optimized PSU boundaries** - Iteratively minimizes area variation using balanced clustering
- **Spatial balance** - Uses GRTS for PSU selection and BAS for secondary sampling
- **Legacy site integration** - Prioritize PSUs containing existing monitoring sites
- **Flexible optimization** - Choose between CV, range, or SD metrics
- **NZ master sample integration** - Seamlessly works with North and South Island frameworks

### Advanced Usage

```r
# With legacy sites and selective sampling
legacy_sites <- st_read("existing_plots.shp")

results <- makePSU(
  sf_poly = study_area,
  psu_area = 500,
  n_samples = 20,
  island = "South",
  legacy_sites = legacy_sites,
  base_n = 10,           # Select 10 base PSUs
  over_n = 5,            # Plus 5 oversample PSUs
  min.dist = 5000        # Minimum 5km between PSU centroids
)
```

See `?makePSU` for complete documentation and additional examples.

---
