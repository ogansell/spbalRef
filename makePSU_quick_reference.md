# makePSU Quick Reference

## Basic Usage

```r
library(spbalRef)

results <- makePSU(
  sf_poly = study_area,      # Your study area polygon (sf object)
  psu_area = 300,            # Target PSU size in hectares
  n_samples = 25,            # Sample points per PSU
  island = "North"           # "North" or "South" for NZ master sample
)
```

## Essential Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `sf_poly` | Study area polygon (sf object) | *Required* |
| `psu_area` | Target PSU size (hectares) | *Required* |
| `n_samples` | Points per PSU | *Required* |
| `island` | "North" or "South" | *Required* |
| `base_n` | Number of base PSUs to select | All PSUs |
| `over_n` | Number of oversample PSUs | 0 |
| `optimization_metric` | "cv", "range", or "sd" | "cv" |

## Optimization Control

```r
results <- makePSU(
  sf_poly = study_area,
  psu_area = 500,
  n_samples = 20,
  island = "South",
  cluster_population = 200000,    # Total points for optimization
  optimization_start = 10000,     # Starting iteration
  optimization_step = 3000,       # Step size
  optimization_metric = "cv"      # Optimization metric
)
```

## Legacy Sites Integration

```r
results <- makePSU(
  sf_poly = study_area,
  psu_area = 400,
  n_samples = 20,
  island = "North",
  legacy_sites = existing_plots,   # sf object with existing sites
  legacy_id_col = "PlotID",        # ID column name
  base_n = 10,                     # Select 10 base PSUs
  over_n = 5                       # Plus 5 oversample
)
```

## GRTS Options

```r
results <- makePSU(
  sf_poly = study_area,
  psu_area = 300,
  n_samples = 15,
  island = "North",
  min.dist = 5000,              # Min distance between PSU centroids (m)
  grts_seltype = "equal",       # Selection type
  grts_caty_var = "habitat"     # Stratification variable
)
```

## Output Structure

```r
names(results)
# [1] "primary_sample_units"   - PSU polygons with selection order
# [2] "secondary_sample_units" - Sample points within PSUs
# [3] "optimization_summary"   - Optimization statistics
# [4] "optimization_data"      - Full optimization data
# [5] "optimal_n_points"       - Optimal iteration value
# [6] "pielou_evenness"        - Spatial balance metric
# [7] "psu_area_stats"         - PSU area statistics
# [8] "optimization_metric"    - Metric used
# [9] "settings"               - Function settings
```

## Accessing Results

```r
# Get PSUs
psus <- results$primary_sample_units

# Get sample points
points <- results$secondary_sample_units

# Get base sample only
base_psus <- psus %>% filter(siteuse == "Base")
base_points <- points %>% filter(siteuse == "Base")

# View statistics
print(results)
results$psu_area_stats
```

## Export Results

```r
# Shapefiles
st_write(results$primary_sample_units, "PSUs.shp")
st_write(results$secondary_sample_units, "points.shp")

# GPS format
st_write(base_points, "field_points.gpx", driver = "GPX")

# Field datasheet
field_data <- results$secondary_sample_units %>%
  st_set_geometry(NULL) %>%
  filter(siteuse == "Base") %>%
  select(PSU_sample_order, Transect_sample_order, uniqueID)

write.csv(field_data, "field_sheet.csv", row.names = FALSE)
```

## Common Recipes

### Minimize computation time
```r
results <- makePSU(
  sf_poly = study_area,
  psu_area = 1000,               # Fewer, larger PSUs
  cluster_population = 50000,    # Fewer points
  optimization_step = 5000,      # Larger steps
  n_samples = 10,
  island = "North"
)
```

### Maximize PSU uniformity
```r
results <- makePSU(
  sf_poly = study_area,
  psu_area = 300,
  cluster_population = 300000,   # More points
  optimization_step = 1000,      # Smaller steps
  optimization_metric = "cv",    # Use CV
  n_samples = 25,
  island = "North"
)
```

### Field sampling plan
```r
# Create PSUs
results <- makePSU(
  sf_poly = study_area,
  psu_area = 500,
  n_samples = 20,
  island = "North",
  base_n = 5,      # Visit 5 PSUs this year
  over_n = 3       # Keep 3 in reserve
)

# Get first PSU to visit
first_psu <- results$primary_sample_units %>%
  filter(PSU_sample_order == 1)

first_psu_points <- results$secondary_sample_units %>%
  filter(PSU_sample_order == 1) %>%
  arrange(Transect_sample_order)
```

## Troubleshooting

| Issue | Solution |
|-------|----------|
| Too few PSUs created | Decrease `psu_area` |
| High CV in PSU areas | Increase `cluster_population`, decrease `optimization_step` |
| Slow computation | Decrease `cluster_population`, increase `optimization_step` |
| PSUs not spatially balanced | Check `pielou_evenness` value, adjust `min.dist` |
| Not enough sample points | Increase `n_samples` or check PSU size |

## Tips

- **Minimum 4-5 PSUs** recommended for meaningful design
- **CV < 0.15** indicates good PSU uniformity
- **Pielou > 0.007** indicates good spatial balance for NZ
- Use **legacy sites** to integrate existing monitoring
- **base_n** controls immediate sampling effort
- **over_n** provides backup sites if access issues arise

## Performance Benchmarks

Approximate run times on typical desktop (4 cores, 16GB RAM):

| Area (ha) | PSUs | Population | Time |
|-----------|------|------------|------|
| 10,000 | 10 | 100,000 | ~2 min |
| 50,000 | 20 | 200,000 | ~8 min |
| 100,000 | 30 | 300,000 | ~20 min |

## See Also

- `?BAS` - Balanced acceptance sampling
- `?spsurvey::grts` - GRTS documentation
- `?sp_balance` - Spatial balance assessment
- Vignette: `vignette("makePSU", package = "spbalRef")`
