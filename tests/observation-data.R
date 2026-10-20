library(terra)
library(sf)
source("01_Data_Prep/data_generation.R")

dem <- rast(nrows = 3, ncols = 3, xmin = 1000, xmax = 1300,
            ymin = 1000, ymax = 1300, crs = "EPSG:3857")
values(dem) <- 101:109
water <- dem
values(water) <- 0.5
values(water)[9] <- NA
permanent <- dem
values(permanent) <- NA
values(permanent)[c(1, 9)] <- 1
center <- xyFromCell(dem, 5)
dam <- st_as_sf(data.frame(x = center[1], y = center[2]),
                coords = c("x", "y"), crs = 3857)
dam <- st_transform(dam, 4326)
dam_path <- tempfile(fileext = ".gpkg")
st_write(dam, dam_path, quiet = TRUE)

observations <- build_observation_data(
  list(dem = dem, water_freq = water, perm_water = permanent),
  use_shore_elev = FALSE, dam_point_path = dam_path,
  dam_point_weight = 2, use_cache = FALSE
)
stopifnot(nrow(observations$elev_df) == 2)
stopifnot(all(abs(observations$elev_df$x - center[1]) < 1e-6))
stopifnot(all(abs(observations$elev_df$y - center[2]) < 1e-6))
stopifnot(all(observations$elev_df$elev == 105))
water_cells <- cellFromXY(dem, observations$water_df[, c("x", "y")])
stopifnot(length(water_cells) == 9, !anyDuplicated(water_cells))
stopifnot(all(observations$water_df$p[water_cells %in% c(1, 9)] == 1))
stopifnot(all(observations$water_df$p[!water_cells %in% c(1, 9)] == 0.5))
stopifnot(all(observations$water_df$y_water %in% c(50, 100)))
unlink(dam_path)
cat("Observation-data regression checks passed.\n")
