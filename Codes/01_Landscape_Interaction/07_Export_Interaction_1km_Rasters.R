data_dir <- "data/"

library(data.table)
library(terra)

result_path <- paste0(data_dir, "GlobalLISMap1km/Results/6CCM/")
out_path <- paste0(data_dir, "GlobalLISMap1km/Results/7LIS_Directions/")
merged_path <- paste0(data_dir, "GlobalLISMap1km/Results/Merged_LIS_map/")
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
dir.create(merged_path, recursive = TRUE, showWarnings = FALSE)
files <- list.files(result_path, pattern = "^CCM_.*\\.csv$", full.names = FALSE)
Tile_list <- sort(as.numeric(sub("^CCM_(.*)\\.csv$", "\\1", files)))

for (tile in Tile_list) {
    img <- rast(paste0(data_dir, "GlobalLISMap1km/Results/1Clip/Rawdata/Tile_", tile, ".tif"))
    ccm <- as.data.frame(fread(paste0(result_path, "CCM_", tile, ".csv")))
    directions <- c("east", "north", "south", "west")
    for (direction in directions) {
        significant <- ccm[[paste0(direction, "_sig")]]
        ccm[[direction]][is.na(significant) | significant == 0] <- NA_real_
    }
    interaction <- rowMeans(ccm[, directions], na.rm = TRUE)
    interaction[!is.finite(interaction)] <- NA_real_
    xyz <- data.frame(x = ccm[[1]], y = ccm[[2]], value = interaction)
    image <- rast(xyz, type = "xyz", crs = crs(img), digits = 3)
    image <- project(image, "+proj=longlat", res = 0.008524211)
    writeRaster(image, paste0(out_path, "LIS_raster_", tile, ".tif"), overwrite = TRUE)
}

images <- lapply(Tile_list, function(tile) rast(paste0(out_path, "LIS_raster_", tile, ".tif")))
if (length(images) > 0L) {
    global <- if (length(images) == 1L) images[[1]] else do.call(terra::mosaic, images)
    writeRaster(global, paste0(merged_path, "LIS_1km.tif"), overwrite = TRUE)
}
