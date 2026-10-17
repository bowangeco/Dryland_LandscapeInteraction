data_dir <- "data/"

library(dplyr)
library(parallel)
library(data.table)
library(terra)

LIS <- rast(paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/LIS_10km.tif"))
LIS <- mask(LIS, LIS)
roi <- vect(paste0(data_dir, "Image/DrylandRegion/Fishnet.shp"))
roi_arid <- terra::as.polygons(LIS > 0)
roi <- intersect(roi_arid, roi)
result_path <- paste0(data_dir, "Coexistence_GlobalAridCutHyper_10km/Results/ClipedLIS/")

for (i in 1:length(roi)) {
    roi_i <- roi[i]
    LIS_i <- crop(LIS, roi_i, mask = TRUE)
    writeRaster(LIS_i, paste0(result_path, "LIS_clip_", i, ".tif"))
}
