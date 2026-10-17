data_dir <- "data/"

library(dplyr)
library(ggplot2)
library(parallel)
library(data.table)
library(terra)

FileList <- list.files(paste0(data_dir, "Image/Global_MODIS_NDVI_1km/QAmaskedStack_final/"), pattern = paste0("^.*NDVI.*\\.tif$"), 
    full.names = TRUE)
tile_ID_list <- sub(".*(h\\d*v\\d*)_.*\\.tif$", "\\1", FileList)
LIS_temp <- rast(paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/LIS_10km.tif"))

for (i in 1:length(FileList)) {
    Total_rast <- rast(FileList[i])
    Img_10km <- terra::project(Total_rast, LIS_temp, method = "average", align_only = TRUE)
    LIS_i <- crop(LIS_temp, Img_10km, extend = TRUE)
    Img_10km_masked <- mask(Img_10km, LIS_i)
    outpath <- paste0(data_dir, "Image/MODIS_NDVI_Projected_10km_Arid_CutHyper/")
    writeRaster(Img_10km_masked, paste0(outpath, "NDVI_QAmask_tile_", i, ".tif"))
}
