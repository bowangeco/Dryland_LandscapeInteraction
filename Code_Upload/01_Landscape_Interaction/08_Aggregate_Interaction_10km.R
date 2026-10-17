data_dir <- "data/"

library(dplyr)
library(data.table)
library(terra)

clip_roi <- vect(paste0(data_dir, "Image/DrylandRegion/GlobalAridRegion_cutHyperArid.shp"))
for (folder in c("LIS_QA_Band/Reliability_Numb/", "10kmProject/")) {
    dir.create(paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/", folder), recursive = TRUE, showWarnings = FALSE)
}

{
    out_path <- paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/LIS_QA_Band/Reliability_Numb/")
    Imgpath <- paste0(data_dir, "Image/Global_MODIS_NDVI_1km/SelectedandStack/")
    FileList <- sort(list.files(Imgpath, pattern = "^.*_Reliability\\.tif$", full.names = TRUE))
    if (length(FileList) == 0L) stop("No Reliability TIFF inputs were found.")
    count_paths <- vapply(FileList, function(path) paste0(out_path, basename(path)), character(1))
    for (i in seq_along(FileList)) {
        rasti <- rast(FileList[i])
        ISValid_rast <- rasti %in% c(0, 1)
        Numb <- sum(ISValid_rast, na.rm = TRUE)
        writeRaster(Numb, count_paths[i], overwrite = TRUE)
    }

    count_tiles <- lapply(count_paths, rast)
    counts_native <- if (length(count_tiles) == 1L) count_tiles[[1]] else do.call(terra::merge, count_tiles)
    LIS_1km <- rast(paste0(data_dir, "GlobalLISMap1km/Results/Merged_LIS_map/LIS_1km.tif"))
    project(counts_native, LIS_1km, method = "near",
        filename = paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/LIS_QA_Band/LIS_Valid_Count_1km.tif"),
        overwrite = TRUE)
}

LIS_1km <- rast(paste0(data_dir, "GlobalLISMap1km/Results/Merged_LIS_map/LIS_1km.tif"))
Landcover <- rast(paste0(data_dir, "Image/Global_landcover/Landcove_ESA_1km.tif")) %>% project(., LIS_1km) %>% 
    crop(., LIS_1km)
LIS_1km[Landcover == 80] <- NA
LIS_10km_NonQAmasked <- project(LIS_1km, crs(LIS_1km), res = res(LIS_1km) * 10)
LIS_QA_Numb <- rast(paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/LIS_QA_Band/LIS_Valid_Count_1km.tif"))
QA_mask <- LIS_QA_Numb < (529 * 0.8)
LIS_1km_QAmask <- LIS_1km
LIS_1km_QAmask[QA_mask] <- NA
LIS_1km_QAmask_To10km <- project(LIS_1km_QAmask, crs(LIS_1km_QAmask), res = res(LIS_1km_QAmask) * 10, 
    method = "average")
QA_mask_10km <- project(QA_mask, crs(LIS_1km_QAmask), res = res(LIS_1km_QAmask) * 10, method = "average")
LIS_QA_band <- 1 - QA_mask_10km
LIS_QA_band <- mask(LIS_QA_band, LIS_10km_NonQAmasked)
clip_roi <- vect(paste0(data_dir, "Image/DrylandRegion/GlobalAridRegion_cutHyperArid.shp"))
LIS_QA_band_cliped <- crop(LIS_QA_band, clip_roi) %>% mask(., clip_roi)
writeRaster(LIS_QA_band_cliped, paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/LIS_Valid_Fraction_10km.tif"), 
    overwrite = TRUE)

LIS_1km_QAmask_To10km_QAmasked <- LIS_1km_QAmask_To10km
LIS_1km_QAmask_To10km_QAmasked[QA_mask_10km > 0.2] <- NA
LIS_1km_QAmask_To10km_QAmasked <- crop(LIS_1km_QAmask_To10km_QAmasked, clip_roi) %>%
    mask(., clip_roi)
writeRaster(LIS_1km_QAmask_To10km_QAmasked, paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/LIS_10km.tif"), overwrite = TRUE)
