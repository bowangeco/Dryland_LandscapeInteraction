data_dir <- "data/"

library(terra)
library(dplyr)
library(data.table)
library(raster)

temp_dir <- tempdir()

if (!dir.exists(temp_dir)) dir.create(temp_dir)

Sys.setenv(TEMP = temp_dir, TMPDIR = temp_dir)
rasterOptions(tmpdir = tempdir())
terraOptions(tempdir = tempdir())
Imgpath <- paste0(data_dir, "Coexistence_GlobalAridCutHyper_10km/Results/ClipedLIS/")
FileList <- list.files(Imgpath, pattern = paste0("^LIS_clip_.*\\.tif$"), full.names = TRUE)
ID_list <- sub("^.*LIS_clip_(.*)\\.tif$", "\\1", FileList) %>% as.numeric()
ID_list <- ID_list[order(ID_list)]
Sentinel_Ext_df <- fread(paste0(data_dir, "Coexistence_GlobalAridCutHyper_10km/Results/1_Preprocessing/Ext_df.csv")) %>% 
    data.frame()
Sentinel_Ext_list <- lapply(1:nrow(Sentinel_Ext_df), FUN = function(i) {
    Ext_i <- ext(Sentinel_Ext_df[i, "xmin"], Sentinel_Ext_df[i, "xmax"], Sentinel_Ext_df[i, "ymin"], 
        Sentinel_Ext_df[i, "ymax"])
    return(Ext_i)
})
result_path <- paste0(data_dir, "Coexistence_GlobalAridCutHyper_10km/Results/1_Preprocessing/Rdata/")

for (ID in ID_list) {
    LIS_grid <- rast(paste0(data_dir, "Coexistence_GlobalAridCutHyper_10km/Results/ClipedLIS/", "LIS_clip_", 
        ID, ".tif"))
    if (sum(!is.na(values(LIS_grid))) == 0) {
        write.csv(NA, paste0(result_path, "NDVIValuesL_ID", ID, "Nodata.csv"))
    }
    else {
        Ext_LIS <- ext(LIS_grid)
        Overlap_Sentinel_IDs <- lapply(1:length(Sentinel_Ext_list), FUN = function(i) {
            Ext_i <- Sentinel_Ext_list[[i]]
            Overlap <- !is.null(terra::intersect(Ext_LIS, Ext_i))
        }) %>% do.call(c, .)
        Overlap_Sentinel_imgs <- lapply(Sentinel_Ext_df$id[Overlap_Sentinel_IDs], FUN = function(i) {
            Sentinel_i <- rast(paste0(data_dir, "Image/GlobalArid_SentinelNDVI/Sentinel_NDVI_masked_global_mean_", 
                i, ".tif"))
            return(Sentinel_i)
        })
        if (length(Overlap_Sentinel_imgs) == 0) {
            write.csv(NA, paste0(result_path, "NDVIValuesL_ID", ID, "Nodata.csv"))
        }
        else {
            if (length(Overlap_Sentinel_imgs) > 1) {
                Overlap_Sentinel_img <- do.call(mosaic, Overlap_Sentinel_imgs)
            }
            else {
                Overlap_Sentinel_img <- Overlap_Sentinel_imgs[[1]]
            }
            writeRaster(Overlap_Sentinel_img, filename = paste0(data_dir, "Coexistence_GlobalAridCutHyper_10km/Results/1_Preprocessing/Temp/", 
                "MergedSentinelNDVI_", ID, ".tif"), overwrite = TRUE)
            Overlap_Sentinel_img_crop <- crop(Overlap_Sentinel_img, Ext_LIS, filename = paste0(data_dir, 
                "Coexistence_GlobalAridCutHyper_10km/Results/1_Preprocessing/Temp/", "ReClipedSentinelNDVI_", 
                ID, ".tif"), overwrite = TRUE)
            rm(Overlap_Sentinel_imgs, Overlap_Sentinel_img)
            {
                NDVI_masked <- raster(Overlap_Sentinel_img_crop)
                LIS_grid <- raster(LIS_grid)
                grid_df <- rasterToPoints(LIS_grid) %>% data.frame()
                colnames(grid_df) <- c("x", "y", "LIS")
                resy <- yres(LIS_grid)
                resx <- xres(LIS_grid)
                r_row_list <- lapply(1:nrow(grid_df), FUN = function(i) {
                  row <- grid_df[i, ]
                  lon <- row["x"] %>% unlist() %>% unname()
                  lat <- row["y"] %>% unlist() %>% unname()
                  LIS <- row["LIS"] %>% unlist() %>% unname()
                  Bin_extent <- extent(lon - resx/2, lon + resx/2, lat - resy/2, lat + resy/2)
                  if (!is.null(terra::intersect(Bin_extent, NDVI_masked))) {
                    Bin_img <- crop(NDVI_masked, Bin_extent)
                    NDVI_values <- getValues(Bin_img)
                    NDVI_values_narm <- NDVI_values[!is.na(NDVI_values)]
                    if (length(NDVI_values_narm) > 10000) {
                      NDVI_values_narm <- sample(NDVI_values_narm, 10000)
                    }
                  }
                  else {
                    NDVI_values_narm <- NA
                  }
                  r_row <- c(lon, lat, LIS, NDVI_values_narm)
                  return(r_row)
                })
                df_list <- r_row_list
                if (length(df_list) > 10000) {
                  df_list_split <- split(df_list, rep(1:ceiling(length(df_list)/10000), each = 10000, 
                    length.out = length(df_list)))
                  lapply(1:length(df_list_split), FUN = function(i) {
                    df_list_i <- df_list_split[[i]]
                    saveRDS(df_list_i, file = paste0(result_path, "NDVIValuesListSplit_ID", ID, "_", 
                      i, "Nrow", length(df_list_i), ".rds"))
                  })
                }
                else {
                  saveRDS(df_list, file = paste0(result_path, "NDVIValuesListSplit_ID", ID, "_1", "Nrow", 
                    length(df_list), ".rds"))
                }
            }
            file.remove(paste0(data_dir, "Coexistence_GlobalAridCutHyper_10km/Results/1_Preprocessing/Temp/", 
                "MergedSentinelNDVI_", ID, ".tif"))
        }
    }
}
