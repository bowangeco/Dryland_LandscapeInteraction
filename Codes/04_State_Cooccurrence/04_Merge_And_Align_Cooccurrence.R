data_dir <- "data/"

library(dplyr)
library(data.table)
library(terra)

Imgpath <- paste0(data_dir, "Coexistence_GlobalAridCutHyper_10km/Results/ClipedLIS/")
FileList <- list.files(Imgpath, pattern = paste0("^LIS_clip_.*\\.tif$"), full.names = TRUE)
ID_list <- sub("^.*LIS_clip_(.*)\\.tif$", "\\1", FileList) %>% as.numeric()
ID_list <- ID_list[order(ID_list)]

{
    Temp_path <- paste0(data_dir, "Coexistence_GlobalAridCutHyper_10km/Results/2_Bistable_test_Temp/")
    out_path <- paste0(data_dir, "Coexistence_GlobalAridCutHyper_10km/Results/2_Bistable_test/")
    for (ID in ID_list) {
        Temp_files <- list.files(Temp_path, pattern = paste0("^Coexistence_ID", ID, "_id.*_10km\\.csv$"), 
            full.names = TRUE)
        if (length(Temp_files) != 0) {
            ids <- sub(paste0("^.*Coexistence_ID", ID, "_id(.*)_10km\\.csv$"), "\\1", Temp_files) %>% 
                as.numeric()
            Temp_files <- Temp_files[order(ids)]
            Total_df <- lapply(Temp_files, FUN = function(i) {
                csv_i <- fread(i) %>% as.data.frame()
                return(csv_i)
            }) %>% do.call(rbind, .)
            fwrite(Total_df, paste0(out_path, "Coexistence_ID", ID, "_10km.csv"))
        }
    }
}

LIS <- rast(paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/LIS_10km.tif"))
Coex_path <- paste0(data_dir, "Coexistence_GlobalAridCutHyper_10km/Results/4_Coexistence_plot_Bymodmean/", 
    "Global_Cooccurrence.tif")
Coex <- rast(Coex_path)
Coex_10km <- project(Coex, LIS, method = "near")
writeRaster(Coex_10km, paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/Cooccurrence_10km.tif"))
