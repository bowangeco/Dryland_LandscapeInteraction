data_dir <- "data/"

library(dplyr)
library(ggplot2)
library(parallel)
library(data.table)
library(terra)
library(imputeTS)

inpath <- paste0(data_dir, "GlobalLISMap1km/Results/1Clip/Rawdata/")
file_list <- list.files(inpath, pattern = paste0("^Tile_.*\\.tif$"), full.names = TRUE)
id_list <- sub("^.*Tile_(.*)\\.tif$", "\\1", file_list) %>% as.numeric()
file_list <- file_list[order(id_list)]
out_path <- paste0(data_dir, "GlobalLISMap1km/Results/2Interpolation/")
file_list_finished <- list.files(out_path, pattern = paste0("^Tile_.*_.*\\.csv$"), full.names = TRUE)
id_list_finished <- sub("^.*Tile_.*_(.*)\\.csv$", "\\1", file_list_finished) %>% as.numeric()
id_list_finished <- id_list_finished[order(id_list_finished)]
undo <- id_list[!id_list %in% id_list_finished]
undo[-1][(undo[-1] - undo[-length(undo)]) > 1]
undo[-length(undo)][(undo[-1] - undo[-length(undo)]) > 1]

for (img_path in file_list) {
    id <- sub("^.*Tile_(.*)\\.tif$", "\\1", img_path) %>% as.numeric()
    P <- 0.2
    rast <- rast(img_path)
    valid_n <- sum(values(sum(is.na(rast)) < (529 * P)))
    if (valid_n > 0) {
        df <- terra::as.data.frame(rast, xy = TRUE, na.rm = FALSE)
        valid_row <- rowSums(is.na(df[, -1:-2])) < (529 * P)
        valid_df <- df[valid_row, ]
        valid_df_filled <- lapply(1:nrow(valid_df), FUN = function(i) {
            ts <- t(valid_df[i, -1:-2])
            ts_fill <- na_interpolation(ts)
            rdf <- cbind(valid_df[i, 1:2], t(ts_fill)) %>% data.frame()
            return(rdf)
        }) %>% do.call(rbind, .)
        df[valid_row, ] <- valid_df_filled
        fwrite(df, paste0(data_dir, "GlobalLISMap1km/Results/2Interpolation/Tile_interpolated_", id, 
            ".csv"))
    }
    else {
        NonValid_df <- data.frame(id = id, valid_n = 0)
        fwrite(NonValid_df, paste0(data_dir, "GlobalLISMap1km/Results/2Interpolation/Tile_NonValid_", 
            id, ".csv"))
    }
}

{
    out_path <- paste0(data_dir, "GlobalLISMap1km/Results/2Interpolation/")
    file_list <- list.files(out_path, pattern = paste0("^Tile_NonValid_.*\\.csv$"), full.names = TRUE)
    id_list <- sub("^.*Tile_NonValid_(.*)\\.csv$", "\\1", file_list) %>% as.numeric()
    file_list <- file_list[order(id_list)]
    Total_result <- lapply(file_list, FUN = function(p) {
        fread(p) %>% data.frame()
    }) %>% do.call(rbind, .)
    fwrite(Total_result, paste0(out_path, "Invalid_Pixels.csv"), row.names = FALSE)
    file.remove(file_list)
}
