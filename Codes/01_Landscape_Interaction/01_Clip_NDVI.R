data_dir <- "data/"

library(dplyr)
library(ggplot2)
library(parallel)
library(data.table)
library(terra)

FileList <- list.files(paste0(data_dir, "Image/Global_MODIS_NDVI_1km/SelectedandStack/"), pattern = paste0("^.*NDVI.*\\.tif$"), 
    full.names = TRUE)
tile_ID_list <- sub(".*(h\\d*v\\d*)_.*\\.tif$", "\\1", FileList)

for (i in seq_along(FileList)) {
    Total_rast <- rast(FileList[i])
    Ext <- ext(Total_rast)
    Res <- res(Total_rast)
    {
        outpath <- paste0(data_dir, "GlobalLISMap1km/Results/1Clip/Rawdata/")
        {
            Width <- 200
            Height <- 200
            Grid_x <- c(seq(Ext[1] + Res[1], Ext[2] - Res[1], by = Width * Res[1]), Ext[2] - Res[1]) %>% 
                unique()
            Grid_y <- c(seq(Ext[3] + Res[2], Ext[4] - Res[1], by = Height * Res[2]), Ext[4] - Res[1]) %>% 
                unique()
            Ext_list <- lapply(1:(length(Grid_x) - 1), FUN = function(i) {
                xmin <- Grid_x[i]
                xmax <- Grid_x[i + 1]
                Ext_list_i <- lapply(1:(length(Grid_y) - 1), FUN = function(j) {
                  ymin <- Grid_y[j]
                  ymax <- Grid_y[j + 1]
                  Ext_ij <- ext(xmin - Res[1], xmax + Res[1], ymin - Res[2], ymax + Res[2])
                  return(Ext_ij)
                })
                return(Ext_list_i)
            }) %>% do.call(c, .)
        }
        n_df <- lapply(1:length(Ext_list), FUN = function(id) {
            Ext_id <- Ext_list[[id]]
            Split_rast_id <- crop(Total_rast, Ext_id)
            tile_id <- (i - 1) * length(Ext_list) + id
            writeRaster(Split_rast_id, paste0(outpath, "Tile_", tile_id, ".tif"))
            rdf <- data.frame(i = i, path = tile_ID_list[i], id = id, tile_id = tile_id)
            return(rdf)
        }) %>% do.call(rbind, .)
        fwrite(n_df, paste0(outpath, "tile_id_df", i, ".csv"))
    }
}

file_list <- list.files(outpath, pattern = paste0("^", "tile_id_df.*\\.csv"), full.names = TRUE)
tile_id <- sub(".*tile_id_df(\\d*)\\.csv$", "\\1", file_list) %>% as.numeric()
file_list <- file_list[order(tile_id)]
Total_result <- lapply(file_list, FUN = function(p) {
    fread(p) %>% data.frame()
}) %>% do.call(rbind, .)
fwrite(Total_result, paste0(outpath, "Total_tile_id_df.csv"), row.names = FALSE)

file.remove(file_list)
