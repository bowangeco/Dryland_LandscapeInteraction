data_dir <- "data/"

library(dplyr)
library(ggplot2)
library(parallel)
library(data.table)
library(terra)

cl <- makeCluster(12)
clusterExport(cl, c("data_dir"))
inpath <- paste0(data_dir, "GlobalLISMap1km/Results/3stl/")
file_list <- list.files(inpath, pattern = paste0("^Residual_.*\\.csv$"), full.names = TRUE)
Tile_list <- sub("^.*Residual_(.*)\\.csv$", "\\1", file_list) %>% as.numeric()
Tile_list <- Tile_list[order(Tile_list)]

multifun <- function(tile) {
    library(terra)
    library(data.table)
    library(dplyr)
    img_path <- paste0(data_dir, "GlobalLISMap1km/Results/3stl/Residual_", tile, ".csv")
    out_path <- paste0(data_dir, "GlobalLISMap1km/Results/4Neighborhood/Neigh_", tile)
    Residual_df <- fread(img_path) %>% as.data.frame()
    rast_id <- rast(paste0(data_dir, "GlobalLISMap1km/Results/1Clip/Rawdata/Tile_", tile, ".tif"))
    tif_list <- lapply(3:ncol(Residual_df), FUN = function(i) {
        xyz <- Residual_df[, c(1, 2, i)]
        rast_i <- rast(xyz, type = "xyz", crs = crs(rast_id))
        return(rast_i)
    })
    rast_tile <- rast(tif_list)
    Ext <- ext(rast_tile)
    Res <- res(rast_tile)
    Ext_crop <- ext(Ext[1] + Res[1], Ext[2] - Res[1], Ext[3] + Res[2], Ext[4] - Res[2])
    rast_tile_crop <- crop(rast_tile, Ext_crop)
    df <- terra::as.data.frame(rast_tile_crop, xy = TRUE, na.rm = FALSE)
    fwrite(df, file = paste0(out_path, ".csv"), row.names = FALSE)
    W_E <- matrix(data = c(NA, NA, NA, NA, NA, NA, NA, 1, NA), nrow = 3)
    rast_tile_E <- focal(rast_tile, w = W_E) %>% crop(., Ext_crop)
    df <- terra::as.data.frame(rast_tile_E, xy = TRUE, na.rm = FALSE)
    fwrite(df, file = paste0(out_path, "_E.csv"), row.names = FALSE)
    W_N <- matrix(data = c(NA, NA, NA, 1, NA, NA, NA, NA, NA), nrow = 3)
    rast_tile_N <- focal(rast_tile, w = W_N) %>% crop(., Ext_crop)
    df <- terra::as.data.frame(rast_tile_N, xy = TRUE, na.rm = FALSE)
    fwrite(df, file = paste0(out_path, "_N.csv"), row.names = FALSE)
    W_W <- matrix(data = c(NA, 1, NA, NA, NA, NA, NA, NA, NA), nrow = 3)
    rast_tile_W <- focal(rast_tile, w = W_W) %>% crop(., Ext_crop)
    df <- terra::as.data.frame(rast_tile_W, xy = TRUE, na.rm = FALSE)
    fwrite(df, file = paste0(out_path, "_W.csv"), row.names = FALSE)
    W_S <- matrix(data = c(NA, NA, NA, NA, NA, 1, NA, NA, NA), nrow = 3)
    rast_tile_S <- focal(rast_tile, w = W_S) %>% crop(., Ext_crop)
    df <- terra::as.data.frame(rast_tile_S, xy = TRUE, na.rm = FALSE)
    fwrite(df, file = paste0(out_path, "_S.csv"), row.names = FALSE)
}

x <- Tile_list
parLapplyLB(cl = cl, X = x, fun = multifun)
stopCluster(cl)
