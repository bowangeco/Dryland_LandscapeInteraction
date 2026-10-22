data_dir <- "data/"

library(dplyr)
library(parallel)
library(data.table)

stl_fun <- function(row, freq = 23) {
    ts <- c(t(row)) %>% tail(., -2) %>% ts(., frequency = freq)
    if (anyNA(ts) == TRUE) {
        remainder <- rep(NA, length(ts))
    }
    else {
        s <- stl(ts, s.window = 23, t.window = 529, t.jump = 529, l.window = 23, robust = TRUE)
        remainder <- s$time.series[, "remainder"]
    }
    return_row <- c(row[c("x", "y")] %>% t(), remainder)
    return(return_row)
}

multifun <- function(id) {
    library(dplyr)
    library(data.table)
    df <- df_list[[id]]
    result_df <- apply(X = df, MARGIN = 1, FUN = function(x) {
        stl_fun(row = x)
    }) %>% t() %>% data.frame()
    fwrite(result_df, file = paste0(out_path, "OneTile_", tile, "_stl_", id, ".csv"), row.names = FALSE)
}

cores <- 20
inpath <- paste0(data_dir, "GlobalLISMap1km/Results/2Interpolation/")
file_list <- list.files(inpath, pattern = paste0("^Tile_interpolated_.*\\.csv"), full.names = TRUE)
Tile_list <- sub("^.*Tile_interpolated_(.*)\\.csv$", "\\1", file_list) %>% as.numeric()
Tile_list <- Tile_list[order(Tile_list)]
out_path <- paste0(data_dir, "GlobalLISMap1km/Results/3stl/")

for (tile in Tile_list) {
    img_path <- paste0(inpath, "Tile_interpolated_", tile, ".csv")
    cl <- makeCluster(cores)
    clusterExport(cl, c("data_dir"))
    NDVI_df <- fread(img_path)
    Complete_row <- complete.cases(NDVI_df)
    nr <- sum(Complete_row)
    if (nr < cores) {
        Complete_row[1:cores] <- TRUE
        nr <- sum(Complete_row)
    }
    NDVI_df_complete <- NDVI_df[Complete_row, ]
    n <- cores
    df_list <- split(NDVI_df_complete, rep(1:n, each = ceiling(nr/n), length.out = nr))
    clusterExport(cl, varlist = c("tile", "out_path", "df_list", "stl_fun", "multifun"))
    parLapplyLB(cl = cl, X = 1:n, fun = multifun)
    stopCluster(cl)
    {
        file_list_tile <- list.files(out_path, pattern = paste0("^OneTile_", tile, "_stl_"), full.names = TRUE)
        id <- sub(paste0(".*OneTile_", tile, "_stl_(.*)\\.csv$"), "\\1", file_list_tile) %>% as.numeric()
        file_list_tile <- file_list_tile[order(id)]
        Total_result <- lapply(file_list_tile, FUN = function(p) {
            fread(p) %>% data.frame()
        }) %>% do.call(rbind, .)
        file.remove(file_list_tile)
        Residual_df <- NDVI_df
        Residual_df[Complete_row, ] <- Total_result
        fwrite(Residual_df, paste0(out_path, "Residual_", tile, ".csv"))
    }
}
