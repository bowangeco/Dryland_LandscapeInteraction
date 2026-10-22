data_dir <- "data/"

library(dplyr)
library(ggplot2)
library(terra)
library(data.table)

result_path <- paste0(data_dir, "GlobalLISMap1km_Trend/Results/5_Project10km_trend/Trend/")

kendall_fun <- function(ccm_ts) {
    test_result <- cor.test(ccm_ts, 1:19, method = "kendall")
    p <- test_result$p.value
    tau <- test_result$estimate
    return(data.frame(p = p, tau = tau, row.names = NULL))
}

ValidPixels <- fread(paste0(data_dir, "GlobalLISMap1km/Results/Trend_Valid_Tiles.csv"))
x <- ValidPixels$tile

for (tile in x) {
    img <- rast(paste0(data_dir, "GlobalLISMap1km/Results/1Clip/Rawdata/Tile_", tile, ".tif"))
    years_df <- lapply(0:18, function(year) {
        LIS_raster_lonlat_10km <- rast(paste0(data_dir, "GlobalLISMap1km_Trend/Results/5_Project10km_trend/LIS_tif_10km/", 
            "LIS_raster_tile", tile, "Year", year, ".tif"))
        year_df <- terra::as.data.frame(LIS_raster_lonlat_10km, xy = TRUE, na.rm = FALSE)
        colnames(year_df) <- c("x", "y", "mean")
        if (year == 0) {
            year_rdf <- year_df[, c("x", "y", "mean")]
            colnames(year_rdf) <- c("x", "y", paste0("LIS", year))
        }
        else {
            year_rdf <- year_df[, c("mean")] %>% data.frame()
            colnames(year_rdf) <- paste0("LIS", year)
        }
        return(year_rdf)
    }) %>% do.call(cbind, .)
    df <- years_df
    valid_row <- (rowSums(!is.na(df[, -1:-2])) >= 17) & (apply(df[, -1:-2], 1, var, na.rm = TRUE) > 0)
    valid_df <- df[valid_row, ]
    if (nrow(valid_df) >= 1) {
        valid_results <- lapply(1:nrow(valid_df), function(rown) {
            ts <- valid_df[rown, -1:-2] %>% unlist() %>% unname()
            kendall_result <- kendall_fun(ts)
            rown_results <- cbind(valid_df[rown, 1:2], kendall_result)
            return(rown_results)
        }) %>% do.call(rbind, .)
    }
    else {
        valid_results <- NA
    }
    df_result <- df[, 1:2]
    df_result[, c("p", "tau")] <- NA
    df_result[valid_row, ] <- valid_results
    img <- rast(paste0(data_dir, "GlobalLISMap1km_Trend/Results/5_Project10km_trend/LIS_tif_10km/", "LIS_raster_tile", 
        tile, "Year", 1, ".tif"))
    crs <- crs(img)
    XYZ <- data.frame(lon = df_result[, 1], lat = df_result[, 2], value = df_result[, "tau"])
    tau_raster <- rast(XYZ, type = "xyz", crs = crs, digits = 3)
    writeRaster(tau_raster, filename = paste0(result_path, "tau/tau_", tile, ".tif"))
    XYZ <- data.frame(lon = df_result[, 1], lat = df_result[, 2], value = df_result[, "p"])
    p_raster <- rast(XYZ, type = "xyz", crs = crs, digits = 3)
    writeRaster(p_raster, filename = paste0(result_path, "p/p_", tile, ".tif"))
    tau_sign <- tau_raster
    tau_sign[p_raster >= 0.05] <- NA
    writeRaster(tau_sign, filename = paste0(result_path, "tau_sign/tau_sign_", tile, ".tif"))
}

{
    result_pathi <- paste0(data_dir, "GlobalLISMap1km_Trend/Results/5_Project10km_trend/LIS_tif_10km/")
    file_list_total <- list.files(result_pathi, pattern = paste0("^LIS_raster_tile.*Year.*\\.tif$"), 
        full.names = TRUE)
    tile_list <- sub(paste0("^.*LIS_raster_tile(.*)Year.*\\.tif$"), "\\1", file_list_total) %>% as.numeric()
    tile_list <- unique(tile_list)
    tile_list <- tile_list[order(tile_list)]
    for (tile in tile_list) {
        LIS_years <- lapply(0:18, function(year) {
            pathi <- paste0(data_dir, "GlobalLISMap1km_Trend/Results/5_Project10km_trend/LIS_tif_10km/", 
                "LIS_raster_tile", tile, "Year", year, ".tif")
            rasti <- rast(pathi)
            return(rasti)
        }) %>% rast()
        names(LIS_years) <- paste0("LIS_", 1:19)
        writeRaster(LIS_years, paste0(data_dir, "GlobalLISMap1km_Trend/Results/5_Project10km_trend/LIS_tif_10km_Years/tile_tif/", 
            "tile_", tile, ".tif"))
    }
    Total_img <- lapply(tile_list, function(tile) {
        img_i <- rast(paste0(data_dir, "GlobalLISMap1km_Trend/Results/5_Project10km_trend/LIS_tif_10km_Years/tile_tif/", 
            "tile_", tile, ".tif"))
        return(img_i)
    })
    batch_size <- 50
    batches <- split(Total_img, ceiling(seq_along(Total_img)/batch_size))
    merged_batches <- lapply(1:length(batches), function(i) {
        batch <- batches[[i]]
        Temp_merge <- do.call(terra::merge, batch)
        writeRaster(Temp_merge, paste0(data_dir, "GlobalLISMap1km_Trend/Results/5_Project10km_trend/LIS_tif_10km_Years/Merge_Temp/", 
            "Merge_Temp_", i, ".tif"))
    })
    batch_files <- list.files(paste0(data_dir, "GlobalLISMap1km_Trend/Results/5_Project10km_trend/LIS_tif_10km_Years/Merge_Temp/"), 
        pattern = "^Merge_Temp_.*\\.tif$", full.names = TRUE)
    stopifnot(length(batch_files) > 0L)
    r_base <- rast(batch_files[1])
    writeRaster(r_base, paste0(data_dir, "GlobalLISMap1km_Trend/Results/5_Project10km_trend/LIS_tif_10km_Years/LIS_Annual_Working_10km.tif"), 
        overwrite = TRUE)
    for (i in seq_along(batch_files)[-1]) {
        r_new <- rast(batch_files[i])
        r_out <- mosaic(r_base, r_new, filename = paste0(data_dir, "GlobalLISMap1km_Trend/Results/5_Project10km_trend/LIS_tif_10km_Years/LIS_Annual_10km.tif"), 
            overwrite = TRUE)
        writeRaster(r_out, paste0(data_dir, "GlobalLISMap1km_Trend/Results/5_Project10km_trend/LIS_tif_10km_Years/LIS_Annual_Working_10km.tif"), 
            overwrite = TRUE)
        r_base <- rast(paste0(data_dir, "GlobalLISMap1km_Trend/Results/5_Project10km_trend/LIS_tif_10km_Years/LIS_Annual_Working_10km.tif"))
    }
}

writeRaster(r_base, paste0(data_dir, "GlobalLISMap1km_Trend/Results/5_Project10km_trend/LIS_tif_10km_Years/LIS_Annual_10km.tif"), 
    overwrite = TRUE)

tau <- rast(paste0(data_dir, "GlobalLISMap1km_Trend/Results/5_Project10km_trend/Trend/tau.tif"))

{
    library(dplyr)
    library(ggplot2)
    library(data.table)
    library(terra)
    library(paletteer)
    library(ggspatial)
    library(sf)
    Total_rast_10km <- tau
    roi <- read_sf(paste0(data_dir, "Image/Boundary/Global_land.shp"))
    roi_fixed <- st_wrap_dateline(roi, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
    clip_roi <- vect(paste0(data_dir, "Image/DrylandRegion/GlobalAridRegion_cutHyperArid.shp"))

    plot_img <- Total_rast_10km %>% crop(clip_roi) %>% mask(., clip_roi)
    LIS_path <- paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/LIS_10km.tif")
    LIS <- rast(LIS_path) %>% crop(., clip_roi) %>% mask(., clip_roi) %>% project(., plot_img, method = "near")
    plot_img <- mask(plot_img, LIS)
    mycolors <- c("#5C0C46", "#9E2072", "white", "#65702E", "#3A4716")[5:1]
    p <- ggplot() + geom_sf(data = roi_fixed, fill = "grey92", color = NA) + layer_spatial(data = plot_img) + 
        scale_fill_gradientn(name = "Trend of landscape interaction strength", colors = mycolors, na.value = NA, limits = c(-1, 1), guide = guide_colorbar(ticks = TRUE, 
            ticks.colour = "black", ticks.linewidth = 0.4, frame.colour = "black", frame.linewidth = 0.2, 
            nbin = 50, draw.ulim = TRUE, draw.llim = TRUE)) + labs(x = NULL, y = NULL) + scale_x_continuous(expand = c(0.03, 
        0.03), limits = c(-154, 168)) + scale_y_continuous(expand = expansion(mult = c(0.02, 0.05)), 
        limits = c(-62.2, 71)) + theme_bw() + theme(plot.title = element_text(size = 16, hjust = 0.5, 
        face = "bold"), text = element_text(family = "sans"), legend.text = element_text(size = 10, hjust = 0.5, 
        vjust = 3), legend.title = element_text(size = 12, face = "plain", hjust = 0.5), legend.key.width = unit(1.2, 
        "cm"), legend.key.height = unit(0.25, "cm"), legend.background = element_rect(colour = "NA", 
        fill = "NA", linewidth = 0), legend.direction = "horizontal", legend.margin = margin(t = -12, 
        unit = "pt"), strip.background = element_rect(colour = "NA", fill = "NA", linewidth = 0), strip.text.x = element_text(size = 10, 
        face = "bold"), strip.text.y = element_text(size = 10, face = "bold"), axis.text = element_blank(), 
        axis.ticks = element_blank(), axis.title = element_blank(), plot.margin = unit(c(0.03, 0.03, 
            0.03, 0.03), "cm"), panel.grid = element_blank(), panel.border = element_rect(fill = NA, 
            linewidth = 0.9), panel.spacing.x = unit(2, "mm"), panel.spacing.y = unit(1, "mm"), legend.position = c(0.572, 
            0.01), )
    ggsave(plot = p, paste0(data_dir, "GlobalLISMap1km_Trend/Results/", "5_Project10km_trend/Trend/Global_LIS_Trend.png"), 
        width = 18, height = 8, dpi = 900, units = "cm")
}
