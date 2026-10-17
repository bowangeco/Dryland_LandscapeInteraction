data_dir <- "data/"

library(data.table)
library(terra)

calc_lambda_ac1 <- function(x) {
    if (anyNA(x)) 
        return(NA_real_)
    ac1 <- as.numeric(acf(x, lag.max = 1, plot = FALSE)$acf[2])
    if (!is.finite(ac1) || ac1 <= 0) 
        return(NA_real_)
    log(ac1)
}

kendall_trend <- function(x) {
    if (sum(!is.na(x)) <= 6L || var(x, na.rm = TRUE) == 0) 
        return(c(p_ac1 = NA_real_, tau_ac1 = NA_real_))
    test <- cor.test(x, seq_along(x), method = "kendall")
    c(p_ac1 = test$p.value, tau_ac1 = unname(test$estimate))
}

export_rasters <- function(df, template, columns, root) {
    for (column in columns) {
        folder <- file.path(root, column)
        dir.create(folder, recursive = TRUE, showWarnings = FALSE)
        xyz <- data.frame(x = df[[1]], y = df[[2]], value = df[[column]])
        image <- rast(xyz, type = "xyz", crs = crs(template), digits = 3)
        image <- project(image, "+proj=longlat", res = 0.08524211, method = "near")
        writeRaster(image, file.path(folder, paste0(column, "_raster_", tile, ".tif")), overwrite = TRUE)
    }
}

merge_rasters <- function(root, columns) {
    for (column in columns) {
        paths <- list.files(file.path(root, column), pattern = "\\.tif$", full.names = TRUE)
        if (length(paths) == 0L) 
            next
        images <- lapply(paths, rast)
        merged <- if (length(images) == 1L) 
            images[[1]]
        else do.call(terra::merge, images)
        writeRaster(merged, file.path(root, paste0(column, ".tif")), overwrite = TRUE)
    }
}

residual_dir <- paste0(data_dir, "CriticalSlowingDown_10km/Results/3stl/")
csv_dir <- paste0(data_dir, "CriticalSlowingDown_10km/Results/lambda_residual/")
trend_csv_dir <- paste0(data_dir, "CriticalSlowingDown_10km/Results/lambda_residual_trend/")
raster_dir <- paste0(data_dir, "CriticalSlowingDown_10km/Results/7RecoverRateTif/residual/")
trend_raster_dir <- paste0(data_dir, "CriticalSlowingDown_10km/Results/7RecoverRateTif/residual_kendall/")
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(trend_csv_dir, recursive = TRUE, showWarnings = FALSE)
files <- list.files(residual_dir, pattern = "^Residual_.*\\.csv$", full.names = TRUE)

for (path in files) {
    tile <- sub("^Residual_(.*)\\.csv$", "\\1", basename(path))
    df <- as.data.frame(fread(path))
    stopifnot(ncol(df) == 531L)
    valid <- complete.cases(df)
    columns <- paste0("lambda_ac1_", 1:19)
    result <- df[, 1:2, drop = FALSE]
    result[columns] <- NA_real_
    if (any(valid)) {
        for (year in 1:19) {
            window <- (1:115) + 23L * (year - 1L)
            result[valid, columns[year]] <- apply(df[valid, window + 2L, drop = FALSE], 1, calc_lambda_ac1)
        }
    }
    fwrite(result, file.path(csv_dir, paste0("Lambda_Window_Tile_", tile, ".csv")))
    template <- rast(paste0(data_dir, "Image/MODIS_NDVI_Projected_10km_Arid_CutHyper/NDVI_QAmask_tile_", 
        tile, ".tif"))
    images <- lapply(columns, function(column) {
        xyz <- data.frame(x = result[[1]], y = result[[2]], value = result[[column]])
        image <- rast(xyz, type = "xyz", crs = crs(template), digits = 3)
        project(image, "+proj=longlat", res = 0.08524211, method = "near")
    })
    folder <- file.path(raster_dir, "lambda_ac1")
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    stack <- rast(images)
    names(stack) <- columns
    writeRaster(stack, file.path(folder, paste0("lambda_ac1_raster_", tile, ".tif")), overwrite = TRUE)
    trend <- result[, 1:2, drop = FALSE]
    values <- t(apply(result[, columns, drop = FALSE], 1, kendall_trend))
    trend <- cbind(trend, as.data.frame(values))
    fwrite(trend, file.path(trend_csv_dir, paste0("Lambda_Trend_Tile_", tile, ".csv")))
    export_rasters(trend, template, c("p_ac1", "tau_ac1"), trend_raster_dir)
}

merge_rasters(trend_raster_dir, c("p_ac1", "tau_ac1"))

library(dplyr)
library(ggplot2)
library(data.table)
library(terra)
library(ggspatial)
library(sf)

out_path <- paste0(data_dir, "CriticalSlowingDown_10km/Results/7RecoverRateTif/residual_kendall/")
Total_rast_10km <- rast(paste0(out_path, "tau_ac1.tif"))
roi <- read_sf(paste0(data_dir, "Image/Boundary/Global_land.shp"))
roi_fixed <- st_wrap_dateline(roi, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
clip_roi <- vect(paste0(data_dir, "Image/DrylandRegion/GlobalAridRegion_cutHyperArid.shp"))

plot_img <- Total_rast_10km %>% crop(clip_roi) %>% mask(., clip_roi)
LIS_path <- paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/LIS_10km.tif")
LIS <- rast(LIS_path) %>% crop(., clip_roi) %>% mask(., clip_roi) %>% project(., plot_img, method = "near")
plot_img <- mask(plot_img, LIS)
mycolors <- c("#223D6C", "#4075BA", "white", "#DF2011", "#B20000")

{
    p <- ggplot() + geom_sf(data = roi_fixed, fill = "grey92", color = NA) + layer_spatial(data = plot_img * 
        (-1)) + scale_fill_gradientn(name = "Trend of recovery rate", colors = mycolors %>% rev(), na.value = NA, limits = c(-1, 
        1), guide = guide_colorbar(ticks = TRUE, ticks.colour = "black", ticks.linewidth = 0.4, frame.colour = "black", 
        frame.linewidth = 0.2, nbin = 50, draw.ulim = FALSE, draw.llim = FALSE)) + labs(x = NULL, y = NULL) + 
        scale_x_continuous(expand = c(0.03, 0.03), limits = c(-154, 168)) + scale_y_continuous(expand = expansion(mult = c(0.02, 
        0.05)), limits = c(-62.2, 71)) + theme_bw() + theme(plot.title = element_text(size = 16, hjust = 0.5, 
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
    ggsave(plot = p, paste0(data_dir, "CriticalSlowingDown_10km/Results/", "7RecoverRateTif/residual_kendall/Global_Recovery_Rate_Trend.png"), 
        width = 18, height = 8, dpi = 900, units = "cm")
}
