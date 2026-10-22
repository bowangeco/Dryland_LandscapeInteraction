data_dir <- "data/"

library(data.table)
library(terra)

calc_lambda_ac1 <- function(x) {
    if (anyNA(x)) 
        return(NA_real_)
    ac1 <- cor(head(x, -1), tail(x, -1), use = "complete.obs")
    if (!is.finite(ac1) || ac1 <= 0) 
        return(NA_real_)
    log(ac1)
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
csv_dir <- paste0(data_dir, "CriticalSlowingDown_10km/Results/8_CSDofTotalTime/lambda_csv/")
raster_dir <- paste0(data_dir, "CriticalSlowingDown_10km/Results/8_CSDofTotalTime/lambda_tif/")
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
files <- list.files(residual_dir, pattern = "^Residual_.*\\.csv$", full.names = TRUE)

for (path in files) {
    tile <- sub("^Residual_(.*)\\.csv$", "\\1", basename(path))
    df <- as.data.frame(fread(path))
    valid <- complete.cases(df)
    result <- df[, 1:2, drop = FALSE]
    result$lambda_ac1 <- NA_real_
    if (any(valid)) {
        result$lambda_ac1[valid] <- apply(df[valid, -(1:2), drop = FALSE], 1, calc_lambda_ac1)
    }
    fwrite(result, file.path(csv_dir, paste0("Lambda_Tile_", tile, ".csv")))
    template <- rast(paste0(data_dir, "Image/MODIS_NDVI_Projected_10km_Arid_CutHyper/NDVI_QAmask_tile_", 
        tile, ".tif"))
    export_rasters(result, template, "lambda_ac1", raster_dir)
}

merge_rasters(raster_dir, "lambda_ac1")

library(dplyr)
library(ggplot2)
library(data.table)
library(terra)
library(ggspatial)
library(sf)

roi <- vect(paste0(data_dir, "Image/Boundary/Global_land.shp"))
clip_roi <- vect(paste0(data_dir, "Image/DrylandRegion/GlobalAridRegion_cutHyperArid.shp"))
out_path <- paste0(data_dir, "CriticalSlowingDown_10km/Results/8_CSDofTotalTime/lambda_tif/plot/")
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
plot_img <- rast(paste0(data_dir, "CriticalSlowingDown_10km/Results/8_CSDofTotalTime/lambda_tif/lambda_ac1_plotimg.tif"))
plot_img[plot_img < (-1)] <- -1
LIS_path <- paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/LIS_10km.tif")
LIS <- rast(LIS_path) %>% crop(., clip_roi) %>% mask(., clip_roi) %>% project(., plot_img, method = "near")
plot_img <- mask(plot_img, LIS)
mycolors <- c("#C2523B", "#D16631", "#E08221", "#ECA412", "#F2BF0C", "#FBE504", "#D3F700", "#78EB00", 
    "#1EDE00", "#0AC92D", "#17B368", "#1E9B89", "#1A808E", "#125382", "#001760")
label_breaks <- rev(c(-1, -0.75, -0.5, -0.25, 0) * (-1))

p <- ggplot() + geom_sf(data = st_as_sf(roi), fill = "grey92", color = NA) + layer_spatial(data = plot_img * 
    (-1)) + scale_fill_gradientn(name = "Recovery rate", limits = c(0, 1), oob = scales::squish, colors = mycolors, na.value = NA, values = scales::rescale((c(-1, 
    -0.89, -0.68, -0.53, -0.47, -0.37, -0.32, -0.32, -0.26, -0.21, -0.21, -0.16, -0.11, -0.11, 0) * -1) %>% 
    rev()), breaks = label_breaks, labels = c("0", "0.25", "0.5", "0.75", "≥1"), guide = guide_colorbar(ticks = TRUE, ticks.colour = "black", 
    ticks.linewidth = 0.4, frame.colour = "black", frame.linewidth = 0.2, nbin = 50, draw.ulim = TRUE, 
    draw.llim = TRUE, reverse = FALSE)) + labs(x = NULL, y = NULL) + scale_x_continuous(expand = c(0.03, 
    0.03), limits = c(-154, 168)) + scale_y_continuous(expand = expansion(mult = c(0.02, 0.05)), limits = c(-56, 
    75.5)) + theme_bw() + theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), text = element_text(family = "sans"), 
    legend.text = element_text(size = 10, hjust = 0.5, vjust = 3), legend.title = element_text(size = 12, 
        face = "plain", hjust = 0.5), legend.key.width = unit(1.2, "cm"), legend.key.height = unit(0.25, 
        "cm"), legend.background = element_rect(colour = "NA", fill = "NA", linewidth = 0), legend.direction = "horizontal", 
    legend.margin = margin(t = -12, unit = "pt"), strip.background = element_rect(colour = "NA", fill = "NA", 
        linewidth = 0), strip.text.x = element_text(size = 10, face = "bold"), strip.text.y = element_text(size = 10, 
        face = "bold"), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), 
    plot.margin = unit(c(0.03, 0.03, 0.03, 0.03), "cm"), panel.grid = element_blank(), panel.border = element_rect(fill = NA, 
        linewidth = 0.9), panel.spacing.x = unit(2, "mm"), panel.spacing.y = unit(1, "mm"), legend.position = c(0.572, 
        0.01), )
ggsave(plot = p, paste0(out_path, "Global_Recovery_Rate.jpeg"), width = 18, 
    height = 7.5, dpi = 900, units = "cm")
