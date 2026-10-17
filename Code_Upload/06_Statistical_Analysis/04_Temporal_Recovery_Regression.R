data_dir <- "data/"
output_dir <- "output/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

library(dplyr)
library(ggplot2)
library(parallel)
library(data.table)
library(terra)
library(ggpointdensity)

clip_roi <- vect(paste0(data_dir, "Image/DrylandRegion/GlobalAridRegion_cutHyperArid.shp"))
LIS_path <- paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/LIS_10km.tif")
LIS <- rast(LIS_path) %>% crop(., clip_roi) %>% mask(., clip_roi)
LIS_trend <- rast(paste0(data_dir, "GlobalLISMap1km_Trend/Results/5_Project10km_trend/Trend/tau.tif"))
LIS_trend <- project(LIS_trend, LIS, method = "near") %>% crop(., LIS) %>% mask(., LIS)
writeRaster(LIS_trend, paste0(output_dir, "LIS_Kendall_Tau_10km.tif"), overwrite = TRUE)

LIS_trend_p <- rast(paste0(data_dir, "GlobalLISMap1km_Trend/Results/5_Project10km_trend/Trend/p.tif"))
LIS_trend_p <- project(LIS_trend_p, LIS, method = "near") %>% crop(., LIS) %>% mask(., LIS)
writeRaster(LIS_trend_p, paste0(output_dir, "LIS_Kendall_P_10km.tif"), overwrite = TRUE)

lambda_ac1_trend <- rast(paste0(data_dir, "CriticalSlowingDown_10km/Results/7RecoverRateTif/residual_kendall/tau_ac1.tif"))
lambda_ac1_trend_p <- rast(paste0(data_dir, "CriticalSlowingDown_10km/Results/7RecoverRateTif/residual_kendall/p_ac1.tif"))
lambda_ac1_trend <- project(lambda_ac1_trend, LIS, method = "near") %>% crop(., LIS) %>% mask(., LIS)
lambda_ac1_trend_p <- project(lambda_ac1_trend_p, LIS, method = "near") %>% crop(., LIS) %>% mask(., 
    LIS)
writeRaster(lambda_ac1_trend * (-1), paste0(output_dir, "Recovery_Kendall_Tau_10km.tif"), overwrite = TRUE)
writeRaster(lambda_ac1_trend_p, paste0(output_dir, "Recovery_Kendall_P_10km.tif"), overwrite = TRUE)

imgs_df <- as.data.frame(rast(list(LIS_trend, lambda_ac1_trend)))
colnames(imgs_df) <- c("LIS_trend", "lambda_ac1_trend")
outpath <- paste0(data_dir, "Summary_AridCutHyper/Results/4_RecoverRate/PointPlot/")
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)

library(lmodel2)
library(dplyr)

rma_df <- imgs_df %>% transmute(x = LIS_trend, y = lambda_ac1_trend) %>% filter(is.finite(x), is.finite(y))
rma_fit <- lmodel2(y ~ x, data = rma_df, range.x = "interval", range.y = "interval", nperm = 0)

rma_result <- rma_fit$regression.results %>% filter(Method == "RMA")
rma_slope <- -rma_result$Slope[1]
rma_intercept <- -rma_result$Intercept[1]
x_grid <- seq(min(rma_df$x, na.rm = TRUE), max(rma_df$x, na.rm = TRUE), length.out = 200)
set.seed(123)
n_boot <- 1000

boot_predictions <- replicate(n_boot, {
    index <- sample(seq_len(nrow(rma_df)), size = nrow(rma_df), replace = TRUE)
    d <- rma_df[index, ]
    boot_r <- cor(d$x, d$y)
    boot_slope <- sign(boot_r) * sd(d$y)/sd(d$x)
    boot_intercept <- mean(d$y) - boot_slope * mean(d$x)
    -(boot_intercept + boot_slope * x_grid)
})

rma_ci <- data.frame(x = x_grid, fit = rma_intercept + rma_slope * x_grid, lower = apply(boot_predictions, 
    1, quantile, probs = 0.025, na.rm = TRUE), upper = apply(boot_predictions, 1, quantile, probs = 0.975, 
    na.rm = TRUE))

p <- ggplot(imgs_df, aes(LIS_trend, -lambda_ac1_trend)) + geom_hline(yintercept = 0, color = "grey45", 
    size = 0.2, linetype = 2) + geom_vline(xintercept = 0, color = "grey45", size = 0.2, linetype = 2) + 
    geom_pointdensity(size = 0.02) + geom_ribbon(data = rma_ci, aes(x = x, ymin = lower, ymax = upper), 
    inherit.aes = FALSE, fill = "red", alpha = 0.3) + geom_line(data = rma_ci, aes(x = x, y = fit), inherit.aes = FALSE, 
    linewidth = 0.5, color = "red") + scale_color_gradient(low = "gray85", high = "grey20") + guides(color = "none") + 
    labs(x = "Trend of landscape interaction strength", y = "Trend of recovery rate") +
    coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) + theme_bw() + theme(plot.title = element_text(size = 10, hjust = 0.5), legend.key.size = unit(0.3, 
    "cm"), legend.text = element_text(size = 8, hjust = 0), legend.title = element_text(size = 8, face = "bold", 
    hjust = 0), legend.background = element_rect(colour = "NA", fill = "NA", size = 0), strip.text.x = element_text(size = 10, 
    face = "bold", vjust = 0.5), strip.text.y = element_text(size = 10, face = "plain"), axis.title.x = element_text(size = 10, 
    vjust = 0.5, color = "black", margin = margin(r = 0.15)), axis.title.y = element_text(size = 10, 
    hjust = 0.5, vjust = 0.5, color = "black", margin = margin(r = 0.6)), axis.text.x = element_text(size = 6, 
    vjust = 0.5, hjust = 0.5, color = "black", margin = margin(r = 0.3)), axis.text.y = element_text(size = 6, 
    color = "black", hjust = 1, margin = margin(r = 0.3)), axis.ticks = element_line(colour = "black", 
    size = 0.25), axis.ticks.length = unit(0.05, "cm"), plot.margin = unit(c(0.2, 0.02, 0.05, 0.02), 
    "cm"), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), 
    panel.grid.minor.y = element_blank(), panel.spacing.x = unit(0.1, "mm"), panel.spacing.y = unit(0.1, 
        "mm"), legend.position = "top", legend.margin = margin(t = -0, b = -8, unit = "pt"), legend.spacing.y = unit(0.002, 
        "cm"), strip.background = element_blank())
ggsave(plot = p, paste0(outpath, "LIS_and_Recovery_Rate_Trends.png"), width = 5, 
    height = 4, dpi = 900, units = "cm")

{
    library(dplyr)
    library(ggplot2)
    library(parallel)
    library(data.table)
    library(terra)
    clip_roi <- vect(paste0(data_dir, "Image/DrylandRegion/GlobalAridRegion_cutHyperArid.shp"))
    LIS_path <- paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/LIS_10km.tif")
    LIS <- rast(LIS_path) %>% crop(., clip_roi) %>% mask(., clip_roi)
    LIS_years <- rast(paste0(data_dir, "GlobalLISMap1km_Trend/Results/5_Project10km_trend/LIS_tif_10km_Years/LIS_Annual_10km.tif"))
    LIS_years <- project(LIS_years, LIS, method = "near") %>% crop(., LIS, extend = TRUE) %>% mask(., 
        LIS)
    out_path <- paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/")
    Annual_Tem <- rast(paste0(out_path, "Annual_Tem", "_MovingWindows_10km.tif")) %>% 
        crop(., LIS, extend = TRUE) %>% mask(., LIS)
    Annual_Pre <- rast(paste0(out_path, "Annual_Pre", "_MovingWindows_10km.tif")) %>% 
        crop(., LIS, extend = TRUE) %>% mask(., LIS)
    Annual_Wind <- rast(paste0(out_path, "Annual_Wind", "_MovingWindows_10km.tif")) %>% 
        crop(., LIS, extend = TRUE) %>% mask(., LIS)
    result_path <- paste0(data_dir, "CriticalSlowingDown_10km/Results/7RecoverRateTif/residual/")
    for (indexi in c("lambda_ac1")) {
        result_pathi <- paste0(result_path, indexi, "/")
        file_list <- list.files(result_pathi, pattern = paste0("^.*_raster_.*\\.tif$"), full.names = TRUE)
        id_list <- sub("^.*_raster_(.*)\\.tif$", "\\1", file_list) %>% as.numeric()
        id_list <- id_list[order(id_list)]
        for (id in id_list) {
            index_img <- rast(paste0(result_pathi, indexi, "_raster_", id, ".tif"))
            index_img <- project(index_img, LIS, method = "near", align_only = TRUE)
            LIS_trend_id <- crop(LIS_years, index_img, extend = TRUE)
            Annual_Tem_id <- crop(Annual_Tem, index_img, extend = TRUE)
            Annual_Pre_id <- crop(Annual_Pre, index_img, extend = TRUE)
            Annual_Wind_id <- crop(Annual_Wind, index_img, extend = TRUE)
            id_imgs <- rast(list(index_img, Annual_Tem_id, Annual_Pre_id, Annual_Wind_id, LIS_trend_id))
            id_df <- terra::as.data.frame(id_imgs, na.rm = FALSE, xy = TRUE)
            colnames(id_df) <- c("x", "y", paste0("Index_", 1:19), paste0("Tem_", 1:19), paste0("Pre_", 
                1:19), paste0("Wind_", 1:19), paste0("LIS_", 1:19))
            valid_row <- complete.cases(id_df)
            valid_df <- id_df[valid_row, ]
            if (nrow(valid_df) > 2) {
                valid_results <- lapply(1:nrow(valid_df), function(rown) {
                  ts <- valid_df[rown, -1:-2] %>% unlist() %>% unname()
                  Index_ts <- ts[1:19]
                  Index_ts[Index_ts < (-998)] <- NA
                  Tem_ts <- ts[1:19 + 19 * 1]
                  Pre_ts <- ts[1:19 + 19 * 2]
                  Wind_ts <- ts[1:19 + 19 * 3]
                  LIS_ts <- ts[1:19 + 19 * 4]
                  ts_df <- data.frame(Index = Index_ts, Tem = Tem_ts, Pre = Pre_ts, Wind = Wind_ts, LIS = LIS_ts)
                  m1 <- lm(Index ~ Tem + Pre + Wind, data = ts_df, na.action = na.omit)
                  m2 <- lm(Index ~ Tem + Pre + Wind + LIS, data = ts_df, na.action = na.omit)
                  summ1 <- summary(m1)
                  summ2 <- summary(m2)
                  N <- complete.cases(ts_df) %>% sum()
                  rown_results <- cbind(valid_df[rown, 1:2], data.frame(r1 = summ1$r.squared, adjr1 = summ1$adj.r.squared, 
                    AIC1 = AIC(m1), r2 = summ2$r.squared, adjr2 = summ2$adj.r.squared, AIC2 = AIC(m2), 
                    N = N))
                  return(rown_results)
                }) %>% do.call(rbind, .)
                dir.create(paste0(data_dir, "Summary_AridCutHyper/Results/4_RecoverRate/PixelRegression_lambda/csv/", 
                  indexi, "/"), recursive = TRUE, showWarnings = FALSE)
                fwrite(valid_results, paste0(data_dir, "Summary_AridCutHyper/Results/4_RecoverRate/PixelRegression_lambda/csv/", 
                  indexi, "/tile_", id, ".csv"))
                df_result <- id_df[, 1:2]
                df_result[, c("r1", "adjr1", "AIC1", "r2", "adjr2", "AIC2", "N")] <- NA
                df_result[valid_row, ] <- valid_results
                for (coli in c("r1", "adjr1", "AIC1", "r2", "adjr2", "AIC2", "N")) {
                  XYZ <- data.frame(lon = df_result[, 1], lat = df_result[, 2], value = df_result[, coli])
                  img_raster <- rast(XYZ, type = "xyz", crs = crs(index_img), digits = 3)
                  dir.create(paste0(data_dir, "Summary_AridCutHyper/Results/4_RecoverRate/PixelRegression_lambda/tif/", 
                    indexi, "/", coli, "/"), recursive = TRUE, showWarnings = FALSE)
                  writeRaster(img_raster, paste0(data_dir, "Summary_AridCutHyper/Results/4_RecoverRate/PixelRegression_lambda/tif/", 
                    indexi, "/", coli, "/tile_", id, ".tif"), overwrite = TRUE)
                }
            }
        }
    }
    {
        library(terra)
        base_dir <- paste0(data_dir, "Summary_AridCutHyper/Results/4_RecoverRate/PixelRegression_lambda/tif/")
        output_dir <- paste0(data_dir, "Summary_AridCutHyper/Results/4_RecoverRate/PixelRegression_lambda/merged_tif/")
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        index_list <- c("lambda_ac1")
        stat_cols <- c("r1", "adjr1", "AIC1", "r2", "adjr2", "AIC2", "N")
        for (indexi in index_list) {
            for (coli in stat_cols) {
                tif_path <- file.path(base_dir, indexi, coli)
                tif_files <- list.files(tif_path, pattern = "\\.tif$", full.names = TRUE)
                if (length(tif_files) > 0) {
                  rast_list <- lapply(tif_files, rast)
                  merged_rast <- do.call(merge, rast_list)
                  output_file <- file.path(output_dir, paste0(indexi, "_", coli, "_merged.tif"))
                  writeRaster(merged_rast, output_file, overwrite = TRUE)
                }
                else {
                  warning("No tif files found for ", indexi, " - ", coli)
                }
            }
        }
    }
}

library(dplyr)
library(ggplot2)
library(parallel)
library(data.table)
library(terra)
library(sf)
library(ggspatial)
library(scales)
library(tidyr)

{
    clip_roi <- vect(paste0(data_dir, "Image/DrylandRegion/GlobalAridRegion_cutHyperArid.shp"))
    img_path <- paste0(data_dir, "Summary_AridCutHyper/Results/4_RecoverRate/PixelRegression_lambda/merged_tif/")
    LIS_path <- paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/LIS_10km.tif")
    LIS <- rast(LIS_path) %>% crop(., clip_roi) %>% mask(., clip_roi)
    indexi <- "lambda_ac1"
    i <- "adjr"
    m1 <- rast(paste0(img_path, indexi, "_", i, "1_merged.tif")) %>% crop(., LIS, extend = TRUE) %>% 
        mask(., LIS)
    m2 <- rast(paste0(img_path, indexi, "_", i, "2_merged.tif")) %>% crop(., LIS, extend = TRUE) %>% 
        mask(., LIS)
    N <- rast(paste0(img_path, indexi, "_N_merged.tif")) %>% crop(., LIS, extend = TRUE) %>% mask(., 
        LIS)
    dif <- m2 - m1
    m1[N < 19] <- NA
    m2[N < 19] <- NA
    dif[N < 19] <- NA
    plot_df <- rbind(data.frame(value = values(m1, na.rm = TRUE), type = "Without landscape interaction"), 
        data.frame(value = values(m2, na.rm = TRUE) %>% unname(), type = "With landscape interaction"))
    mycolors <- c(`With landscape interaction` = "#FF4D4D", `Without landscape interaction` = "grey75")
    p <- ggplot(plot_df, aes(x = value, fill = type, group = type)) + geom_histogram(data = plot_df[which(plot_df$type == 
        "Without landscape interaction"), ], aes(y = after_stat(count/1000)), alpha = 0.45, color = "white", 
        bins = 48, linewidth = 0.1) + geom_histogram(data = plot_df[which(plot_df$type == "With landscape interaction"), 
        ], aes(y = after_stat(count/1000)), alpha = 0.45, color = "white", bins = 48, linewidth = 0.1) + 
        scale_fill_manual(values = mycolors, name = "") + labs(x = expression(Adjusted~R^2), y = expression(Frequency~(10^3))) + scale_x_continuous(expand = c(0.02, 
        0.02)) + scale_y_continuous(limits = c(NA, 36)) + guides(fill = guide_legend(title = NULL)) + theme_bw() + theme(plot.title = element_text(size = 10, 
        hjust = 0.5), legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 8, hjust = 0), 
        legend.title = element_text(size = 8, face = "bold", hjust = 0), legend.background = element_rect(colour = "NA", 
            fill = "NA", size = 0), strip.text.x = element_text(size = 10, face = "bold", vjust = 0.5), 
        strip.text.y = element_text(size = 10, face = "plain"), axis.title.x = element_text(size = 10, 
            vjust = 0.5, color = "black", margin = margin(r = 0.5)), axis.title.y = element_text(size = 10, 
            hjust = 0.5, vjust = 0.5, color = "black", margin = margin(r = 3)), axis.text.x = element_text(size = 6, 
            vjust = 0.5, hjust = 0.5, color = "black", margin = margin(r = 0.3)), axis.text.y = element_text(size = 6, 
            color = "black", hjust = 1, margin = margin(r = 0.3)), axis.ticks = element_line(colour = "black", 
            size = 0.25), axis.ticks.length = unit(0.05, "cm"), plot.margin = unit(c(0.2, 0.02, 0.05, 
            0.02), "cm"), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.spacing.x = unit(0.1, 
            "mm"), panel.spacing.y = unit(0.1, "mm"), legend.position = "top", legend.margin = margin(t = -0, 
            b = -8, unit = "pt"), legend.spacing.y = unit(0.002, "cm"), strip.background = element_blank())
    ggsave(plot = p, paste0(data_dir, "Summary_AridCutHyper/Results/4_RecoverRate/PixelRegression_lambda/", 
        "Temporal_Recovery_Model_Comparison.png"), width = 4.2, height = 4.4, dpi = 900, units = "cm")
}
