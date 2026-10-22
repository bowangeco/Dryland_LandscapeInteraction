data_dir <- "data/"
output_dir <- "output/"

library(dplyr)
library(ggplot2)
library(parallel)
library(data.table)
library(terra)
library(ggpointdensity)

clip_roi <- vect(paste0(data_dir, "Image/DrylandRegion/GlobalAridRegion_cutHyperArid.shp"))
LIS_path <- paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/LIS_10km.tif")
LIS <- rast(LIS_path) %>% crop(., clip_roi) %>% mask(., clip_roi)
Pre <- rast(paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/Pre_10km.tif")) %>% 
    crop(., LIS, extend = TRUE) %>% mask(., LIS)
Tem <- rast(paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/Tem_10km.tif")) %>% 
    crop(., LIS, extend = TRUE) %>% mask(., LIS)
NDVI <- rast(paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/NDVI_10km.tif")) %>% 
    crop(., LIS, extend = TRUE) %>% mask(., LIS)
Wind <- rast(paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/Wind_10km.tif")) %>% 
    crop(., LIS, extend = TRUE) %>% mask(., LIS)
Human <- rast(paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/Human_10km.tif")) %>% 
    crop(., LIS, extend = TRUE) %>% mask(., LIS)
DEM <- rast(paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/DEM_10km.tif")) %>% 
    crop(., LIS, extend = TRUE) %>% mask(., LIS)
DEM_sd <- rast(paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/DEM_sd_10km.tif")) %>% 
    crop(., LIS, extend = TRUE) %>% mask(., LIS)
sand <- rast(paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/sand_10km.tif")) %>% 
    crop(., LIS, extend = TRUE) %>% mask(., LIS)
sand_sd <- rast(paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/sand_sd_10km.tif")) %>% 
    crop(., LIS, extend = TRUE) %>% mask(., LIS)

water <- rast(paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/water_10km.tif")) %>% 
    crop(., LIS, extend = TRUE) %>% mask(., LIS)
SOC <- rast(paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/SOC_10km.tif")) %>% 
    crop(., LIS, extend = TRUE) %>% mask(., LIS)
GroundWater_rast <- rast(paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/GroundWater.tif")) %>% 
    crop(., LIS, extend = TRUE) %>% mask(., LIS)
MSI_Pre <- rast(paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/MSI_Pre_10km.tif")) %>%
    crop(., LIS, extend = TRUE) %>% mask(., LIS)
MSI_Tem <- rast(paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/MSI_Tem_10km.tif")) %>%
    crop(., LIS, extend = TRUE) %>% mask(., LIS)
MSI_Wind <- rast(paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/MSI_Wind_10km.tif")) %>%
    crop(., LIS, extend = TRUE) %>% mask(., LIS)
Biome2017 <- rast(paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/Biome2017.tif")) %>% 
    crop(., LIS, extend = TRUE) %>% mask(., LIS)
Biome2017[Biome2017 == "N/A"] <- NA
Biome2017 <- mask(Biome2017, Biome2017 %in% c("Deserts & Xeric Shrublands", "Tropical & Subtropical Grasslands, Savannas & Shrublands", 
    "Temperate Grasslands, Savannas & Shrublands", "Mediterranean Forests, Woodlands & Scrub", "Montane Grasslands & Shrublands"), 
    maskvalues = c(NA, 0))
ac1_img <- rast(paste0(data_dir, "CriticalSlowingDown_10km/Results/8_CSDofTotalTime/lambda_tif/lambda_ac1.tif")) %>% 
    crop(., LIS, extend = TRUE) %>% terra::project(., LIS, method = "near") %>% mask(., LIS)
imgs <- rast(list(LIS, Pre, Tem, NDVI, Wind, Human, DEM, DEM_sd, sand, sand_sd, water, SOC, GroundWater_rast, 
    MSI_Pre, MSI_Tem, MSI_Wind, ac1_img, Biome2017))
imgs_df <- as.data.frame(imgs)
colnames(imgs_df) <- c("LIS", "Pre", "Tem", "NDVI", "Wind", "Human", "DEM", "DEM_sd", "sand", "sand_sd", 
    "water", "SOC", "GroundWater", "MSI_Pre", "MSI_Tem", "MSI_Wind", "lambda_ac1", "Ecoregion")
predictor_columns <- setdiff(names(imgs_df), c("lambda_ac1", "Ecoregion"))
imgs_df[predictor_columns] <- lapply(imgs_df[predictor_columns], function(x) as.numeric(scale(x)))

imgs_df_global <- imgs_df
imgs_df_global$Ecoregion <- "Global"
imgs_df <- rbind(imgs_df, imgs_df_global)
imgs_df <- imgs_df[complete.cases(imgs_df), ]
Econames <- unique(imgs_df$Ecoregion)
imgs_df$Ecoregion <- factor(imgs_df$Ecoregion, levels = c("Global", "Deserts & Xeric Shrublands", "Tropical & Subtropical Grasslands, Savannas & Shrublands", 
    "Temperate Grasslands, Savannas & Shrublands", "Mediterranean Forests, Woodlands & Scrub", "Montane Grasslands & Shrublands"), 
    labels = c("Global", "Desert", "Tropical", "Temperate", "Mediterranean", "Montane"))

library(dplyr)
library(tidyr)
library(ggplot2)

imgs_df_lm <- imgs_df[imgs_df$Ecoregion == "Global", ]

{
    library(MuMIn)
    options(na.action = "na.fail")
    full_model <- lm(lambda_ac1 ~ Pre + Tem + Wind + NDVI + DEM + DEM_sd + sand + sand_sd + water + SOC + 
        GroundWater + Human + MSI_Pre + MSI_Tem + MSI_Wind + LIS, data = imgs_df_lm)
    # Select the environmental baseline without LIS for model comparison.
    environment_model <- update(full_model, . ~ . - LIS)
    model_set <- dredge(environment_model, rank = "AIC")
    outpath <- paste0(data_dir, "Summary_AridCutHyper/Results/4_RecoverRate/OnTotalTime/")
    dir.create(outpath, recursive = TRUE, showWarnings = FALSE)
    fwrite(model_set, paste0(outpath, "Environmental_Model_Selection.csv"))
    best_model <- model_set[1, ]
    selected_baseline <- get.models(model_set, 1)[[1]]
    included_vars <- data.frame(term = names(coef(selected_baseline)), estimate = unname(coef(selected_baseline)))
    included_vars$AIC <- best_model$AIC
    included_vars$adjR2 <- summary(get.models(model_set, 1)[[1]])$adj.r.squared
    # Fit lambda first; convert standardized effects to recovery rate below.
    coefficient_model <- full_model
    extract_ci <- function(model, model_name) {
        frame <- model.frame(model)
        predictor_sd <- apply(model.matrix(model), 2, sd)
        response_sd <- sd(model.response(frame))
        scale_factor <- predictor_sd/response_sd
        scale_factor["(Intercept)"] <- NA_real_
        estimates <- coef(model) * scale_factor
        intervals <- confint(model) * scale_factor
        data.frame(term = names(estimates), Standardized_Beta = estimates, CI_low = intervals[, 1], CI_high = intervals[, 
            2], model = model_name, stringsAsFactors = FALSE)
    }
    coef_full <- extract_ci(coefficient_model, "full")
    fwrite(coef_full, paste0(outpath, "Lambda_Standardized_Coefficients.csv"))
    coef_full$Standardized_Beta <- -coef_full$Standardized_Beta
    coef_full[c("CI_low", "CI_high")] <- -coef_full[c("CI_high", "CI_low")]
    fwrite(coef_full, paste0(outpath, "Recovery_Standardized_Coefficients.csv"))
    coef_df <- bind_rows(coef_full) %>% filter(term != "(Intercept)")
    name_map <- c(Human = "Human activity", LIS = "Landscape interaction", Pre = "Precipitation", Tem = "Temperature", 
        NDVI = "NDVI", Wind = "Wind speed", DEM = "Elevation", DEM_sd = "Elevation SD", sand = "Sand fraction", 
        sand_sd = "Sand SD", water = "Soil Moisture", SOC = "Soil organic carbon", GroundWater = "Groundwater",
        MSI_Pre = "MSI (precipitation)", MSI_Tem = "MSI (temperature)", MSI_Wind = "MSI (wind speed)")
    coef_df$Variable_label <- name_map[coef_df$term]
    sort_order <- coef_df %>% filter(model == "full") %>% arrange(desc(Standardized_Beta)) %>% pull(Variable_label)
    coef_df$Variable_label <- factor(coef_df$Variable_label, levels = sort_order %>% rev())
    p <- ggplot(coef_df[which(coef_df$model == "full"), ], aes(x = Variable_label, y = Standardized_Beta, 
        fill = model)) + geom_col(position = position_dodge(width = 0.7), width = 0.7, alpha = 0.9) + 
        geom_errorbar(aes(ymin = CI_low, ymax = CI_high), position = position_dodge(width = 0.7), width = 0.25, 
            linewidth = 0.3, color = "black") + geom_hline(yintercept = 0, color = "gray40", linetype = "dashed", 
        linewidth = 0.3) + scale_fill_manual(values = c(full = "#7091CB")) + scale_x_discrete(expand = expansion(mult = c(0.04, 
        0.04))) + labs(x = NULL, y = "Standardized coefficient", fill = NULL) + guides(fill = "none") + 
        coord_flip() + theme_bw(base_size = 14) + theme(axis.title.x = element_text(size = 8, vjust = 2.5, 
        color = "black"), axis.title.y = element_text(size = 8, hjust = 0.5, vjust = 0.5, color = "black", 
        margin = margin(r = 0.1, unit = "mm")), axis.text.x = element_text(size = 8, angle = 0, hjust = 1, 
        color = "black", margin = margin(t = -0.1, unit = "mm")), axis.text.y = element_text(size = 8, 
        color = "black"), axis.ticks = element_line(colour = "black", size = 0.25), axis.ticks.length = unit(0.05, 
        "cm"), plot.margin = unit(c(0.02, 0.02, 0.02, 0.4), "cm"), legend.position = c(0.55, 0.85))
    outpath <- paste0(data_dir, "Summary_AridCutHyper/Results/4_RecoverRate/OnTotalTime/")
    ggsave(plot = p, paste0(outpath, "Recovery_Standardized_Coefficients.jpeg"), 
        width = 6.5, height = 5, dpi = 900, units = "cm")
}

results <- data.frame(Ecoregion = character(), R2_m0 = numeric(), AIC_m0 = numeric(), R2_m1 = numeric(), 
    AIC_m1 = numeric(), stringsAsFactors = FALSE)

for (eco in levels(imgs_df$Ecoregion)) {
    sub_df <- subset(imgs_df, Ecoregion == eco)
    environment_model <- lm(lambda_ac1 ~ Pre + Tem + Wind + NDVI + DEM + DEM_sd + sand + sand_sd + 
        water + SOC + GroundWater + Human + MSI_Pre + MSI_Tem + MSI_Wind, data = sub_df, na.action = na.fail)
    environment_set <- MuMIn::dredge(environment_model, rank = "AIC")
    m0 <- MuMIn::get.models(environment_set, 1)[[1]]
    m1 <- update(m0, . ~ . + LIS)
    results <- rbind(results, data.frame(Ecoregion = eco, R2_m0 = round(summary(m0)$adj.r.squared, 3), 
        AIC_m0 = round(AIC(m0), 1), R2_m1 = round(summary(m1)$adj.r.squared, 3), AIC_m1 = round(AIC(m1), 
            1)))
}

results$R2_dif <- ((results$R2_m1 - results$R2_m0)/results$R2_m0 * 100) %>% round()
results$AIC_dif <- ((results$AIC_m1 - results$AIC_m0)/results$AIC_m0 * 100) %>% round(.)
fwrite(results, paste0(data_dir, "Summary_AridCutHyper/Results/4_RecoverRate/OnTotalTime/Recovery_Model_Comparison.csv"))

library(ggplot2)
library(dplyr)
library(tidyr)

results <- fread(paste0(data_dir, "Summary_AridCutHyper/Results/4_RecoverRate/OnTotalTime/Recovery_Model_Comparison.csv")) %>% 
    as.data.frame()
plot_df <- results %>% dplyr::select(Ecoregion, R2_m0, R2_m1, AIC_m0, AIC_m1) %>% tidyr::pivot_longer(cols = -Ecoregion, 
    names_to = c("Metric", "Model"), names_pattern = "(.*)_m(\\d)") %>% dplyr::mutate(Model = factor(Model, 
    levels = c("0", "1"), labels = c("m0", "m1")), fill_color = ifelse(Model == "m1", "red", "grey65"))
plot_df$Ecoregion <- factor(plot_df$Ecoregion, levels = c("Global", "Desert", "Tropical", "Temperate", 
    "Mediterranean", "Montane"), labels = c("Global", "Desert", "Tropical", "Temperate", "Mediterranean", 
    "Montane"))

{
    aic_min <- min(plot_df$value[plot_df$Metric == "AIC"])
    aic_scaler <- abs(aic_min) * 2
    plot_df2 <- plot_df %>% mutate(value_scaled = case_when(Metric == "AIC" ~ value/aic_scaler, TRUE ~ 
        value), label_value = round(value, 2))
    df_r2 <- plot_df2 %>% filter(Metric == "R2")
    df_aic <- plot_df2 %>% filter(Metric == "AIC")
    df_aic_diff <- df_aic %>% dplyr::select(Ecoregion, Model, label_value) %>% pivot_wider(names_from = Model, 
        values_from = label_value) %>% mutate(diff = m1 - m0)
    df_aic_diff$diff_zoom <- sign(df_aic_diff$diff) * log10(abs(df_aic_diff$diff))
    df_aic_diff$diff_zoom <- df_aic_diff$diff_zoom/12
    p <- ggplot() + geom_col(data = df_r2 %>% filter(Model == "m1"), aes(x = Ecoregion, y = value_scaled, fill = Model), 
        width = 0.6, alpha = 1) + geom_col(data = df_r2 %>% filter(Model == "m0"), 
        aes(x = Ecoregion, y = value_scaled, fill = Model), width = 0.6) + geom_col(data = df_aic_diff, 
        aes(x = Ecoregion, y = diff_zoom), width = 0.6, fill = "#7091CB", alpha = 1) + geom_hline(yintercept = 0, 
        color = "black", linewidth = 0.5) + scale_y_continuous(breaks = c(-log10(1e+05)/12, -log10(10000)/12, 
        -log10(1000)/12, -log10(100)/12, -log10(10)/12, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), limits = c(-0.43, 
        0.8), labels = function(x) {
        sapply(x, function(v) {
            if (v < 0) {
                -10^(-v * 12)
            }
            else {
                round(v, 2)
            }
        })
    }) + scale_fill_manual(values = c(m1 = "#FF4D4D", m0 = "grey75"),
        breaks = c("m1", "m0"), labels = c("With landscape interaction", "Without landscape interaction"), name = NULL) +
        labs(x = NULL, y = expression(Adjusted~R^2~"/"~Delta*AIC)) + theme_bw(base_size = 14) + theme(axis.text.x = element_text(size = 8, 
        angle = 30, hjust = 1, margin = margin(t = -0.1, unit = "mm"), color = "black"), axis.text.y = element_text(size = 6, 
        color = "black"), axis.title.y = element_text(size = 8, color = "black"), axis.ticks = element_line(size = 0.25), 
        axis.ticks.length = unit(0.05, "cm"), plot.margin = unit(c(0.04, 0.04, 0.03, 0.04), "cm"), panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), legend.position = "top", legend.direction = "vertical")
    outpath <- paste0(data_dir, "Summary_AridCutHyper/Results/4_RecoverRate/OnTotalTime/")
    ggsave(plot = p, paste0(outpath, "Recovery_Model_Comparison.jpeg"), 
        width = 6, height = 5.2, dpi = 900, units = "cm")
}
