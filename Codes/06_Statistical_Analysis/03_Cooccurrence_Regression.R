data_dir <- "data/"
output_dir <- "output/"

library(dplyr)
library(ggplot2)
library(parallel)
library(data.table)
library(terra)
library(pROC)

{
    library(dplyr)
    library(ggplot2)
    library(parallel)
    library(data.table)
    library(terra)
    library(pROC)
    clip_roi <- vect(paste0(data_dir, "Image/DrylandRegion/GlobalAridRegion_cutHyperArid.shp"))
    LIS_path <- paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/LIS_10km.tif")
    LIS <- rast(LIS_path) %>% crop(., clip_roi) %>% mask(., clip_roi)
    Coex_10km <- rast(paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/Cooccurrence_10km.tif")) %>% 
        crop(., LIS, extend = TRUE) %>% mask(., LIS)
    data_path <- paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/")
    Biome2017 <- rast(paste0(data_path, "Biome2017.tif")) %>% crop(., LIS, extend = TRUE) %>% mask(., 
        LIS)
    Biome2017[Biome2017 == "N/A"] <- NA
    Biome2017 <- mask(Biome2017, Biome2017 %in% c("Deserts & Xeric Shrublands", "Tropical & Subtropical Grasslands, Savannas & Shrublands", 
        "Temperate Grasslands, Savannas & Shrublands", "Mediterranean Forests, Woodlands & Scrub", "Montane Grasslands & Shrublands"), 
        maskvalues = c(NA, 0))
    imgs <- rast(list(LIS, Coex_10km, Biome2017))
    imgs_df <- as.data.frame(imgs)
    colnames(imgs_df) <- c("LIS", "Coex", "Ecoregion")
    imgs_df_global <- imgs_df
    imgs_df_global$Ecoregion <- "Global"
    imgs_df <- rbind(imgs_df, imgs_df_global)
    imgs_df <- imgs_df[complete.cases(imgs_df), ]
    Econames <- unique(imgs_df$Ecoregion)
    Total_hist_all <- lapply(1:length(Econames), function(i) {
        namei <- Econames[i]
        sub_df <- imgs_df[imgs_df$Ecoregion == namei, 1:2]
        if (nrow(sub_df) > 1000) {
            hist1 <- hist(sub_df$LIS[which(sub_df$Coex == 1)], breaks = seq(0, 1, 0.04), plot = FALSE)
            hist0 <- hist(sub_df$LIS[which(sub_df$Coex == 0)], breaks = seq(0, 1, 0.04), plot = FALSE)
            Total_hist <- data.frame(Coex = "Coexistence", LIS = hist1$mids, Prob = hist1$counts/(hist1$counts + 
                hist0$counts), count = hist1$counts + hist0$counts) %>% rbind(., data.frame(Coex = "Non-coexistence", 
                LIS = hist0$mids, Prob = hist0$counts/(hist1$counts + hist0$counts), count = hist0$counts + 
                  hist1$counts))
            Total_hist$Prob[which(Total_hist$count < 10)] <- NA
            Total_hist$Ecoregion <- namei
        }
        else {
            return(NULL)
        }
        return(Total_hist)
    }) %>% do.call(rbind, .)
    Total_hist_all$Ecoregion <- factor(Total_hist_all$Ecoregion, levels = c("Global", "Deserts & Xeric Shrublands", 
        "Tropical & Subtropical Grasslands, Savannas & Shrublands", "Temperate Grasslands, Savannas & Shrublands", 
        "Mediterranean Forests, Woodlands & Scrub", "Montane Grasslands & Shrublands"), labels = c("Global", 
        "Desert", "Tropical", "Temperate", "Mediterranean", "Montane"))
    imgs_df$Ecoregion <- factor(imgs_df$Ecoregion, levels = c("Global", "Deserts & Xeric Shrublands", 
        "Tropical & Subtropical Grasslands, Savannas & Shrublands", "Temperate Grasslands, Savannas & Shrublands", 
        "Mediterranean Forests, Woodlands & Scrub", "Montane Grasslands & Shrublands"), labels = c("Global", 
        "Desert", "Tropical", "Temperate", "Mediterranean", "Montane"))
    library(scales)
    Total_hist_all$Prob[Total_hist_all$count < 30] <- NA
    Total_hist_all$Prob[Total_hist_all$Ecoregion == "Global" & Total_hist_all$count < 100] <- NA
    p <- ggplot(imgs_df, aes(LIS, Coex, group = Ecoregion)) + geom_col(data = Total_hist_all[which(Total_hist_all$Coex == 
        "Coexistence"), ], aes(x = LIS, y = Prob, fill = Coex), position = "dodge") + facet_wrap(~Ecoregion, 
        nrow = 2, scales = "free_y") + scale_fill_manual(breaks = c("Coexistence", "Non-coexistence"), 
        values = c("#7091CB", "grey65")) + stat_smooth(geom = "line", method = "glm", method.args = list(family = binomial), 
        formula = "y~x", fullrange = TRUE, fill = "grey85", color = "red", se = TRUE, size = 1, alpha = 0.9) + 
        labs(x = "Landscape interaction strength", y = "Local co-occurrence probability (%)") + scale_x_continuous(breaks = c(0, 
        0.2, 0.4, 0.6, 0.8), expand = c(0.005, 0.005), limits = c(-0.01 + 0.3, 0.91)) + scale_y_continuous(expand = expansion(mult = c(0.02, 
        0.1), add = c(0, 0)), limits = c(0, NA), labels = function(x) x * 100, breaks = pretty_breaks(n = 4)) + 
        guides(fill = "none") + theme_bw() + theme(plot.title = element_text(size = 8, hjust = 0.5), 
        legend.text = element_text(size = 8, hjust = 0), legend.title = element_text(size = 6, hjust = 0), 
        legend.background = element_rect(colour = "NA", fill = "NA", size = 0), strip.text.x = element_text(size = 7, 
            vjust = 0.5, margin = margin(r = 0.13)), strip.text.y = element_text(size = 8, face = "plain"), 
        axis.title.x = element_text(size = 8, vjust = 0.5, color = "black", margin = margin(r = 0.3)), 
        axis.title.y = element_text(size = 8, hjust = 0.5, vjust = 0.5, color = "black", margin = margin(r = 0.3)), 
        axis.text.x = element_text(size = 6, vjust = 0.5, hjust = 0.5, color = "black", margin = margin(r = 0.2)), 
        axis.text.y = element_text(size = 6, color = "black", hjust = 0.5, margin = margin(r = 0.15)), 
        axis.ticks = element_line(colour = "black", size = 0.25), axis.ticks.length = unit(0.05, "cm"), 
        plot.margin = unit(c(0.02, 0.02, 0.03, 0.02), "cm"), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.spacing.x = unit(1, 
            "mm"), panel.spacing.y = unit(0.8, "mm"), legend.position = c(0.08, 0.78), strip.background = element_blank())
    outpath <- paste0(data_dir, "Summary_AridCutHyper/Results/3_Coexistence/Correlation/Logist/CalculateIn10km/")
    ggsave(plot = p, paste0(outpath, "Cooccurrence_Probability_by_Biome.png"), width = 10.8, height = 5.5, 
        dpi = 900, units = "cm")
}

{
    outpath <- paste0(data_dir, "Summary_AridCutHyper/Results/3_Coexistence/EnhanceModel/Logist/Calculate10km_rep99/")
    Pre <- rast(paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/Pre_10km.tif")) %>% 
        crop(., LIS, extend = TRUE) %>% mask(., LIS)
    Tem <- rast(paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/Tem_10km.tif")) %>% 
        crop(., LIS, extend = TRUE) %>% mask(., LIS)
    Wind <- rast(paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/Wind_10km.tif")) %>% 
        crop(., LIS, extend = TRUE) %>% mask(., LIS)
    imgs <- rast(list(LIS, Coex_10km, Pre, Wind, Tem))
    imgs_df <- as.data.frame(imgs)
    colnames(imgs_df) <- c("LIS", "Coex", "Pre", "Wind", "Tem")
    imgs_df <- imgs_df[complete.cases(imgs_df), ]
    n <- nrow(imgs_df)
    AUC0_list <- c()
    AUC1_list <- c()
    ROC_df <- data.frame(Specificity = numeric(), Sensitivity = numeric(), Model = numeric())
    set.seed(428)
    for (j in seq_len(1000)) {
        samplie_index <- sample(1:n, 0.8 * n)
        train_df <- imgs_df[samplie_index, c("Coex", "LIS", "Pre", "Wind", "Tem")]
        test_df <- imgs_df[-samplie_index, c("Coex", "LIS", "Pre", "Wind", "Tem")]
        if (sum(test_df$Coex) == 0) {
            AUC0 <- NA
            AUC1 <- NA
            df0 <- data.frame(Specificity = NA, Sensitivity = NA, Model = 0)
            df1 <- data.frame(Specificity = NA, Sensitivity = NA, Model = 1)
        }
        else {
            mod0 <- glm(Coex ~ Pre + Wind + Tem + I(Pre^2) + I(Wind^2) + I(Tem^2), train_df, family = binomial(link = "logit"))
            fitted0 <- predict(mod0, newdata = test_df, type = "response")
            ROC0 <- roc(response = test_df$Coex, predictor = fitted0)
            AUC0 <- auc(response = test_df$Coex, predictor = fitted0)
            df0 <- data.frame(Specificity = ROC0$specificities, Sensitivity = ROC0$sensitivities, Model = 0)
            mod1 <- glm(Coex ~ Pre + Wind + Tem + I(Pre^2) + I(Wind^2) + I(Tem^2) + LIS, train_df, family = binomial(link = "logit"))
            fitted1 <- predict(mod1, newdata = test_df, type = "response")
            ROC1 <- roc(response = test_df$Coex, predictor = fitted1)
            AUC1 <- auc(response = test_df$Coex, predictor = fitted1)
            df1 <- data.frame(Specificity = ROC1$specificities, Sensitivity = ROC1$sensitivities, Model = 1)
        }
        AUC0_list <- c(AUC0_list, AUC0)
        AUC1_list <- c(AUC1_list, AUC1)
        df <- rbind(df0, df1)
        ROC_df <- rbind(ROC_df, df)
    }
    fwrite(ROC_df, paste0(outpath, "ROC_Global.csv"))
    meanAUC0 <- mean(AUC0_list, na.rm = TRUE)
    meanAUC1 <- mean(AUC1_list, na.rm = TRUE)
    imgs_df_scaled <- scale(imgs_df[, c("LIS", "Pre", "Wind", "Tem")]) %>% cbind(data.frame(Coex = imgs_df$Coex), 
        .) %>% as.data.frame()
    m0 <- glm(data = imgs_df_scaled, formula = Coex ~ Pre + Wind + Tem + I(Pre^2) + I(Wind^2) + I(Tem^2), 
        family = binomial(link = "logit"))
    m1 <- glm(data = imgs_df_scaled, formula = Coex ~ Pre + Wind + Tem + I(Pre^2) + I(Wind^2) + I(Tem^2) + 
        LIS, family = binomial(link = "logit"))
    aov <- anova(m0, m1, test = "LR")
    aov_p <- aov$`Pr(>Chi)`
    m1_coef <- m1$coefficients
    m0_coef <- m0$coefficients
    m1_coef <- c(m1_coef, c(AIC = AIC(m1), AUC = meanAUC1, Model = "m1", p = aov_p))
    m0_coef <- c(m0_coef, c(LIS = NA, AIC = AIC(m0), AUC = meanAUC0, Model = "m0", p = aov_p))
    df <- rbind(m1_coef, m0_coef)
    Out_path <- paste0(outpath, "Cooccurrence_Regression_Global.csv")
    fwrite(df, Out_path)
    {
        library(parallel)
        result_path <- paste0(outpath, "Plot_ROC/")
        if (!dir.exists(result_path)) {
            dir.create(result_path)
        }
        for (Model in 0:1) {
            Dune_df <- ROC_df[which(ROC_df$Model == Model), ]
            mean_fun <- function(x) {
                x_min <- x
                x_max <- x + by
                y <- Dune_df$Sensitivity[which(Dune_df$Specificity >= x_min & Dune_df$Specificity < x_max)]
                mean <- mean(y, na.rm = TRUE)
                sd <- sd(y, na.rm = TRUE)
                df <- data.frame(Model = Model, Specificity = x, mean = mean, sd = sd)
                return(df)
            }
            multifun <- function(id) {
                library(dplyr)
                library(data.table)
                sub_seq <- df_list[[id]]
                result_df <- lapply(X = sub_seq, FUN = mean_fun) %>% do.call(rbind, .)
                fwrite(result_df, file = paste0(result_path, "ROC_mean_", id, ".csv"), row.names = FALSE)
            }
            cores <- 2
            cl <- makeCluster(cores)
            clusterExport(cl, c("data_dir", "output_dir"))
            by <- 0.001
            x_seq <- seq(0, 1, by = by)
            n <- cores
            nr <- length(x_seq)
            if (ceiling(nr/n) < n) {
                break
            }
            df_list <- split(x_seq, rep(1:n, each = ceiling(nr/n), length.out = nr))
            clusterExport(cl, varlist = c("Dune_df", "mean_fun", "multifun", "Model", "result_path", 
                "by", "df_list"))
            parLapply(cl = cl, X = 1:n, fun = multifun)
            stopCluster(cl)
            file_list <- list.files(result_path, pattern = paste0("^ROC_mean_\\d+"), full.names = TRUE)
            id <- sub(".*ROC_mean_(\\d*)\\.csv$", "\\1", file_list) %>% as.numeric()
            file_list <- file_list[order(id)]
            Total_result <- lapply(file_list, FUN = function(p) {
                fread(p) %>% data.frame()
            }) %>% do.call(rbind, .)
            fwrite(Total_result, paste0(result_path, "All_ROC_", Model, ".csv"), row.names = FALSE)
            file.remove(file_list)
        }
        file_list <- list.files(result_path, pattern = paste0("^All_ROC_.+"), full.names = TRUE)
        Total_result <- lapply(file_list, FUN = function(p) {
            fread(p) %>% data.frame()
        }) %>% do.call(rbind, .)
        fwrite(Total_result, paste0(result_path, "ROC_Summary.csv"), row.names = FALSE)
        Total_result <- fread(paste0(result_path, "ROC_Summary.csv"))
        Total_result <- rbind(Total_result, data.frame(Model = 0, Specificity = 1, mean = 0, sd = NA), 
            data.frame(Model = 1, Specificity = 1, mean = 0, sd = NA))
        Total_result$Model <- factor(Total_result$Model, levels = 0:1, labels = c("Without landscape interaction", "With landscape interaction"))
        colors <- c("grey65", "red")
        P_ROC <- ggplot(Total_result, aes(x = 1 - Specificity, y = mean, group = Model)) + geom_area(data = Total_result[which(Total_result$Model == 
            "With landscape interaction"), ], fill = colors[2], position = "identity", alpha = 0.7) + geom_area(data = Total_result[which(Total_result$Model == 
            "Without landscape interaction"), ], fill = "grey90", position = "identity", alpha = 1) + geom_abline(slope = 1, 
            intercept = 0, color = "grey45", size = 0.2, linetype = 2) + geom_line(aes(color = Model, 
            group = Model), linewidth = 0.2) + coord_equal() + scale_x_continuous(breaks = seq(0, 1, 
            0.2), expand = c(0.01, 0.01), limits = c(0, 1)) + scale_y_continuous(breaks = seq(0, 1, 0.2), 
            expand = c(0.01, 0.01), limits = c(0, 1)) + guides(fill = "none", color = "none") + scale_fill_manual(values = colors) + 
            scale_color_manual(values = colors) + labs(x = "FPR", y = "TPR") + theme_bw() + theme(plot.title = element_text(size = 10, 
            hjust = 0.5), legend.text = element_text(size = 10, hjust = 0), legend.title = element_text(size = 10, 
            face = "bold", hjust = 0), legend.background = element_rect(colour = "NA", fill = "NA", size = 0), 
            legend.direction = "horizontal", strip.background = element_rect(colour = "NA", fill = "NA", 
                size = 0), strip.text.x = element_text(size = 10, face = "bold"), strip.text.y = element_text(size = 10, 
                face = "bold"), axis.title.x = element_text(size = 8, vjust = 0.5, hjust = 0.5, color = "black", 
                margin = margin(r = 0.15)), axis.title.y = element_text(size = 8, hjust = 0.5, vjust = 0.5, 
                color = "black", margin = margin(r = 0.2)), axis.text.x = element_blank(), axis.text.y = element_blank(), 
            axis.ticks = element_line(colour = "black", size = 0.15), axis.ticks.length = unit(0.01, 
                "cm"), plot.margin = unit(c(0.03, 0.03, 0.03, 0.03), "cm"), panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), 
            panel.spacing.x = unit(2, "mm"), panel.spacing.y = unit(1, "mm"), panel.background = element_rect(colour = "NA", 
                fill = "NA", size = 0), legend.position = c(0.8, -0.11), plot.background = element_rect(fill = "transparent", 
                colour = NA, size = 0.02), axis.line = element_line(colour = "black", linewidth = 0.2), 
            panel.grid = element_blank(), )
        ggsave(plot = P_ROC, paste0(outpath, "ROC_Global.png"), width = 2.2, height = 2.1, dpi = 900, units = "cm")
    }
    {
        for (regions in "Biome2017") {
            outpath <- paste0(data_dir, "Summary_AridCutHyper/Results/3_Coexistence/EnhanceModel/Logist/", 
                regions, "/")
            if (!dir.exists(outpath)) {
                dir.create(outpath)
            }
            region_img <- get(regions)
            imgs <- rast(list(LIS, Coex_10km, Pre, Wind, Tem, region_img))
            imgs_df <- as.data.frame(imgs)
            colnames(imgs_df) <- c("LIS", "Coex", "Pre", "Wind", "Tem", "Ecoregion")
            imgs_df <- imgs_df[complete.cases(imgs_df), ]
            Econames <- unique(imgs_df$Ecoregion)
            for (i in 1:length(Econames)) {
                namei <- Econames[i]
                sub_df <- imgs_df[imgs_df$Ecoregion == namei, 1:5]
                if (nrow(sub_df) > 1000) {
                  n <- nrow(sub_df)
                  AUC0_list <- c()
                  AUC1_list <- c()
                  ROC_df <- data.frame(Specificity = numeric(), Sensitivity = numeric(), Model = numeric())
                  set.seed(428)
                  for (j in seq_len(1000)) {
                    samplie_index <- sample(1:n, 0.8 * n)
                    train_df <- sub_df[samplie_index, c("Coex", "LIS", "Pre", "Wind", "Tem")]
                    test_df <- sub_df[-samplie_index, c("Coex", "LIS", "Pre", "Wind", "Tem")]
                    if (sum(test_df$Coex) == 0) {
                      AUC0 <- NA
                      AUC1 <- NA
                      df0 <- data.frame(Specificity = NA, Sensitivity = NA, Model = 0)
                      df1 <- data.frame(Specificity = NA, Sensitivity = NA, Model = 1)
                    }
                    else {
                      mod0 <- glm(Coex ~ Pre + Wind + Tem + I(Pre^2) + I(Wind^2) + I(Tem^2), train_df, 
                        family = binomial(link = "logit"))
                      fitted0 <- predict(mod0, newdata = test_df, type = "response")
                      ROC0 <- roc(response = test_df$Coex, predictor = fitted0)
                      AUC0 <- auc(response = test_df$Coex, predictor = fitted0)
                      df0 <- data.frame(Specificity = ROC0$specificities, Sensitivity = ROC0$sensitivities, 
                        Model = 0)
                      mod1 <- glm(Coex ~ Pre + Wind + Tem + I(Pre^2) + I(Wind^2) + I(Tem^2) + LIS, train_df, 
                        family = binomial(link = "logit"))
                      fitted1 <- predict(mod1, newdata = test_df, type = "response")
                      ROC1 <- roc(response = test_df$Coex, predictor = fitted1)
                      AUC1 <- auc(response = test_df$Coex, predictor = fitted1)
                      df1 <- data.frame(Specificity = ROC1$specificities, Sensitivity = ROC1$sensitivities, 
                        Model = 1)
                    }
                    AUC0_list <- c(AUC0_list, AUC0)
                    AUC1_list <- c(AUC1_list, AUC1)
                    df <- rbind(df0, df1)
                    ROC_df <- rbind(ROC_df, df)
                  }
                  fwrite(ROC_df, paste0(outpath, namei, "_ROC_df.csv"))
                  meanAUC0 <- mean(AUC0_list, na.rm = TRUE)
                  meanAUC1 <- mean(AUC1_list, na.rm = TRUE)
                  sub_df_scaled <- scale(sub_df[, c("LIS", "Pre", "Wind", "Tem")]) %>% cbind(data.frame(Coex = sub_df$Coex), 
                    .) %>% as.data.frame()
                  m0 <- glm(data = sub_df_scaled, formula = Coex ~ Pre + Wind + Tem + I(Pre^2) + I(Wind^2) + 
                    I(Tem^2), family = binomial(link = "logit"))
                  m1 <- glm(data = sub_df_scaled, formula = Coex ~ Pre + Wind + Tem + I(Pre^2) + I(Wind^2) + 
                    I(Tem^2) + LIS, family = binomial(link = "logit"))
                  aov <- anova(m0, m1, test = "LR")
                  aov_p <- aov$`Pr(>Chi)`
                  m1_coef <- m1$coefficients
                  m0_coef <- m0$coefficients
                  m1_coef <- c(m1_coef, c(AIC = AIC(m1), AUC = meanAUC1, Model = "m1", name = namei, 
                    p = aov_p))
                  m0_coef <- c(m0_coef, c(LIS = NA, AIC = AIC(m0), AUC = meanAUC0, Model = "m0", name = namei, 
                    p = aov_p))
                  df <- rbind(m1_coef, m0_coef)
                  Out_path <- paste0(outpath, namei, "_Regression_df.csv")
                  fwrite(df, file = Out_path)
                  {
                    library(parallel)
                    result_path <- paste0(outpath)
                    for (Model in 0:1) {
                      Dune_df <- ROC_df[which(ROC_df$Model == Model), ]
                      mean_fun <- function(x) {
                        x_min <- x
                        x_max <- x + by
                        y <- Dune_df$Sensitivity[which(Dune_df$Specificity >= x_min & Dune_df$Specificity < 
                          x_max)]
                        mean <- mean(y, na.rm = TRUE)
                        sd <- sd(y, na.rm = TRUE)
                        df <- data.frame(Model = Model, Specificity = x, mean = mean, sd = sd)
                        return(df)
                      }
                      multifun <- function(id) {
                        library(dplyr)
                        library(data.table)
                        sub_seq <- df_list[[id]]
                        result_df <- lapply(X = sub_seq, FUN = mean_fun) %>% do.call(rbind, .)
                        fwrite(result_df, file = paste0(result_path, "ROC_mean_", id, ".csv"), row.names = FALSE)
                      }
                      cores <- 2
                      cl <- makeCluster(cores)
                      clusterExport(cl, c("data_dir", "output_dir"))
                      by <- 0.001
                      x_seq <- seq(0, 1, by = by)
                      n <- cores
                      nr <- length(x_seq)
                      if (ceiling(nr/n) < n) {
                        break
                      }
                      df_list <- split(x_seq, rep(1:n, each = ceiling(nr/n), length.out = nr))
                      clusterExport(cl, varlist = c("Dune_df", "mean_fun", "multifun", "Model", "result_path", 
                        "by", "df_list"))
                      parLapply(cl = cl, X = 1:n, fun = multifun)
                      stopCluster(cl)
                      file_list <- list.files(result_path, pattern = paste0("^ROC_mean_\\d+"), full.names = TRUE)
                      id <- sub(".*ROC_mean_(\\d*)\\.csv$", "\\1", file_list) %>% as.numeric()
                      file_list <- file_list[order(id)]
                      Total_result <- lapply(file_list, FUN = function(p) {
                        fread(p) %>% data.frame()
                      }) %>% do.call(rbind, .)
                      fwrite(Total_result, paste0(result_path, "All_ROC_", Model, ".csv"), row.names = FALSE)
                      file.remove(file_list)
                    }
                    file_list <- list.files(result_path, pattern = paste0("^All_ROC_.+"), full.names = TRUE)
                    Total_result <- lapply(file_list, FUN = function(p) {
                      fread(p) %>% data.frame()
                    }) %>% do.call(rbind, .)
                    fwrite(Total_result, paste0(result_path, namei, "_ROC_Summary.csv"), row.names = FALSE)
                    file.remove(file_list)
                  }
                }
            }
        }
    }
    {
        for (regions in "Biome2017") {
            outpath <- paste0(data_dir, "Summary_AridCutHyper/Results/3_Coexistence/EnhanceModel/Logist/", 
                regions, "/")
            if (!dir.exists(outpath)) {
                dir.create(outpath)
            }
            region_img <- get(regions)
            imgs <- rast(list(LIS, Coex_10km, Pre, Wind, Tem, region_img))
            imgs_df <- as.data.frame(imgs)
            colnames(imgs_df) <- c("LIS", "Coex", "Pre", "Wind", "Tem", "Ecoregion")
            imgs_df <- imgs_df[complete.cases(imgs_df), ]
            Econames <- unique(imgs_df$Ecoregion)
            for (i in 1:length(Econames)) {
                namei <- Econames[i]
                sub_df <- imgs_df[imgs_df$Ecoregion == namei, 1:5]
                if (nrow(sub_df) > 1000) {
                  {
                    library(parallel)
                    result_path <- paste0(outpath)
                    Total_result <- fread(paste0(result_path, namei, "_ROC_Summary.csv"))
                    Total_result <- rbind(Total_result, data.frame(Model = 0, Specificity = 1, mean = 0, 
                      sd = NA), data.frame(Model = 1, Specificity = 1, mean = 0, sd = NA))
                    Total_result$Model <- factor(Total_result$Model, levels = 0:1, labels = c("Without landscape interaction", 
                      "With landscape interaction"))
                    colors <- c("grey45", "red")
                    P_ROC <- ggplot(Total_result, aes(x = 1 - Specificity, y = mean, group = Model)) + 
                      geom_area(data = Total_result[which(Total_result$Model == "With landscape interaction"), ], fill = colors[2], 
                        position = "identity", alpha = 0.7) + geom_area(data = Total_result[which(Total_result$Model == 
                      "Without landscape interaction"), ], fill = "grey90", position = "identity", alpha = 1) + geom_abline(slope = 1, 
                      intercept = 0, color = "grey45", size = 0.2, linetype = 2) + geom_line(aes(color = Model, 
                      group = Model), linewidth = 0.2) + coord_equal() + scale_x_continuous(breaks = seq(0, 
                      1, 0.2), expand = c(0.01, 0.01), limits = c(0, 1)) + scale_y_continuous(breaks = seq(0, 
                      1, 0.2), expand = c(0.01, 0.01), limits = c(0, 1)) + guides(fill = "none", color = "none") + 
                      scale_fill_manual(values = colors) + scale_color_manual(values = colors) + labs(x = "FPR", 
                      y = "TPR") + theme_bw() + theme(plot.title = element_text(size = 10, hjust = 0.5), 
                      legend.text = element_text(size = 10, hjust = 0), legend.title = element_text(size = 10, 
                        face = "bold", hjust = 0), legend.background = element_rect(colour = "NA", fill = "NA", 
                        size = 0), legend.direction = "horizontal", strip.background = element_rect(colour = "NA", 
                        fill = "NA", size = 0), strip.text.x = element_text(size = 10, face = "bold"), 
                      strip.text.y = element_text(size = 10, face = "bold"), axis.title.x = element_text(size = 8, 
                        vjust = 0.5, hjust = 0.5, color = "black", margin = margin(r = 0.15)), axis.title.y = element_text(size = 8, 
                        hjust = 0.5, vjust = 0.5, color = "black", margin = margin(r = 0.2)), axis.text.x = element_blank(), 
                      axis.text.y = element_blank(), axis.ticks = element_line(colour = "black", size = 0.15), 
                      axis.ticks.length = unit(0.01, "cm"), plot.margin = unit(c(0.03, 0.03, 0.03, 0.03), 
                        "cm"), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
                      panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.spacing.x = unit(2, 
                        "mm"), panel.spacing.y = unit(1, "mm"), panel.background = element_rect(colour = "NA", 
                        fill = "NA", size = 0), legend.position = c(0.8, -0.11), plot.background = element_rect(fill = "transparent", 
                        colour = NA, size = 0.02), axis.line = element_line(colour = "black", linewidth = 0.2), 
                      panel.grid = element_blank(), )
                    ggsave(plot = P_ROC, paste0(outpath, "ROC_", namei, ".png"), 
                      width = 2.2, height = 2.1, dpi = 900, units = "cm")
                  }
                }
            }
        }
    }
}

data_path <- paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/")

{
    outpath <- paste0(data_dir, "Summary_AridCutHyper/Results/3_Coexistence/EnhanceModel/Logist/Calculate10km_rep99/")
    clip_roi <- vect(paste0(data_dir, "Image/DrylandRegion/GlobalAridRegion_cutHyperArid.shp"))
    LIS_path <- paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/LIS_10km.tif")
    LIS <- rast(LIS_path) %>% crop(., clip_roi) %>% mask(., clip_roi)
    Coex_10km <- rast(paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/Cooccurrence_10km.tif")) %>% 
        project(., LIS) %>% crop(., LIS, extend = TRUE) %>% mask(., LIS)
    Pre <- rast(paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/Pre_10km.tif")) %>% 
        project(., LIS) %>% crop(., LIS, extend = TRUE) %>% mask(., LIS)
    Tem <- rast(paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/Tem_10km.tif")) %>% 
        project(., LIS) %>% crop(., LIS, extend = TRUE) %>% mask(., LIS)
    Wind <- rast(paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/Wind_10km.tif")) %>% 
        project(., LIS) %>% crop(., LIS, extend = TRUE) %>% mask(., LIS)
    Biome2017 <- rast(paste0(data_path, "Biome2017.tif")) %>% crop(., LIS, extend = TRUE) %>% mask(., 
        LIS)
    Biome2017[Biome2017 == "N/A"] <- NA
    Biome2017 <- mask(Biome2017, Biome2017 %in% c("Deserts & Xeric Shrublands", "Tropical & Subtropical Grasslands, Savannas & Shrublands", 
        "Temperate Grasslands, Savannas & Shrublands", "Mediterranean Forests, Woodlands & Scrub", "Montane Grasslands & Shrublands"), 
        maskvalues = c(NA, 0))
    imgs <- rast(list(LIS, Coex_10km, Pre, Wind, Tem, Biome2017))
    imgs_df <- as.data.frame(imgs)
    colnames(imgs_df) <- c("LIS", "Coex", "Pre", "Wind", "Tem", "Ecoregion")
    imgs_df2 <- imgs_df
    imgs_df2$Ecoregion <- "Global"
    imgs_df <- rbind(imgs_df, imgs_df2)
    imgs_df <- imgs_df[complete.cases(imgs_df), ]
    imgs_df$Ecoregion <- factor(imgs_df$Ecoregion, levels = c("Global", "Deserts & Xeric Shrublands", 
        "Tropical & Subtropical Grasslands, Savannas & Shrublands", "Temperate Grasslands, Savannas & Shrublands", 
        "Mediterranean Forests, Woodlands & Scrub", "Montane Grasslands & Shrublands"), labels = c("Global", 
        "Desert", "Tropical", "Temperate", "Mediterranean", "Montane"))
    Total_results <- lapply(c("Global", "Desert", "Tropical", "Temperate", "Mediterranean", "Montane"), 
        FUN = function(Ecoregion) {
            imgs_df_selected <- imgs_df[imgs_df$Ecoregion == Ecoregion, ]
            imgs_df_scaled <- scale(imgs_df_selected[, c("LIS", "Pre", "Wind", "Tem")]) %>% cbind(data.frame(Coex = imgs_df_selected$Coex), 
                .) %>% as.data.frame()
            m_full <- glm(Coex ~ Pre + I(Pre^2) + Wind + I(Wind^2) + Tem + I(Tem^2) + LIS, data = imgs_df_scaled, 
                family = binomial())
            m_nopre <- glm(Coex ~ Wind + I(Wind^2) + Tem + I(Tem^2) + LIS, data = imgs_df_scaled, family = binomial(link = "logit"))
            m_nowind <- glm(Coex ~ Pre + I(Pre^2) + Tem + I(Tem^2) + LIS, data = imgs_df_scaled, family = binomial(link = "logit"))
            m_notem <- glm(Coex ~ Pre + I(Pre^2) + Wind + I(Wind^2) + LIS, data = imgs_df_scaled, family = binomial(link = "logit"))
            m_nosis <- glm(Coex ~ Pre + I(Pre^2) + Wind + I(Wind^2) + Tem + I(Tem^2), data = imgs_df_scaled, 
                family = binomial(link = "logit"))
            lrt_pre <- anova(m_nopre, m_full, test = "Chisq")
            lrt_wind <- anova(m_nowind, m_full, test = "Chisq")
            lrt_tem <- anova(m_notem, m_full, test = "Chisq")
            lrt_sis <- anova(m_nosis, m_full, test = "Chisq")
            results <- data.frame(Group = c("Precipitation", "Wind", "Temperature", "LIS"), Df = c(lrt_pre$Df[2], 
                lrt_wind$Df[2], lrt_tem$Df[2], lrt_sis$Df[2]), Deviance_Change = c(lrt_pre$Deviance[2], 
                lrt_wind$Deviance[2], lrt_tem$Deviance[2], lrt_sis$Deviance[2]), p_value = c(lrt_pre$`Pr(>Chi)`[2], 
                lrt_wind$`Pr(>Chi)`[2], lrt_tem$`Pr(>Chi)`[2], lrt_sis$`Pr(>Chi)`[2]))
            results$Percent_Contribution <- round(results$Deviance_Change/sum(results$Deviance_Change) * 
                100, 1)
            results$Percent_Label <- paste0(results$Percent_Contribution, "%")
            results$Group <- factor(results$Group, levels = c("LIS", "Precipitation", "Wind", "Temperature"), 
                labels = c("Spatial interaction", "Precipitation parameters", "Wind parameters", "Temperature parameters"))
            library(ggplot2)
            results$Significance <- cut(results$p_value, breaks = c(-Inf, 0.001, 0.01, 0.05, 1), labels = c("***", 
                "**", "*", ""), right = FALSE)
            results$Ecoregion <- Ecoregion
            return(results)
        }) %>% do.call(rbind, .)
    Total_results$Group <- factor(Total_results$Group, levels = c("Spatial interaction", "Precipitation parameters", 
        "Wind parameters", "Temperature parameters"), labels = c("Spatial interaction", "Precipitation", 
        "Wind speed", "Temperature"))
    Total_results$Ecoregion <- factor(Total_results$Ecoregion, levels = c("Global", "Desert", "Tropical", 
        "Temperate", "Mediterranean", "Montane"), labels = c("Global", "Desert", "Tropical", "Temperate", 
        "Mediterranean", "Montane"))
    {
        library(dplyr)
        library(ggplot2)
        plot_df <- Total_results %>% mutate(Ecoregion = factor(Ecoregion, levels = c("Global", "Desert", 
            "Tropical", "Temperate", "Mediterranean", "Montane") %>% rev()), Group = factor(Group, levels = c("Spatial interaction", 
            "Precipitation", "Wind speed", "Temperature"), labels = c("Landscape interaction", "Precipitation", 
            "Wind speed", "Temperature"), )) %>% group_by(Ecoregion) %>% mutate(Percent_Contribution = Percent_Contribution/sum(Percent_Contribution, 
            na.rm = TRUE) * 100) %>% ungroup()
        library(ggplot2)
        library(dplyr)
        library(grid)
        p <- ggplot(plot_df, aes(Percent_Contribution, Ecoregion, fill = Group)) + geom_col(width = 0.75, 
            color = "white", linewidth = 0.3) + scale_x_continuous(limits = c(0, 100), breaks = c(0, 
            25, 50, 75, 100), expand = c(0, 0)) + scale_fill_manual(values = c(`Landscape interaction` = "#C94C4C", 
            Precipitation = "#8FA9BC", `Wind speed` = "#9A8FB3", Temperature = "#CDBA91")) + labs(x = "Relative contribution to explained deviance (%)", 
            y = NULL, fill = NULL) + guides(fill = guide_legend(title = NULL)) + theme_classic(base_size = 8) + theme(axis.title.x = element_text(size = 7.5, 
            color = "black", margin = margin(t = 1)), axis.text.x = element_text(size = 6.5, color = "black"), 
            axis.text.y = element_text(size = 7, color = "black"), axis.ticks = element_line(color = "black", 
                linewidth = 0.25), axis.ticks.length = unit(0.05, "cm"), axis.line = element_line(color = "black", 
                linewidth = 0.35), legend.position = "right", legend.direction = "vertical", legend.justification = "center", 
            legend.text = element_text(size = 7, color = "black"), legend.key = element_blank(), legend.spacing.x = unit(0.12, 
                "cm"), legend.margin = margin(b = 1), panel.border = element_blank(), plot.margin = margin(3, 
                8, 2, 3))
        ggsave(plot = p, paste0(outpath, "Cooccurrence_Predictor_Contributions.jpeg"), 
            width = 8.5, height = 2.9, dpi = 900, units = "cm")
    }
}
