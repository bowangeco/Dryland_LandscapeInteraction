# Based on R version 4.1.3.
model_dir <- "05_Spatial_Models/01_Sediment_Vegetation_Model/"

library(dplyr)
library(ggplot2)
library(data.table)
library(parallel)

source(file.path(model_dir, "01_Define_Model_Functions.R"))
dir.create(file.path(model_dir, "results/Recovery/Plots/"),
    recursive = TRUE, showWarnings = FALSE)
DP <- 600
p <- 300
dt <- 0.01
noise_sd <- 1
noise_prob <- 0.1
rep_n <- 9
meanV <- 0.5

simulation_fun <- function(D, meanV, DP, p) {
    NX <- 100
    NY <- 100
    LX <- 100
    LY <- 100
    Time <- 500
    simulation_PDE(LX = LX, LY = LY, NX = NX, NY = NY, dt = dt, D = D, meanV = meanV, 
        DP = DP, p = p, Time = Time, noise_sd = noise_sd, noise_prob = noise_prob)
}

multifun <- function(D) {
    library(dplyr)
    library(data.table)
    source(file.path(model_dir, "01_Define_Model_Functions.R"))
    set.seed(428)
    repetition_result <- data.frame()
    for (i in 1:rep_n) {
        result_df <- simulation_fun(D, meanV, DP, p)
        meanV_df <- lapply(result_df, FUN = function(list_i) {
            mean(list_i, na.rm = TRUE)
        }) %>% unlist()
        if (nrow(repetition_result) == 0) {
            repetition_result <- data.frame(meanV = meanV_df)
        }
        else {
            repetition_result <- cbind(repetition_result, data.frame(meanV = meanV_df))
        }
    }
    repetition_result$D <- D
    repetition_result$time <- 1:nrow(repetition_result) - 1
    fwrite(repetition_result, file = paste0(result_path, "Simulation_D_", D, ".csv"), row.names = FALSE)
}

{
    cores <- 9
    cl <- makeCluster(cores)
    clusterExport(cl, c("model_dir"))
    result_path <- file.path(model_dir, "results/Recovery/Fixed_D/")
    dir.create(result_path, recursive = TRUE, showWarnings = FALSE)
    D_seq <- c(0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28, 2.56)
    clusterExport(cl, varlist = c("result_path", "simulation_fun", "multifun", "dt", "noise_sd", "noise_prob", 
        "meanV", "DP", "p", "rep_n"))
    parLapply(cl = cl, X = D_seq, fun = multifun)
    stopCluster(cl)
    file_list <- paste0(result_path, "Simulation_D_", D_seq, ".csv")
    DP_id <- sub("^.*Simulation_D_(.*)\\.csv$", "\\1", file_list) %>% as.numeric()
    file_list <- file_list[order(DP_id)]
    Total_result <- lapply(file_list, FUN = function(p) {
        fread(p) %>% data.frame()
    }) %>% do.call(rbind, .)
    fwrite(Total_result, paste0(result_path, "Recovery_Fixed_D_Series.csv"), row.names = FALSE)
    file.remove(file_list)
}

{
    Total_result <- fread(paste0(result_path, "Recovery_Fixed_D_Series.csv")) %>% as.data.frame()
}

{
    calc_ar1 <- function(x) {
        ac1 <- acf(x, lag.max = 1, plot = FALSE)$acf[2]
        return(ac1)
    }
    Lambda_df <- lapply(unique(Total_result$D), FUN = function(D) {
        result_D <- Total_result[which(Total_result$D == D), ]
        ts_df <- result_D[, 1:(ncol(result_D) - 2)]
        rdf_D <- lapply(1:ncol(ts_df), FUN = function(col_id) {
            x <- ts_df[, col_id] %>% c() %>% tail(., -100)
            dt <- 1
            if (anyNA(x) == TRUE) {
                lambda_ac1 <- NA
            }
            else {
                ac1 <- calc_ar1(x)
                if (ac1 <= 0) {
                  lambda_ac1 <- NA
                }
                else {
                  lambda_ac1 <- log(ac1)/dt
                }
            }
            rdf_D_col <- data.frame(lambda_ac1 = lambda_ac1, D = D)
            return(rdf_D_col)
        }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .)
    p <- ggplot(Lambda_df, aes(x = as.factor(D), y = -lambda_ac1)) + stat_boxplot(geom = "errorbar", 
        width = 0.35, position = position_dodge(0.8), color = "black") + geom_boxplot(outlier.shape = 16, 
        fill = "#7091CB") + labs(x = "Model landscape interaction strength, D", y = "Recovery rate") + 
        theme_bw() + theme(plot.title = element_text(size = 10, hjust = 0.5), legend.text = element_text(size = 10, 
        hjust = 0), legend.title = element_text(size = 10, face = "bold", hjust = 0), legend.background = element_rect(colour = "NA", 
        fill = "NA", size = 0), strip.text.x = element_text(size = 10, face = "bold", vjust = 0.3), strip.text.y = element_text(size = 10, 
        face = "plain"), axis.title.x = element_text(size = 10, vjust = 0, color = "black", face = "bold"), 
        axis.title.y = element_text(size = 10, hjust = 0.5, vjust = 1.3, color = "black", face = "bold"), 
        axis.text.x = element_text(size = 8, vjust = 1.5, hjust = 0.5, color = "black"), axis.text.y = element_text(size = 8, 
            color = "black", hjust = 0.5), axis.ticks = element_line(colour = "black", size = 0.25), 
        axis.ticks.length = unit(0.05, "cm"), plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"), panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), 
        panel.spacing.x = unit(1, "mm"), panel.spacing.y = unit(0.2, "mm"), legend.position = c(0.28, 
            0.85), strip.background = element_blank())
    ggsave(plot = p, file.path(model_dir, "results/Recovery/Plots/Recovery_Rate_by_D.jpeg"), 
        width = 8, height = 6, dpi = 900, units = "cm")
}

{
    simulation_PDE_changingD <- function(LX = 100, LY = 100, NX = 100, NY = 100, dt = 0.01, D, meanV, 
        DP, p, Time, noise_sd, noise_prob) {
        dx <- 2 * LX/NX
        dy <- 2 * LY/NY
        V <- randm_01(meanV, NX, NY)
        simulation_result <- list(V)
        if (class(Time) == "numeric") {
            for (time in seq(0, Time, dt)) {
                dV_dt <- PDE(V, DP = DP, p = p, D = D, dx = dx, dy = dy) + noise_mat(noise_sd, 
                  noise_prob, NX, NY)/sqrt(dt)
                V <- (V + dV_dt * dt) %>% limit_mat()
                if (abs(as.integer(time) - time) < dt * 0.1) {
                  simulation_result <- append(simulation_result, list(V))
                }
                D_change_speed <- 0.0002
                if (time > 2000) {
                    D <- D + D_change_speed * dt
                }
            }
        }
        return(simulation_result)
    }
    set.seed(428)
    DP <- 600
    p <- 300
    dt <- 0.01
    noise_sd <- 1
    noise_prob <- 0.1
    rep_n <- 1
    meanV <- 0.5
    simulation_fun <- function(D, meanV, DP, p) {
        NX <- 100
        NY <- 100
        LX <- 100
        LY <- 100
        dt <- 0.01
        Time <- 12000
        simulation_PDE_changingD(LX = LX, LY = LY, NX = NX, NY = NY, dt = dt, D = D, meanV = meanV, 
            DP = DP, p = p, Time = Time, noise_sd = noise_sd, noise_prob = noise_prob)
    }
    repetition_result <- data.frame()
    for (i in 1:rep_n) {
        result_df <- simulation_fun(D = 0, meanV, DP, p)
        meanV_df <- lapply(result_df, FUN = function(list_i) {
            mean(list_i, na.rm = TRUE)
        }) %>% unlist()
        if (nrow(repetition_result) == 0) {
            repetition_result <- data.frame(meanV = meanV_df)
        }
        else {
            repetition_result <- cbind(repetition_result, data.frame(meanV = meanV_df))
        }
    }
    repetition_result$time <- 1:nrow(repetition_result) - 1
    result_path <- file.path(model_dir, "results/Recovery/Dynamic_D/")
    dir.create(result_path, recursive = TRUE, showWarnings = FALSE)
    fwrite(repetition_result, paste0(result_path, "Recovery_Dynamic_D_Series.csv"))
    Total_result <- repetition_result
    Total_result$Type <- Total_result$time >= 2000
    {
        calc_ar1 <- function(x) {
            ac1 <- acf(x, lag.max = 1, plot = FALSE)$acf[2]
            return(ac1)
        }
        win_size <- 100
        Lambda_df <- lapply(unique(as.integer(round(seq(from = 0, to = nrow(Total_result) - win_size, length.out = 200)))), 
            FUN = function(time_win) {
                result_win <- Total_result[time_win + 1:win_size, ]
                ts_df <- result_win[, setdiff(names(result_win), c("time", "Type")), drop = FALSE]
                rdf_D <- lapply(1:ncol(ts_df), FUN = function(col_id) {
                  x <- select(ts_df, col_id) %>% unlist() %>% unname()
                  dt <- 1
                  if (anyNA(x) == TRUE) {
                    lambda_ac1 <- NA
                  }
                  else {
                    ac1 <- calc_ar1(x)
                    if (ac1 <= 0) {
                      lambda_ac1 <- NA
                    }
                    else {
                      lambda_ac1 <- log(ac1)/dt
                    }
                  }
                  rdf_D_col <- data.frame(lambda_ac1 = lambda_ac1, time = time_win)
                  return(rdf_D_col)
                }) %>% do.call(rbind, .)
            }) %>% do.call(rbind, .)
        {
            trend_data <- subset(Lambda_df, time >= 2000 & is.finite(lambda_ac1))
            trend_data$plot_time <- trend_data$time - 2000
            trend_model <- lm(lambda_ac1 ~ plot_time, data = trend_data)
            trend_curve <- data.frame(plot_time = seq(min(trend_data$plot_time),
                max(trend_data$plot_time), length.out = 80))
            trend_prediction <- predict(trend_model, newdata = trend_curve, interval = "confidence")
            trend_curve$recovery <- -trend_prediction[, "fit"]
            trend_curve$lower <- -trend_prediction[, "upr"]
            trend_curve$upper <- -trend_prediction[, "lwr"]
            p <- ggplot(trend_data, aes(plot_time, -lambda_ac1)) + geom_point(size = 0.5) + 
                geom_ribbon(data = trend_curve, aes(x = plot_time, ymin = lower, ymax = upper),
                    inherit.aes = FALSE, fill = "red", alpha = 0.3) +
                geom_line(data = trend_curve, aes(x = plot_time, y = recovery),
                    inherit.aes = FALSE, color = "red", linewidth = 0.5) + labs(x = "Time (with increasing D)", 
                y = "Recovery rate") + scale_x_continuous(limits = c(0, 10000),
                breaks = seq(0, 10000, 2500)) + theme_bw() + theme(plot.title = element_text(size = 10, hjust = 0.5), 
                legend.text = element_text(size = 10, hjust = 0), legend.title = element_text(size = 10, 
                  face = "bold", hjust = 0), legend.background = element_rect(colour = "NA", fill = "NA", 
                  size = 0), strip.text.x = element_text(size = 10, face = "bold", vjust = 0.3), strip.text.y = element_text(size = 10, 
                  face = "plain"), axis.title.x = element_text(size = 10, vjust = 0, color = "black", 
                  face = "bold"), axis.title.y = element_text(size = 10, hjust = 0.5, vjust = 1.3, color = "black", 
                  face = "bold"), axis.text.x = element_text(size = 8, vjust = 1.5, hjust = 0.5, color = "black"), 
                axis.text.y = element_text(size = 8, color = "black", hjust = 0.5), axis.ticks = element_line(colour = "black", 
                  size = 0.25), axis.ticks.length = unit(0.05, "cm"), plot.margin = unit(c(0.2, 0.3, 
                  0.2, 0.2), "cm"), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
                panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.spacing.x = unit(1, 
                  "mm"), panel.spacing.y = unit(0.2, "mm"), legend.position = c(0.28, 0.85), strip.background = element_blank())
            ggsave(plot = p, file.path(model_dir, "results/Recovery/Plots/Recovery_Rate_Over_Time.jpeg"), 
                width = 8, height = 6, dpi = 900, units = "cm")
        }
    }
}
