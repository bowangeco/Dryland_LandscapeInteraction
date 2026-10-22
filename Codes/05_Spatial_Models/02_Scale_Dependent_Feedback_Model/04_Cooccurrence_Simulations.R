# Based on R version 4.1.3.
model_dir <- "05_Spatial_Models/02_Scale_Dependent_Feedback_Model/"

library(dplyr)
library(ggplot2)
library(data.table)
library(parallel)

source(file.path(model_dir, "01_Define_Model_Functions.R"))

simulation_fun <- function(D, meanV, a, m) {
    simulation_PDE(LX = LX, LY = LY, NX = NX, NY = NY, dt = dt, D = D, meanV = meanV, 
        a = a, m = m, Dw = Dw, Time = Time, noise_sd = noise_sd, noise_prob = noise_prob)
}

bimodal_fun <- function(result_list) {
    library(multimode)
    library(diptest)
    Eq_time <- attr(result_list, "stop_time")
    Eq_pattern <- result_list[[length(result_list)]]
    values <- c(Eq_pattern)
    numb <- sum(!is.na(values))
    dip_p <- dip.test(values)$p
    if (length(unique(values)) == 1) {
        ACR_p <- 1
        mod1mean <- NA
        mod1height <- NA
        mod2mean1 <- NA
        mod2mean2 <- NA
        mod2height1 <- NA
        mod2height2 <- NA
        coexistence <- 0
    }
    else {
        if (length(values) > 10000) {
            ACR_p <- multimode::modetest(sample(values, 10000), method = "ACR")$p.value
        }
        else {
            ACR_p <- multimode::modetest(values, method = "ACR")$p.value
        }
        mod1 <- locmodes(data = values, mod0 = 1, lowsup = 0, uppsup = max(values))
        mod2 <- locmodes(data = values, mod0 = 2, lowsup = 0, uppsup = max(values))
        mod1mean <- mod1$locations
        mod1height <- mod1$fvalue
        mod2mean1 <- mod2$locations[1]
        mod2mean2 <- mod2$locations[3]
        mod2height1 <- mod2$fvalue[1]
        mod2height2 <- mod2$fvalue[3]
        coexistence <- (dip_p < 0.05 & ACR_p < 0.05 & mod2height1/mod2height2 > 
            0.05 & mod2height1/mod2height2 < 20)
    }
    result <- data.frame(numb, dip_p, ACR_p, mod1mean, mod1height, mod2mean1, mod2mean2, mod2height1, 
        mod2height2, coexistence, Eq_time)
    colnames(result) <- c("numb", "dip_p", "ACR_p", "mod1mean", "mod1height", "mod2mean1", "mod2mean2", 
        "mod2height1", "mod2height2", "coexistence", "Eq_time")
    return(result)
}

multifun_meanV <- function(meanV) {
    library(dplyr)
    library(data.table)
    source(file.path(model_dir, "01_Define_Model_Functions.R"))
    set.seed(428)
    repetition_result <- data.frame()
    for (i in 1:rep_n) {
        simulation_result <- simulation_fun(D, meanV, a, m)
        coexistence_result <- bimodal_fun(simulation_result)
        coexistence_result$meanV <- meanV
        coexistence_result$D <- D
        coexistence_result$rep_id <- i
        repetition_result <- rbind(repetition_result, coexistence_result)
    }
    fwrite(repetition_result, file = paste0(result_path, "Simulation_D_", D, "_meanV_", meanV, ".csv"), 
        row.names = FALSE)
}

{
    cores <- 9
    cl <- makeCluster(cores)
    clusterExport(cl, c("model_dir"))
    result_path <- file.path(model_dir, "results/Cooccurrence/")
    if (!dir.exists(result_path)) {
        dir.create(result_path, recursive = TRUE)
    }
    D_seq <- c(0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28, 2.56)
    meanV_seq <- seq(0, 1, by = 1/100/100)
    NX <- 100
    NY <- 100
    LX <- 100
    LY <- 100
    dt <- 0.001
    noise_sd <- 0
    noise_prob <- 0
    Time <- "Eq"
    rep_n <- 1
    Dw <- 20
    a <- 2
    m_mean <- 0.45
    set.seed(428)
    m <- rnorm(NX * NY, m_mean, m_mean * 0.3) %>% matrix(., NX, NY) %>% limit_mat()
    for (D in D_seq) {
        clusterExport(cl, varlist = c("result_path", "simulation_fun", "multifun_meanV", "bimodal_fun", 
            "NX", "NY", "LX", "LY", "dt", "noise_sd", "noise_prob", "Time", "a", "m", "Dw", "rep_n", 
            "D"))
        parLapply(cl, X = meanV_seq, fun = multifun_meanV)
        file_list <- paste0(result_path, "Simulation_D_", D, "_meanV_", meanV_seq, ".csv")
        Total_result <- lapply(file_list, fread) %>% do.call(rbind, .)
        fwrite(Total_result, paste0(result_path, "Cooccurrence_D_", D, ".csv"), row.names = FALSE)
        file.remove(file_list)
    }
    stopCluster(cl)
    final_list <- paste0(result_path, "Cooccurrence_D_", D_seq, ".csv")
    All_result <- lapply(final_list, fread) %>% do.call(rbind, .)
    fwrite(All_result, paste0(result_path, "Cooccurrence_Results.csv"), row.names = FALSE)
    file.remove(final_list)
}

{
    plot_df <- as.data.frame(All_result)
    plot_df$coexistence[is.na(plot_df$coexistence)] <- 0
    plot_df$coexistence[plot_df$mod2mean1 > 5 & plot_df$mod2mean2 > 5] <- 0
    plot_df$coexistence[plot_df$mod2mean1 < 5 & plot_df$mod2mean2 < 5] <- 0
    {
        p1 <- ggplot(plot_df, aes(log(D), coexistence)) + geom_jitter(width = 0.2, 
            height = 0.015, alpha = 0.04, size = 0.5) + stat_smooth(geom = "line", method = "glm", method.args = list(family = binomial), 
            formula = "y~x", fullrange = TRUE, fill = "grey85", color = "red", se = TRUE, size = 1, alpha = 0.9) + 
            scale_x_continuous(breaks = log(D_seq), labels = D_seq) + labs(x = "Model landscape interaction strength, D", y = "Co-occurrence proportion") + theme_bw() + theme(plot.title = element_text(size = 10, 
            hjust = 0.5), legend.text = element_text(size = 10, hjust = 0), legend.title = element_text(size = 10, 
            face = "bold", hjust = 0), legend.background = element_rect(colour = "NA", fill = "NA", size = 0), 
            strip.text.x = element_text(size = 10, face = "bold", vjust = 0.3), strip.text.y = element_text(size = 10, 
                face = "plain"), axis.title.x = element_text(size = 10, vjust = 0, color = "black", face = "bold"), 
            axis.title.y = element_text(size = 10, hjust = 0.5, vjust = 1.3, color = "black", face = "bold"), 
            axis.text.x = element_text(size = 8, vjust = 1.5, hjust = 0.5, color = "black"), axis.text.y = element_text(size = 8, 
                color = "black", hjust = 0.5), axis.ticks = element_line(colour = "black", size = 0.25), 
            axis.ticks.length = unit(0.05, "cm"), plot.margin = unit(c(0.2, 0.3, 0.2, 0.2), "cm"), panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), 
            panel.spacing.x = unit(1, "mm"), panel.spacing.y = unit(0.2, "mm"), legend.position = c(0.28, 
                0.85), strip.background = element_blank())
        ggsave(plot = p1, paste0(result_path, "Cooccurrence_by_D.png"), width = 8, height = 6, 
            dpi = 900, units = "cm")
    }
}
