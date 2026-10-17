# Based on R version 4.1.3.
model_dir <- "05_Spatial_Models/01_Sediment_Vegetation_Model/"

library(dplyr)
library(ggplot2)
library(data.table)
library(parallel)

source(file.path(model_dir, "01_Define_Model_Functions.R"))

simulation_fun <- function(D, meanV, DP, p) {
    simulation_PDE(LX = LX, LY = LY, NX = NX, NY = NY, dt = dt, D = D, meanV = meanV, 
        DP = DP, p = p, Time = Time, noise_sd = noise_sd, noise_prob = 1)
}

bimodal_fun <- function(result_list) {
    library(multimode)
    library(diptest)
    Eq_time <- length(result_list) * dt
    Eq_pattern <- result_list[[length(result_list)]]
    values <- c(Eq_pattern)
    numb <- sum(!is.na(values))
    dip_p <- dip.test(values)$p
    if (length(unique(values)) == 1L) {
        return(data.frame(numb = numb, dip_p = 1, ACR_p = 1, mod1mean = values[1], mod1height = NA_real_, 
            mod2mean1 = NA_real_, mod2mean2 = NA_real_, mod2height1 = NA_real_, mod2height2 = NA_real_, 
            coexistence = FALSE, Eq_time = Eq_time))
    }
    ACR_values <- if (length(values) > 10000) sample(values, 10000) else values
    ACR_p <- multimode::modetest(ACR_values, method = "ACR")$p.value
    mod1 <- locmodes(data = values, mod0 = 1, lowsup = 0, uppsup = 1)
    mod2 <- locmodes(data = values, mod0 = 2, lowsup = 0, uppsup = 1)
    mod1mean <- mod1$locations
    mod1height <- mod1$fvalue
    mod2mean1 <- mod2$locations[1]
    mod2mean2 <- mod2$locations[3]
    mod2height1 <- mod2$fvalue[1]
    mod2height2 <- mod2$fvalue[3]
    coexistence <- (dip_p < 0.05 & ACR_p < 0.05 & mod2height1/mod2height2 > 
        0.05 & mod2height1/mod2height2 < 20)
    result <- data.frame(numb, dip_p, ACR_p, mod1mean, mod1height, mod2mean1, mod2mean2, mod2height1, 
        mod2height2, coexistence, Eq_time)
    colnames(result) <- c("numb", "dip_p", "ACR_p", "mod1mean", "mod1height", "mod2mean1", "mod2mean2", 
        "mod2height1", "mod2height2", "coexistence", "Eq_time")
    return(result)
}

multifun <- function(D) {
    library(dplyr)
    library(data.table)
    source(file.path(model_dir, "01_Define_Model_Functions.R"))
    set.seed(428)
    repetition_result <- data.frame()
    for (i in 1:rep_n) {
        meanV_seq <- seq(0, 1, by = 1/NY/NX)
        result_df <- lapply(meanV_seq, FUN = function(meanV) {
            simulation_result <- simulation_fun(D, meanV, DP, p)
            coexistence_result <- bimodal_fun(simulation_result)
            coexistence_result$meanV <- meanV
            return(coexistence_result)
        }) %>% do.call(rbind, .)
        result_df$rep_id <- i
        result_df$D <- D
        result_df$DP <- DP
        repetition_result <- rbind(repetition_result, result_df)
    }
    fwrite(repetition_result, file = paste0(result_path, "Simulation_D_", D, ".csv"), row.names = FALSE)
}

{
    cores <- 6
    cl <- makeCluster(cores)
    clusterExport(cl, c("model_dir"))
    result_path <- file.path(model_dir, "results/Cooccurrence/")
    dir.create(result_path, recursive = TRUE, showWarnings = FALSE)
    D_seq <- c(0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28, 2.56)
    DP <- 600
    p <- 300
    NX <- 100
    NY <- 100
    LX <- 100
    LY <- 100
    dt <- 0.01
    noise_sd <- 0
    Time <- "Eq"
    rep_n <- 1
    clusterExport(cl, varlist = c("result_path", "simulation_fun", "multifun", "bimodal_fun", "NX", "NY", 
        "LX", "LY", "dt", "noise_sd", "Time", "DP", "p", "rep_n"))
    parLapply(cl = cl, X = D_seq, fun = multifun)
    stopCluster(cl)
    file_list <- paste0(result_path, "Simulation_D_", D_seq, ".csv")
    DP_id <- sub("^.*Simulation_D_(.*)\\.csv$", "\\1", file_list) %>% as.numeric()
    file_list <- file_list[order(DP_id)]
    Total_result <- lapply(file_list, FUN = function(p) {
        fread(p) %>% data.frame()
    }) %>% do.call(rbind, .)
    fwrite(Total_result, paste0(result_path, "Cooccurrence_Results.csv"), row.names = FALSE)
    file.remove(file_list)
}

{
    result_path <- file.path(model_dir, "results/Cooccurrence/")
    plot_df <- as.data.frame(Total_result)
    plot_df$coexistence <- 0
    plot_df$coexistence[which(plot_df$dip_p < 0.05 & plot_df$ACR_p < 0.05 & plot_df$mod2mean1 < 0.5 & 
        plot_df$mod2mean2 > 0.5)] <- 1
    plot_df$coexistence <- as.numeric(plot_df$coexistence)
    {
        p1 <- ggplot(plot_df, aes(log(D), coexistence)) + geom_jitter(width = 0.2, 
            height = 0.015, alpha = 0.02, size = 0.5) + stat_smooth(geom = "line", method = "glm", method.args = list(family = binomial), 
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
