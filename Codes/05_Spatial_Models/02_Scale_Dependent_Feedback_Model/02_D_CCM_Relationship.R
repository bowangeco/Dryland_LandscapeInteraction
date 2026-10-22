# Based on R version 4.1.3.
model_dir <- "05_Spatial_Models/02_Scale_Dependent_Feedback_Model/"
library(dplyr)
library(ggplot2)
library(rEDM)
library(data.table)
library(parallel)
library(Kendall)

source(file.path(model_dir, "01_Define_Model_Functions.R"))
result_path <- file.path(model_dir, "results/CCM/")
dir.create(result_path, recursive = TRUE, showWarnings = FALSE)
D_seq <- c(0.01, 1:8)
rep_n <- 9
cores <- 9

simulation_fun <- function(D) {
    simulation_PDE(LX = 100, LY = 100, NX = 100, NY = 100, dt = 0.001, D = D,
        a = 2, m = 0.45, Dw = 5, meanV = 0.5, Time = 500, noise_sd = 10, noise_prob = 0.1)
}

CCM_fun <- function(result_list) {
    NY <- nrow(result_list[[1]])
    NX <- ncol(result_list[[1]])
    CCM_df <- do.call(rbind, lapply(result_list, as.vector))
    neighbors <- list(
        do.call(rbind, lapply(result_list, function(V) as.vector(V[, c(NX, seq_len(NX - 1))]))),
        do.call(rbind, lapply(result_list, function(V) as.vector(V[c(2:NY, 1), ]))),
        do.call(rbind, lapply(result_list, function(V) as.vector(V[c(NY, seq_len(NY - 1)), ]))),
        do.call(rbind, lapply(result_list, function(V) as.vector(V[, c(2:NX, 1)])))
    )
    n_time <- nrow(CCM_df)
    split_time <- ceiling(n_time/2)
    lib_sizes <- unique(c(10, 20, 30, 50, 70, seq(100, n_time, 50), n_time))
    sampled_cells <- sample(seq_len(ncol(CCM_df)), ncol(CCM_df)/50)
    rho <- vapply(sampled_cells, function(i) {
        # Preserve the reference MAE-based E selection and zero-lag CCM.
        embedding <- rEDM::simplex(CCM_df[, i], E = 2:10,
            lib = c(1, split_time), pred = c(split_time + 1, n_time))
        embedding$mae <- as.numeric(unlist(embedding$mae, use.names = FALSE))
        if (!any(is.finite(embedding$mae))) return(NA_real_)
        best_E <- embedding$E[which.min(replace(embedding$mae, !is.finite(embedding$mae), Inf))]
        directional_rho <- vapply(neighbors, function(neighbor) {
            test_df <- data.frame(time = seq_len(n_time), lib = CCM_df[, i], target = neighbor[, i])
            cmap <- rEDM::CCM(dataFrame = test_df, E = best_E, Tp = 0,
                columns = "lib", target = "target", libSizes = lib_sizes,
                random = FALSE, showPlot = FALSE)
            skills <- cmap[[2]]
            if (length(skills) < 2 || any(!is.finite(skills))) return(NA_real_)
            kend <- Kendall::MannKendall(skills)
            rho_max <- tail(skills, 1)
            rho_min <- skills[1]
            N <- max(cmap$LibSize)
            z <- abs(0.5 * (log((1 + rho_max)/(1 - rho_max)) -
                log((1 + rho_min)/(1 - rho_min))) * (2/(N - 3))^(-0.5))
            z_p <- 1 - pnorm(z)
            significant <- kend$tau[1] > 0 && kend$sl[1] < 0.05 && z_p < 0.05
            if (is.na(significant) || !significant) return(NA_real_)
            rho_max
        }, numeric(1))
        if (all(is.na(directional_rho))) return(NA_real_)
        mean(directional_rho, na.rm = TRUE)
    }, numeric(1))
    data.frame(rho_mean = if (all(is.na(rho))) NA_real_ else mean(rho, na.rm = TRUE),
        rho_sd = sd(rho, na.rm = TRUE), Numb_sign = sum(!is.na(rho)))
}

multifun <- function(D) {
    library(dplyr)
    library(data.table)
    source(file.path(model_dir, "01_Define_Model_Functions.R"))
    set.seed(428)
    repetition_result <- lapply(seq_len(rep_n), function(i) {
        result <- CCM_fun(simulation_fun(D))
        result$D <- D
        result$i <- i
        result
    })
    fwrite(rbindlist(repetition_result), paste0(result_path, "Simulation_D_", D, ".csv"))
    invisible(NULL)
}

run_simulations <- function() {
    cl <- makeCluster(cores)
    on.exit(stopCluster(cl))
    clusterExport(cl, c("model_dir", "result_path", "simulation_fun", "multifun", "CCM_fun", "rep_n"))
    invisible(parLapply(cl, D_seq, multifun))
}
run_simulations()

file_list <- paste0(result_path, "Simulation_D_", D_seq, ".csv")
Total_result <- rbindlist(lapply(file_list, fread))
fwrite(Total_result, paste0(result_path, "LIS_CCM_by_D.csv"))

p <- ggplot(Total_result, aes(D, rho_mean)) +
    geom_smooth(method = "lm", linewidth = 0.5) +
    geom_point(size = 0.01, shape = 20) +
    labs(x = "Model landscape interaction strength, D", y = expression(paste("Cross map skill, ", italic(rho)))) +
    theme_bw() +
    theme(text = element_text(color = "black"), axis.title = element_text(size = 10),
        axis.text = element_text(size = 10, color = "black"),
        axis.ticks = element_line(color = "black", linewidth = 0.25),
        axis.ticks.length = grid::unit(0.1, "cm"),
        plot.margin = grid::unit(rep(0.2, 4), "cm"), panel.grid = element_blank(),
        legend.position = "none")
ggsave(paste0(result_path, "LIS_CCM_by_D.png"),
    plot = p, width = 8, height = 6.3, dpi = 900, units = "cm")
