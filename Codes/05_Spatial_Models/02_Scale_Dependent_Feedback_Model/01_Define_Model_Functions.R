# Based on R version 4.1.3.
library(dplyr)

d2_dxy2 <- function(M, dx, dy) {
    nrow <- nrow(M)
    ncol <- ncol(M)
    Mleft <- M[, c(ncol, 1:(ncol - 1))]
    Mright <- M[, c(2:ncol, 1)]
    Mup <- M[c(nrow, 1:(nrow - 1)), ]
    Mdown <- M[c(2:nrow, 1), ]
    M2 <- (Mleft + Mright - 2 * M)/dx^2 + (Mup + Mdown - 2 * M)/dy^2
    return(M2)
}

limit_mat <- function(M) {
    M[which(M < 0)] <- 0
    return(M)
}

mean_mat <- function(M) {
    m_val <- mean(M, na.rm = TRUE)
    M2 <- matrix(m_val, nrow = nrow(M), ncol = ncol(M))
    return(M2)
}

noise_mat <- function(sd, prob = 1, NX, NY) {
    nrow <- NY
    ncol <- NX
    M_noise <- rnorm(nrow * ncol, sd = sd) * rbinom(nrow * ncol, 1, prob = prob) %>% matrix(data = ., 
        nrow = nrow, ncol = ncol, dimnames = list(1:nrow, 1:ncol))
    return(M_noise)
}

d_k <- function(M, dx) {
    nrow <- nrow(M)
    ncol <- ncol(M)
    Mleft <- matrix(NA, nrow, ncol, dimnames = list(1:nrow, 1:ncol))
    Mleft[, -1] <- M[, -ncol]
    Mdif <- (M - Mleft)/dx^2
    return(Mdif)
}

PDE <- function(w, p, dx, dy, a = 2, m = 0.45, Dw, D = 1) {
    Laplace_w <- d2_dxy2(w, dx, dy)
    Laplace_p <- d2_dxy2(p, dx, dy)
    dwdt <- a - w - w * p^2 + Dw * Laplace_w
    dpdt <- w * p^2 - m * p + D * Laplace_p
    return(list(dwdt = dwdt, dpdt = dpdt))
}

randm_01 <- function(meanV, NX, NY) {
    matrix(10 * rbinom(NX * NY, size = 1, prob = meanV), nrow = NY, ncol = NX, dimnames = list(seq_len(NY), 
        seq_len(NX)))
}

simulation_PDE <- function(LX = 100, LY = 100, NX = 100, NY = 100, dt = 0.001, a, m, Dw, D, meanV, 
    Time, noise_sd, noise_prob) {
    dx <- 2 * LX/NX
    dy <- 2 * LY/NY
    initial_w <- matrix(runif(NX * NY), nrow = NX, ncol = NY)
    w <- initial_w
    p <- randm_01(meanV, NX, NY)
    simulation_result_w <- list(w)
    simulation_result_p <- list(p)
    if (Time == "Eq") {
        time_max <- 10000
        time_Eq <- 10
        Equilibration <- 0
        for (time in seq(0, time_max, dt)) {
            PDE_result <- PDE(w, p, dx, dy, a, m, Dw, D)
            dwdt <- PDE_result$dwdt
            dpdt <- PDE_result$dpdt
            w <- w + dwdt * dt
            p <- (p + dpdt * dt) %>% limit_mat()
            if (abs(as.integer(time) - time) < dt * 0.1) {
                simulation_result_p <- append(simulation_result_p, list(p))
            }
            max_change <- max(abs(dwdt), abs(dpdt), na.rm = TRUE)
            if (max_change < 1e-06) {
                Equilibration <- Equilibration + 1
            }
            else {
                Equilibration <- 0
            }
            if (Equilibration >= time_Eq/dt) {
                break
            }
        }
    }
    if (class(Time) == "numeric") {
        for (time in seq(0, Time, dt)) {
            PDE_result <- PDE(w, p, dx, dy, a, m, Dw, D)
            dwdt <- PDE_result$dwdt
            dpdt <- PDE_result$dpdt + noise_mat(noise_sd, noise_prob, NX, NY)/sqrt(dt)
            w <- (w + dwdt * dt)
            p <- (p + dpdt * dt) %>% limit_mat()
            if (abs(as.integer(time) - time) < dt * 0.1) {
                simulation_result_p <- append(simulation_result_p, list(p))
            }
        }
    }
    if (Time == "Eq") attr(simulation_result_p, "stop_time") <- time
    return(simulation_result_p)
}
