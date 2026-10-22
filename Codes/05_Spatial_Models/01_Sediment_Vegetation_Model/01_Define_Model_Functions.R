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
    M[which(M > 1)] <- 1
    return(M)
}

noise_mat <- function(sd, prob = 1, NX, NY) {
    nrow <- NY
    ncol <- NX
    M_noise <- rnorm(nrow * ncol, sd = sd) * rbinom(nrow * ncol, 1, prob = prob) %>% matrix(data = ., 
        nrow = nrow, ncol = ncol, dimnames = list(1:nrow, 1:ncol))
    return(M_noise)
}

PDE <- function(V, DP = 400, p = 300, D = 0.1, dx, dy) {
    {
        p_min <- 50
        c <- 100
        a_max <- 0.15
        n <- 0.2
        V_max <- 1
        V_c <- 0.3
        b <- 15
        u <- 0
        e <- 0.001
        r <- 8e-04
        beta <- 0.00067
        k <- 0.0035
        {
            if (p < p_min) {
                a_p <- 0
            }
            else {
                a_p <- a_max * (1 - exp(-(p - p_min)/c))
            }
        }
    }
    g_function <- function(V) {
        return(0.5 * (tanh(b * (V_c - V)) + 1))
    }
    T1 <- a_p * (V + n) * (1 - V/V_max)
    T2 <- -e * DP * g_function(V) * V
    T3 <- -r * DP^(2/3) * V
    T4 <- d2_dxy2(V, dx, dy)
    dV_dt <- T1 + T2 + T3 + T4 * D
    return(dV_dt)
}

randm_01 <- function(meanV, NX, NY) {
    matrix(1 * rbinom(NX * NY, size = 1, prob = meanV), nrow = NY, ncol = NX, dimnames = list(seq_len(NY), 
        seq_len(NX)))
}

simulation_PDE <- function(LX = 100, LY = 100, NX = 100, NY = 100, dt = 0.01, D, meanV, DP, 
    p, Time, noise_sd, noise_prob) {
    dx <- 2 * LX/NX
    dy <- 2 * LY/NY
    V <- randm_01(meanV, NX, NY)
    simulation_result <- list(V)
    if (Time == "Eq") {
        time_max <- 10000
        time_Eq <- 10
        Equilibration <- 0
        for (time in seq(0, time_max, dt)) {
            dV_dt <- PDE(V, DP = DP, p = p, D = D, dx = dx, dy = dy)
            V <- (V + dV_dt * dt) %>% limit_mat()
            if (abs(as.integer(time) - time) < dt * 0.1) {
                simulation_result <- append(simulation_result, list(V))
            }
            if (max(abs(dV_dt), na.rm = TRUE) < 1e-09) {
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
            dV_dt <- PDE(V, DP = DP, p = p, D = D, dx = dx, dy = dy) + noise_mat(noise_sd, 
                noise_prob, NX, NY)/sqrt(dt)
            V <- (V + dV_dt * dt) %>% limit_mat()
            if (abs(as.integer(time) - time) < dt * 0.1) {
                simulation_result <- append(simulation_result, list(V))
            }
        }
    }
    return(simulation_result)
}
