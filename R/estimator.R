# Description: This file contains the functions to estimate the parameters of the model

tilde_k_g <- function(df_long, fixed_res, k,g, df_prob_cluster_g, t_total) {
    treated_df <- df_long[df_long$cluster_gmin <= t_total,]
    n_trt <- nrow(treated_df[treated_df$time==0,]) 
    c_trt <- unique(treated_df[treated_df$time==0,"cluster"]) #all treated clusters
    marginal_prob_gmin_eq_g_trt = mean(df_prob_cluster_g[c_trt,as.character(g)]/(1 - df_prob_cluster_g[c_trt,as.character(t_total + 1)])) #marginal probability of gmin = g for treated clusters
    #C_trt <- length(c_trt)
    n_c_trt <- fixed_res$n_c[c_trt]
    treated_df$w <- rep(0,nrow(treated_df))
    for (i in 1:nrow(treated_df)) {
        time = treated_df$time[i]
        cluster = treated_df$cluster[i]
        n_c = fixed_res$n_c[cluster]
        gmin_eq_g = as.numeric(treated_df$cluster_gmin[i] == g)
        time_eq_g_plus_k = time == g+k
        treated_df$w[i] <-  (1/n_c) * (n_c/n_trt) * (1/marginal_prob_gmin_eq_g_trt) * gmin_eq_g * time_eq_g_plus_k #weights on TE
    }
    return(treated_df)
}

theta_hat_k_g <- function(df_obs, fixed_res, k,g, gmin_frequency, t_total){
    treated_df <- df_obs[df_obs$cluster_gmin <= t_total & df_obs$time == 0,]
    n_trt <- nrow(treated_df) 
    c_trt <- unique(treated_df[,"cluster"]) #all treated clusters
    C_trt <- length(c_trt)
    C = fixed_res$C
    start = g - 1
    end = g + k
    theta_hat <- df_obs[df_obs$time == end, c("Y_obs", "cluster_gmin", "cluster")] #include i for treated and control clusters
    theta_hat$diff <- theta_hat$Y_obs - df_obs[df_obs$time == start, "Y_obs"]
    theta_hat$w_trt <- rep(0,nrow(theta_hat)) # = 0 if i is not treated
    theta_hat$w_ctrl <- rep(0,nrow(theta_hat)) # = 0 if i is treated
    theta_hat$w <- rep(0,nrow(theta_hat)) # total weight
    theta_hat$n_c <- fixed_res$n_c[theta_hat$cluster]

    contional_prob <- gmin_frequency/C_trt #estimated marginal probability of gmin = g conditional on being treated
    gmin_prob_est = contional_prob[as.character(g)]
    gmin_inv_est <- ifelse(gmin_prob_est == 0, 0, 1/gmin_prob_est)
    H_weight <- gmin_inv_est * sum(theta_hat$cluster_gmin == g) * 1/n_trt
    #H_weight <- gmin_inv_est * mean(theta_hat$cluster_gmin == g)

    for (i in 1:nrow(theta_hat)) {
       
        cluster = theta_hat$cluster[i]
        cluster_is_nev_trt <- as.numeric(theta_hat$cluster_gmin[i] > t_total)
        gmin_eq_g = as.numeric(theta_hat$cluster_gmin[i] == g)
        n_c = theta_hat$n_c[i]

        theta_hat$w_trt[i] <-  (1 - cluster_is_nev_trt) * (gmin_eq_g) * gmin_inv_est * 1/n_trt
        theta_hat$w_ctrl[i] <-  cluster_is_nev_trt * (1/(fixed_res$n - n_trt)) * H_weight
        theta_hat$w[i] <- theta_hat$w_trt[i] - theta_hat$w_ctrl[i]
        # if (cluster == 5) {
        #     print(paste("gmin_eq_g", gmin_eq_g, "gmin_inv_est", gmin_inv_est, "C_trt", C_trt, "n_c", n_c, "theta_hat$w_trt[i]", theta_hat$w_trt[i], "theta_hat$w_ctrl[i]", theta_hat$w_ctrl[i], "theta_hat$w[i]", theta_hat$w[i], "theta_hat$cluster_gmin[i]", theta_hat$cluster_gmin[i], "g", g, "i", i))
        # }
    }
    return(theta_hat)
}
