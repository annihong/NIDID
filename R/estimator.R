# Description: This file contains the functions to estimate the parameters of the model

tilde_k_g <- function(df, fixed_res, k,g, df_prob_cluster_g, t_total) {
    treated_df <- df[df$cluster_gmin <= t_total,]
    n_trt <- nrow(treated_df[treated_df$time==0,]) 
    c_trt <- unique(treated_df[treated_df$time==0,"cluster"]) #all treated clusters
    C_trt <- length(c_trt)
    n_c_trt <- fixed_res$n_c[c_trt]
    treated_df$w <- rep(0,nrow(treated_df))
    for (i in 1:nrow(treated_df)) {
        time = treated_df$time[i]
        cluster = treated_df$cluster[i]
        n_c = fixed_res$n_c[cluster]
        prob_gmin_eq_g = df_prob_cluster_g[cluster,as.character(g)]
        gmin_eq_g = treated_df$cluster_gmin[i] == g
        time_eq_g_plus_k = time == g+k
        treated_df$w[i] <-  (1/n_c) * (1/C_trt) * (1/prob_gmin_eq_g) * gmin_eq_g * time_eq_g_plus_k #weights on TE
    }
    return(treated_df)
}

theta_hat_k_g <- function(df, fixed_res, k,g, df_prob_cluster_g_est, t_total){
    treated_df <- df[df$cluster_gmin <= t_total,]
    n_trt <- nrow(treated_df[treated_df$time==0,]) 
    c_trt <- unique(treated_df[treated_df$time==0,"cluster"]) #all treated clusters
    C_trt <- length(c_trt)
    start = g - 1
    end = g + k
    theta_hat <- df[df$time == end, c("Y_obs", "cluster_gmin", "cluster")]
    theta_hat$diff <- theta_hat$Y_obs - df[df$time == start, "Y_obs"]
    theta_hat$w_trt <- rep(0,nrow(theta_hat))
    theta_hat$w_ctrl <- rep(0,nrow(theta_hat))
    theta_hat$w <- rep(0,nrow(theta_hat))
    pmin_est = df_prob_cluster_g_est[as.character(g)]
    pmin_inv_est <- ifelse(pmin_est == 0, 0, 1/pmin_est)
    W_trt <- length(unique(treated_df[treated_df$cluster_gmin == g, "cluster"]))*(1/C_trt)*pmin_inv_est
    C = fixed_res$C
    for (i in 1:nrow(theta_hat)) {
       
        cluster = theta_hat$cluster[i]
        cluster_is_nev_trt <- theta_hat$cluster_gmin[i] > t_total
        gmin_eq_g = theta_hat$cluster_gmin[i] == g
        n_c = fixed_res$n_c[cluster]

        theta_hat$w_trt[i] <-  gmin_eq_g * pmin_inv_est * (1/C_trt) * (1/n_c)
        theta_hat$w_ctrl[i] <-  cluster_is_nev_trt * (1/(C - C_trt)) * (1/n_c) * W_trt
        theta_hat$w[i] <- theta_hat$w_trt[i] - theta_hat$w_ctrl[i]
        # if (cluster == 5) {
        #     print(paste("gmin_eq_g", gmin_eq_g, "pmin_inv_est", pmin_inv_est, "C_trt", C_trt, "n_c", n_c, "theta_hat$w_trt[i]", theta_hat$w_trt[i], "theta_hat$w_ctrl[i]", theta_hat$w_ctrl[i], "theta_hat$w[i]", theta_hat$w[i], "theta_hat$cluster_gmin[i]", theta_hat$cluster_gmin[i], "g", g, "i", i))
        # }
    }
    return(theta_hat)
}
