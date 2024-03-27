# Description: This file contains the functions to estimate the parameters of the model

# tilde_k_g <- function(df_long, fixed_res, k,g, df_prob_cluster_g, t_total) {
#     treated_df <- df_long[df_long$cluster_gmin <= t_total,]
#     n_trt <- nrow(treated_df[treated_df$time==0,]) 
#     c_trt <- unique(treated_df[treated_df$time==0,"cluster"]) #all treated clusters
#     marginal_prob_gmin_eq_g_trt = mean(df_prob_cluster_g[c_trt,as.character(g)]/(1 - df_prob_cluster_g[c_trt,as.character(t_total + 1)])) #marginal probability of gmin = g for treated clusters
#     #C_trt <- length(c_trt)
#     n_c_trt <- fixed_res$n_c[c_trt]
#     treated_df$w <- rep(0,nrow(treated_df))
#     for (i in 1:nrow(treated_df)) {
#         time = treated_df$time[i]
#         cluster = treated_df$cluster[i]
#         n_c = fixed_res$n_c[cluster]
#         gmin_eq_g = as.numeric(treated_df$cluster_gmin[i] == g)
#         time_eq_g_plus_k = time == g+k
#         treated_df$w[i] <-  (1/n_c) * (n_c/n_trt) * (1/marginal_prob_gmin_eq_g_trt) * gmin_eq_g * time_eq_g_plus_k #weights on TE
#     }
#     return(treated_df)
# }
#' @export
theta_tilde_k <- function(df_long, fixed_res, k, df_prob_cluster_g, t_total, result=NULL) {
    g_seq=1:t_total
    treated_df <- df_long[df_long$cluster_gmin <= t_total,]
    n_trt <- nrow(treated_df[treated_df$time==0,]) 
    c_trt <- unique(treated_df[treated_df$time==0,"cluster"]) #all treated clusters
    marginal_prob_gmin_trt = sapply(1:t_total, function(g) mean(df_prob_cluster_g[c_trt,as.character(g)]/(1 - df_prob_cluster_g[c_trt,as.character(t_total + 1)]))) #marginal probability of gmin = g conditional on being treated
    theta_tilde_k_g_ws <- lapply(g_seq, function(g) theta_tilde_k_g(df_long, fixed_res, k, g, marginal_prob_gmin_trt[g], t_total, n_trt))
    theta_tilde_k_g_s <- sapply(theta_tilde_k_g_ws, function(w) sum(w * df_long$TE))
    df_long$theta_tilde_k_g_w_TE <- rep(0,nrow(df_long))
    for (g in g_seq) { #Reduce("+", lapply(theta_tilde_k_g_ws, function(w) w * df_long$TE))
        df_long$theta_tilde_k_g_w_TE <- df_long$theta_tilde_k_g_w_TE + theta_tilde_k_g_ws[[g]] * df_long$TE * marginal_prob_gmin_trt[g]
    }
    theta_tilde_k_g_w_TE_ci <- lapply(0:t_total, function(t) df_long[df_long$time == t, "theta_tilde_k_g_w_TE"])
    theta_tilde_k_ci <-  Reduce("+", theta_tilde_k_g_w_TE_ci)
    theta_tilde_k <- sum(theta_tilde_k_g_s * marginal_prob_gmin_trt)
    res <- list(theta_tilde_k = theta_tilde_k, theta_tilde_k_g_s = theta_tilde_k_g_s, theta_tilde_k_ci = theta_tilde_k_ci)

    if (!is.null(result)) {
        return(res[[result]])
    } else {
        return(res)
    }
}
#' @export
theta_tilde_k_g <- function(df_long, fixed_res, k,g, marginal_prob_gmin_trt_g, t_total, n_trt) {
    w <- rep(0,nrow(df_long)) #weight on all points, but the w = 0 if not treated
    for (i in 1:nrow(df_long)) {
        cluster_gmin = df_long$cluster_gmin[i]
        time = df_long$time[i]
        cluster = df_long$cluster[i]
        n_c = fixed_res$n_c[cluster]
        gmin_eq_g = as.numeric(cluster_gmin == g)
        cluster_is_nev_trt <- as.numeric(cluster_gmin > t_total)
        time_eq_g_plus_k = time == g+k
        w[i] <-  (1/n_c) * (n_c/n_trt) * (1/marginal_prob_gmin_trt_g) * gmin_eq_g * time_eq_g_plus_k #weights on TE
    }
    return(w)
}
#' @export
theta_hat_k <- function(df_obs, fixed_res, k, gmin_frequency, t_total, result=NULL){
    g_seq <- 1:(t_total - k) #we won't be able to estimate theta for g > t_total - k
    treated_df <- df_obs[df_obs$cluster_gmin <= t_total & df_obs$time == 0,]
    n_trt <- nrow(treated_df) 
    c_trt <- unique(treated_df[,"cluster"]) #all treated clusters
    C_trt <- length(c_trt)

    conditional_probs <- (gmin_frequency/C_trt)[g_seq] #estimated marginal probability of gmin = g conditional on being treated
    theta_hat_k_g_w_diffs <- lapply(g_seq, function(g) theta_hat_k_g(df_obs, fixed_res, k, g,t_total, conditional_probs[g], n_trt, c_trt, C_trt))
    theta_hat_k_g_s <- sapply(theta_hat_k_g_w_diffs, sum)
    theta_hat_k_ci <- rep(0, fixed_res$n)
    for (g in g_seq) {
        theta_hat_k_ci <- theta_hat_k_ci + theta_hat_k_g_w_diffs[[g]] * conditional_probs[g]
    }
    theta_hat_k <- sum(theta_hat_k_ci)
    res <- list(theta_hat_k = theta_hat_k, theta_hat_k_g_s = unlist(theta_hat_k_g_s), theta_hat_k_ci = unlist(theta_hat_k_ci))

    if (!is.null(result)) {
        return(res[[result]])
    } else {
        return(res)
    }
}

#' @export
theta_hat_k_g <- function(df_obs, fixed_res, k,g, t_total, gmin_prob_est, n_trt, c_trt, C_trt){
    C = fixed_res$C
    start = g - 1
    end = g + k
    theta_hat <- df_obs[df_obs$time == end, c("Y_obs", "cluster_gmin", "cluster")] #include i for treated and control clusters
    theta_hat$diff <- theta_hat$Y_obs - df_obs[df_obs$time == start, "Y_obs"]
    theta_hat$w_trt <- rep(0,nrow(theta_hat)) # = 0 if i is not treated
    theta_hat$w_ctrl <- rep(0,nrow(theta_hat)) # = 0 if i is treated
    theta_hat$w <- rep(0,nrow(theta_hat)) # total weight
    theta_hat$n_c <- fixed_res$n_c[theta_hat$cluster]

    gmin_inv_est <-  ifelse(gmin_prob_est > 0, 1/gmin_prob_est, 0)
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
    }
    return(theta_hat$w * theta_hat$diff)
}
#' @export
omega_tilde_n <- function(theta_hat_ci, theta_tilde_ci, fixed_res){
    n <- fixed_res$n
    C <- fixed_res$C
    df <- data.frame(theta_hat_ci = theta_hat_ci, theta_tilde_ci = theta_tilde_ci, cluster=fixed_res$cluster)
    res <- 0
    for (c in 1:C) {
        df_c <- df[df$cluster == c,]
        sum_sq <- (sum(df_c$theta_hat_ci*n - df_c$theta_tilde_ci*n))^2
        res <- res + sum_sq
        # if (sum_sq < 1e-8) {
        # print(df_c)
        # }
    }
    res <- res / n
    return(res)
}
#' @export
omega_hat_n <- function(theta_hat_ci, fixed_res){
    n <- fixed_res$n
    df <- data.frame(theta_hat_ci = theta_hat_ci, cluster=fixed_res$cluster)
    theta_bar <- sum(df$theta_hat_ci)
    n <- fixed_res$n
    C <- fixed_res$C
    res <- 0
    for (c in 1:C) {
        df_c <- df[df$cluster == c,]
        sum_sq <- (sum(df_c$theta_hat_ci*n - theta_bar))^2
        res <- res + sum_sq

    }
    res <- res / n
    return(res)
}

