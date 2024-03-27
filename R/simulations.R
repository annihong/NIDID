
#helper functions:

#' Generate individual term based on cluster term and normal distribution
#'
#' This function generates an individual term based on a cluster term and a normal distribution.
#' The cluster term is selected from `mu_c` based on the index `c`, and the individual term is generated
#' by drawing from a normal distribution with mean equal to the cluster term and standard deviation 5.
#' The number of individual terms generated is determined by `n_c[c]`.
#'
#' @param c An index for selecting the cluster term from `mu_c` and the number of terms from `n_c`.
#' @param mu_c A vector of cluster terms.
#' @param n_c A vector of the number of individual terms to generate for each cluster.
#'
#' @return A vector of individual terms.
#' @export
helper_indi_term_c <-  function(c, mu_c, n_c){
        cluster_term = mu_c[c]
        indi_term = rnorm(n_c[c], cluster_term, 5)
        return(indi_term)
}
#functions to help create simulated data:
#1 fixed value simulation (only sim once)

#' Simulate fixed values for a given number of clusters
#'
#' This function simulates fixed values for a given number of clusters (C). It generates
#' the overall cluster mean, number of observations in each cluster, total number of observations,
#' probability of being untreated on the cluster level, cluster level treatment assignment, individual term, cluster,
#' never treated, and treatment.
#'
#' @param C An integer. The number of clusters.
#' @param ... Additional arguments if need to change distributions of the fixed terms
#' 
#' @return A list containing the following elements:
#' * mu_c: A numeric vector. The overall cluster mean. \in R^{C x 1}
#' * n_c: A numeric vector. The number of observations in each cluster. \in R^{C x 1}
#' * n: An integer. The total number of observations. \in R
#' * W_p: A numeric vector. The probability of being untreated. \in R^{C x 1}
#' * W_c: A numeric vector. The cluster level treatment assignment. W_c = 1 => never treated cluster \in R^{C x 1}
#' * indi_term: A numeric vector. The individual term. \in R^{n x 1}
#' * cluster: A numeric vector. The cluster. \in R^{n x 1}
#'
#' @examples
#' \dontrun{
#' fixed_value_simulation(5)
#' }
#' @export
fixed_value_simulation <- function(C, m, ...){
    res <- list()
    res$C <- C
    res$mu_c <- rnorm(C, 3, 1) #overall cluster mean is set to be 3, variance is 1
    res$n_c <- rep(m, C) # m observations in each cluster
    res$n <- sum(res$n_c) #total number of observations
    indi_term <- lapply(1:C, helper_indi_term_c, mu_c = res$mu_c, n_c = res$n_c)
    res$indi_term = unlist(indi_term)
    res$cluster = rep(1:C, res$n_c)

    # res$W_p <- 1 - gtools::inv.logit(res$mu_c, min = 0, max = 0.8)#probability of being untreated
    # res$W_c <- rbinom(C, 1, res$W_p) #cluster level treatment assignment, W_c = 1 means untreated
    # res$never_trt <- res$W_c[rep(1:C, res$n_c)]
    # res$trt <- 1 - res$never_trt
    return(res)
}

#2 simulate potential outcome under control Y_t(inf)
#'
#' This function simulates a dataframe of potential outcomes under control for a given number of time points such that Y_t(inf) = indi_term + v_t, where v_0 = norm(n, 0,1) and v_t = norm(n, v_t_minus_1, 1) for t = 1, 2, ..., t_total.
#'
#' @param t_total The total number of time points (excluding t = 0) for the simulation 
#' @param fixed_res A list returned from the 'fixed_value_simulation' function. 
#' The list 'res' should contain the following elements: "n", "indi_term".
#' @return A data frame where each column represents the Y_inf at a different time point \in R^{n x t_total + 1}]}
#' @export
sim_Y_t_inf <- function(t_total, fixed_res, beta=0){ #beta is the AR(1) coefficient, controls the correlation between time points
    t_seq <- 0:t_total
    v_0 <- rnorm(fixed_res$n, 0, 1)
    v_t <- Reduce(function(v_t_minus_1, t) {
        beta*v_t_minus_1 + rnorm(fixed_res$n,0,1)
    }, 1:t_total, init = v_0, accumulate = TRUE)
    df <- data.frame(do.call(cbind, v_t)) #create a df of time varying errors 
    df <- apply(df, 2, function(x) x + fixed_res$indi_term) #add individual term to each column
    colnames(df) <- t_seq
    return(df)
}


#generic, can be used for simulating cluster or individual level treatment assignment
# size = n or C depending on the level of treatment assignment
treat_sim_fun_uniform <- function(options, size){ #equal prob of being treated at any time point except for t = 0
    treat_prob <- rep(1/length(options), length(options))
    res <- data.frame(t(treat_prob))
    res <- res[rep(1, size),]
    colnames(res) <- options
    return(res)
    }
#3 simulate treatment distribution
# g \in {1, ..., t_total, t_total + 1} where t_total + 1 is the never/not yet treated 
#default that all units follow the same treatment distribution but can be changed via inputting a treat_prob_fun. 
# ... the output of fixed_value_simulation
#returns

#' Simulate treatment distribution
#'
#' This function simulates the distribution of treatment times for a set of units.
#' By default, all units have the same probability of being treated at any time point     #' from 1:t_total,
#' but a custom treatment probability function can be provided.
#'
#' @param t_total The total number of time points (excluding t = 0) for the simulation 
#' @param fixed_res A list containing fixed results from the 'fixed_value_simulation' function.
#' @param ... Additional arguments passed to the function.
#' @param treat_prob_fun A function that takes 't_total' as input and returns a vector of treatment probabilities.
#' If NULL (the default), all time points have equal probability of treatment.
#' @return A data frame with n rows and columns 'G' (the treatment time), 'cluster' (the cluster ID), 'Gmin' (the minimum treatment time in the cluster), and 'G_mod' (the modified treatment time, with 't_total + 1' indicating never treated).
#' @export
treatment_time_G <- function(t_total, fixed_res, df_prob_cluster_g, df_prob_indi_l, ...){
    #first simulate cluster level treatment time:
    gmin <- as.numeric(sapply(1:nrow(df_prob_cluster_g), function(x) sample(colnames(df_prob_cluster_g), size=1, prob=df_prob_cluster_g[x,])))
    #gmin <- treat_sim_fun(1:(t_total + 1), fixed_res$C, fixed_res$n) #number of options = t_total + 1 where t_total + 1 is the never/not yet treated
    gmin_n <- rep(gmin,fixed_res$n_c)
    #l <- treat_sim_fun(0:(t_total + 1), fixed_res$n, fixed_res$n) #simulate individual level treatment time l, where gmin_c + l_ci = g_ci 
    l <-  as.numeric(sapply(1:nrow(df_prob_indi_l), function(x) sample(colnames(df_prob_indi_l), size=1, prob=df_prob_indi_l[x,])))
    df <- data.frame(cluster=fixed_res$cluster, cluster_gmin=gmin_n, l = l, g = gmin_n + l)
    return(df)
}

#' Calculate the treatment effect over time
#'
#' This function calculates the treatment effect over time based on the proportion of treatment,
#' individual treatment, and past treatment effect. The treatment effect at each time point is
#' calculated as the sum of the cluster effect, individual effect, and past effect.
#'
#' @param t_total The total time period (excluding time = 0).
#' @param fixed_res A list containing fixed results from the 'fixed_value_simulation' function.
#' @param treatment_df A data frame containing the treatment info G, cluster, G_min, and G_mod.
#' @param cluster_coef The coefficient for the cluster effect.
#' @param indi_coef The coefficient for the individual effect.
#' @param past_coef The coefficient for the past effect.
#'
#' @return A data frame containing the treatment effect at each time point.
#' @export
treatment_effect <- function(t_total, fixed_res, treatment_df, cluster_coef, indi_coef, past_coef){
    prop_trt <- function(time) {
        return(ave(treatment_df$g <= time, treatment_df$cluster, FUN = mean))
    }
    indi_trt <- function(time) {
        return(treatment_df$g <= time)
    }
    TE_0 <- rep(0, fixed_res$n)
    TE_t <- Reduce(function(TE_t_minus_1, t) {
        cluster_coef*prop_trt(t) + indi_coef*indi_trt(t) + past_coef*TE_t_minus_1
    }, 1:t_total, init = TE_0, accumulate = TRUE)
    df <- data.frame(do.call(cbind, TE_t)) #create a df of of the treatment effect at each time point
    colnames(df) <- 0:t_total
    return(df)
}

#' Generate a data frame of potential outcomes over time
#'
#' This function generates a data frame of potential outcomes over time based on the treatment effect,
#' the simulated Y_t_inf, and the treatment data. The potential outcomes are calculated as the sum of
#' Y_t_inf and the treatment effect. The function can return the data in either long or wide format.
#'
#' @param t_total The total time period (excluding t = 0).
#' @param fixed_res A list containing fixed results from the 'fixed_value_simulation' function.
#' @param treatment_df A data frame containing the treatment info G, cluster, G_min, and G_mod.
#' @param cluster_coef The coefficient for the cluster effect.
#' @param indi_coef The coefficient for the individual effect.
#' @param past_coef The coefficient for the past effect.
#' @param long Logical indicating whether to return the data in long format. If TRUE, the data is returned
#' in long format with one row per time point per individual. If FALSE, the data is returned in wide format
#' with one row per individual and one column per time point. Default is TRUE.
#'
#' @return A data frame of potential outcomes over time. If `long` is TRUE, the data frame also includes
#' the treatment effect and Y_t_inf at each time point.
#' @export
#'
#' @examples
#' t_total <- 10
#' C = 10 
#' fixed_res <- fixed_value_simulation(C)
#' treatment_df <- treatment_time_G(t_total, fixed_res)
#' cluster_coef <- 3
#' indi_coef <- 3
#' past_coef <- 0.5
#' potential_outcome_df(t_total, fixed_res, treatment_df, cluster_coef, indi_coef, past_coef)
potential_outcome_df <- function(t_total, fixed_res, treatment_df, cluster_coef, indi_coef, past_coef,long=TRUE){
    df_TE <- treatment_effect(t_total, fixed_res, treatment_df, cluster_coef, indi_coef, past_coef)
    df_Y_inf <- sim_Y_t_inf(t_total, fixed_res)
    df_Y_G <- df_Y_inf + df_TE
    df_Y_G$id <- 1:fixed_res$n
    df_Y_G$cluster <- fixed_res$cluster
    df_Y_G$g <- treatment_df$g
    df_Y_G$cluster_gmin <- treatment_df$cluster_gmin

    if (long){
        df_Y_G_long <- reshape2::melt(df_Y_G, id.vars = c("id", "cluster", "g", "cluster_gmin"), variable.name = "time", value.name = "Y_G")
        df_Y_G_long$TE <- reshape2::melt(df_TE, id.vars = NULL, variable.name = "time", value.name = "TE")$TE
        df_Y_G_long$Y_inf <- reshape2::melt(df_Y_inf, id.vars = NULL, variable.name = "time", value.name = "Y_inf")$Y_inf
        df_Y_G_long$time <- as.numeric(df_Y_G_long$time) - 1 #time starts from 0, turn factors into numbers that starts at 1
        return(df_Y_G_long)
    } else {
        return(df_Y_G)
    }
}

observed_df <- function(df_long) {
    res <- df_long[, c("time", "id", "cluster", "g", "cluster_gmin", "Y_G")]
    colnames(res) <- c("time", "id", "cluster", "g", "cluster_gmin", "Y_obs")
    res$indi_trt <- res$g <= res$time #Z_t
    return(res)
}

experiment <- function(t_total, C, m, g, k, cluster_coef, indi_coef, past_coef){
    fixed_res <- fixed_value_simulation(C, m)
    df_prob_cluster_g <- treat_sim_fun_uniform(1:(t_total + 1),fixed_res$C)
    df_prob_cluster_g <- as.data.frame(cbind(rep(0, nrow(df_prob_cluster_g)), df_prob_cluster_g)) #adding the column for time = 0, since no cluster is treated at time = 0
    colnames(df_prob_cluster_g) <- 0:(t_total + 1) # rename the columns.
    df_prob_indi_l <- treat_sim_fun_uniform(0:t_total,fixed_res$n) # l \in {0,1,...,t_total}
    treatment_df <- treatment_time_G(t_total, fixed_res, df_prob_cluster_g, df_prob_indi_l)

    df_long <- potential_outcome_df(t_total, fixed_res, treatment_df, cluster_coef, indi_coef, past_coef,long=TRUE)
    df_wide <- potential_outcome_df(t_total, fixed_res, treatment_df, cluster_coef, indi_coef, past_coef, long=FALSE)
    df_obs <- observed_df(df_long)

    theta_tilde_kg <- theta_tilde_k(df_long, fixed_res, k, df_prob_cluster_g, t_total)$theta_tilde_k_g_s[g]
    #theta_tilde_kg <- sum(theta_tilde_kg_df$w * theta_tilde_kg_df$TE)

    cluster_min <- aggregate(cluster_gmin ~ cluster, treatment_df, mean)$cluster_gmin
    gmin_frequency <- table(factor(cluster_min, levels = 1:(t_total + 1)))
    #df_prob_cluster_g_est <- df_prob_cluster_g_est/sum(df_prob_cluster_g_est) #empirical probability of each gmin
    theta_hat_kg <- theta_hat_k(df_obs, fixed_res, k,g, gmin_frequency, t_total)$theta_hat_k_g_s[g]
    error <- theta_hat_kg - theta_tilde_kg
    return(list(theta_tilde_kg = theta_tilde_kg, theta_hat_kg = theta_hat_kg, error = error))
}

experiment_k <- function(t_total, C, m, g, k, cluster_coef, indi_coef, past_coef){
    fixed_res <- fixed_value_simulation(C, m)
    df_prob_cluster_g <- treat_sim_fun_uniform(1:(t_total + 1),fixed_res$C)
    df_prob_cluster_g <- as.data.frame(cbind(rep(0, nrow(df_prob_cluster_g)), df_prob_cluster_g)) #adding the column for time = 0, since no cluster is treated at time = 0
    colnames(df_prob_cluster_g) <- 0:(t_total + 1) # rename the columns.
    df_prob_indi_l <- treat_sim_fun_uniform(0:t_total,fixed_res$n) # l \in {0,1,...,t_total}
    treatment_df <- treatment_time_G(t_total, fixed_res, df_prob_cluster_g, df_prob_indi_l)

    df_long <- potential_outcome_df(t_total, fixed_res, treatment_df, cluster_coef, indi_coef, past_coef,long=TRUE)
    df_wide <- potential_outcome_df(t_total, fixed_res, treatment_df, cluster_coef, indi_coef, past_coef, long=FALSE)
    df_obs <- observed_df(df_long)

    theta_tilde_k_res <- theta_tilde_k(df_long, fixed_res, k, df_prob_cluster_g, t_total)
    theta_tilde_k <- theta_tilde_k_res$theta_tilde_k
    #theta_tilde_kg <- sum(theta_tilde_kg_df$w * theta_tilde_kg_df$TE)

    cluster_min <- aggregate(cluster_gmin ~ cluster, treatment_df, mean)$cluster_gmin
    gmin_frequency <- table(factor(cluster_min, levels = 1:(t_total + 1)))
    #df_prob_cluster_g_est <- df_prob_cluster_g_est/sum(df_prob_cluster_g_est) #empirical probability of each gmin
    theta_hat_k_res <- theta_hat_k(df_obs, fixed_res, k,g, gmin_frequency, t_total)
    theta_hat_k <- theta_hat_k_res$theta_hat_k
    error <- theta_hat_k - theta_tilde_k
    omega_hat_n <- omega_hat_n(theta_hat_k_res$theta_hat_k_ci, fixed_res)
    omega_tilde_n <- omega_tilde_n(theta_hat_k_res$theta_hat_k_ci, theta_tilde_k_res$theta_tilde_k_ci, fixed_res)
    return(list(theta_tilde_k = theta_tilde_k, theta_hat_k = theta_hat_k, error = error, omega_hat= omega_hat_n))
}

#' @export
data_sim_single <-  function(t_total, C, m, cluster_coef, indi_coef, past_coef) {
    fixed_res <- fixed_value_simulation(C, m)
    df_prob_cluster_g <- treat_sim_fun_uniform(1:(t_total + 1),fixed_res$C)
    df_prob_cluster_g <- as.data.frame(cbind(rep(0, nrow(df_prob_cluster_g)), df_prob_cluster_g)) #adding the column for time = 0, since no cluster is treated at time = 0
    colnames(df_prob_cluster_g) <- 0:(t_total + 1) # rename the columns.
    df_prob_indi_l <- treat_sim_fun_uniform(0:t_total,fixed_res$n) # l \in {0,1,...,t_total}
    treatment_df <- treatment_time_G(t_total, fixed_res, df_prob_cluster_g, df_prob_indi_l)

    df_long <- potential_outcome_df(t_total, fixed_res, treatment_df, cluster_coef, indi_coef, past_coef,long=TRUE)
    df_wide <- potential_outcome_df(t_total, fixed_res, treatment_df, cluster_coef, indi_coef, past_coef, long=FALSE)
    df_obs <- observed_df(df_long)
    cluster_min <- aggregate(cluster_gmin ~ cluster, treatment_df, mean)$cluster_gmin
    gmin_frequency <- table(factor(cluster_min, levels = 1:(t_total + 1)))

    return(list(df_long = df_long, df_wide = df_wide, df_obs = df_obs, fixed_res = fixed_res, treatment_df = treatment_df, df_prob_cluster_g = df_prob_cluster_g, df_prob_indi_l = df_prob_indi_l, t_total = t_total, gmin_frequency = gmin_frequency))

}

#' @export
single_run <- function(sim_num, generate_sim, get_results, generate_sim_args = list(), get_results_args = list()) {

  cat("." )
  if(!(sim_num %% 100)) {
      cat("\n Simulation number", sim_num, ": ")
  }
  
  # simulate data
  output <- do.call(generate_sim, generate_sim_args)

  # get results
#   get_results_quietly <- purrr:::quietly(get_results)
#   results <-  get_results_quietly(output, ...)$result
  get_results_quietly <- purrr:::quietly(get_results)
  args <- c(output, get_results_args)
  matched_args <- intersect(names(args), names(formals(get_results)))
  results <-  unlist(do.call(get_results_quietly, args[matched_args])$result)

  #results$sim_num <- sim_num
  # return the results
  return(results)
}

#' Run many simulations in parallel
#' @param n_sims Number of simulations to run
#' @param n_cores Number of cores
#' @param ... Simulation arguments
# generate_sim, get_results, generate_sim_args = list(), get_results_args = list(),
#' @export
run_sim <- function(n_sims, n_cores=1,  ...) {


    cat("Number of simulations", n_sims, "\n")
    ## run n_sim simulations and evaluate
    out <- parallel::mclapply(
        1:n_sims,
        function(j) single_run(j, ...),
        mc.cores=n_cores)

    #out <- out[lapply(out, function(x) !is.null(names(x))) == TRUE]
    cat("\n\n")
    #return(dplyr::bind_rows(out))
    return(out)
}







