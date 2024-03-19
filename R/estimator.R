# Description: This file contains the functions to estimate the parameters of the model

# returns the weights on each ci, k periods after treatment. if observed is true, then df = df_observed containing variable: "time", "id", "cluster", "g", "cluster_min", "Y_obs", it returns the DiD weights for the observed data where Gmin_ci < T are in the control group. If observed is false, then df = df_long, the true data, containing variable: "time", "id", "cluster", "g", "cluster_min", "Y_G", TE, "P_Gmin" it returns the weights on the TE_ci of the random version of the true estimator theta_k_weights. 
# theta_k_weights <- function(df, k, observed=TRUE, ...) {

# }