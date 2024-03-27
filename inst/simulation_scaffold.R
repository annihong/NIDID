
# single_run <- function(sim_num, generate_sim, get_results, generate_sim_args = list(), get_results_args = list()) {

#   cat("." )
#   if(!(sim_num %% 100)) {
#       cat("\n Simulation number", sim_num, ": ")
#   }
  
#   # simulate data
#   output <- do.call(generate_sim, generate_sim_args)

#   # get results
# #   get_results_quietly <- purrr:::quietly(get_results)
# #   results <-  get_results_quietly(output, ...)$result
#   get_results_quietly <- purrr:::quietly(get_results)
#   args <- c(output, get_results_args)
#   matched_args <- intersect(names(args), names(formals(get_results)))
#   results <-  unlist(do.call(get_results_quietly, args[matched_args])$result)

#   #results$sim_num <- sim_num
#   # return the results
#   return(results)
# }

# #' Run many simulations in parallel
# #' @param n_sims Number of simulations to run
# #' @param n_cores Number of cores
# #' @param ... Simulation arguments
# # generate_sim, get_results, generate_sim_args = list(), get_results_args = list(),
# run_sim <- function(n_sims, n_cores=1,  ...) {


#     cat("Number of simulations", n_sims, "\n")
#     ## run n_sim simulations and evaluate
#     out <- parallel::mclapply(
#         1:n_sims,
#         function(j) single_run(j, ...),
#         mc.cores=n_cores)

#     #out <- out[lapply(out, function(x) !is.null(names(x))) == TRUE]
#     cat("\n\n")
#     #return(dplyr::bind_rows(out))
#     return(out)
# }



library(NIDID)


if(!interactive()) {


  # define arguments here, e.g.
  n_sims <- 500
  n_cores <- 4



  # args <- commandArgs(TRUE)
  # if(length(args) != 0) {
  #     for(i in 1:length(args)) {
  #         eval(parse(text = args[i]))
  #     }
  # }
  # print(args)

  # put in other function args here
  res <- run_sim(n_sims, n_cores, generate_sim = data_sim_single, get_results = theta_tilde_k, generate_sim_args = list(t_total=10, C=50, m=25, cluster_coef=3, indi_coef=3, past_coef=1), get_results_args = list(k = 1, result = "theta_tilde_k_ci"))
  
  # save file
  timestamp <- format(Sys.time(), "%d-%b-%Y-%H-%M-%S")
  
  
  sim_name <- "tidle_k_ci_clusterw"

  path <- ""

  res_name <- paste(path,
                    paste(sim_name,
                          timestamp, sep = "_"),
                    ".rds", sep = "")
  save(res, file=res_name)
}

