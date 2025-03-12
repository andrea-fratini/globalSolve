################# Error handling ##################

rm(list=ls())

# Set of functions to obtain the hierarchical structure

dynamic_to_static <- function(variables) {
  gsub("_[fb][0-9]+", "", variables)  # Remove time subscripts like "_f1" or "_b1"
}

# Function to count the number of endogenous and state_endogenous variables

count_vars <- function(foc, list_of_vars) {
  variables <- all.vars(foc)  # Extract all variables
  print(variables) ####
  static_variables <- dynamic_to_static(variables)  # Normalize variable name
  print(static_variables) ####
  # Count how many belong to endogenous and state_endogenous
  sum(static_variables %in% list_of_vars$endogenous)
}

hierarchical_structure <- function(FOCs, list_of_vars){
  
  counts <- sapply(FOCs, count_vars, list_of_vars = list_of_vars)
  
  # Order the equations based on counts
  FOCs_ordered <- FOCs[order(counts)]
  
  # Print ordered equations
  return(list(FOCs=FOCs_ordered, n_endo=counts[order(counts)]))
}

# Function to obtain the general timing of the variables

check_and_timing <- function(FOCs, list_of_vars, params) {
  #' @param FOCs set of First Order Conditions to be checked
  #' @param list_of_vars set of variables defined by the user
  #' @param params set of parameters defined by the user
  #' @description
  #' The function takes as argument the set of First Order Conditions
  #' and checks the correctness of the specification, 
  #' in particular:
  #' Checks if the temporal dimension is correctly specified,
  #' Checks if the variables has been defined by the user
  #' Checks the timing of the variable for each one of the variables and returns it
  #' WHY: The goal is to use this functions to get the temporal dimension of the variables in the FOC's
  #' NOTICE: this function and hierarchical_structure share some operations that could be optimized in just one
  
  all_variables <- unlist(list_of_vars)
  all_parameters <- unlist(parameters)
  names(all_parameters) <- gsub("_\\.", "_", names(all_parameters))
  
  variables <- unique(unlist(lapply(FOCs, function(foc) setdiff(all.vars(foc), names(all_parameters)))))
  
  
  # Function to validate the format of each variable
  validate_var_format <- function(var) {
    # The regular expression that matches the form: "x", "x_ft", "x_bt"
    return(grepl("^([a-z]+)(?:_([fb][0-9]+))?$", var))
  }
  
  parsed <- lapply(variables, function(var) {
    # Check if the variable matches the expected format
    if (!validate_var_format(var)) {
      stop(paste0("The variable \"", var, "\" does not match the expected format", call. = FALSE))
    }
    
    # Extract the base variable and time component
    match <- regmatches(var, regexec("([a-z]+)(?:_([fb][0-9]+))?", var))[[1]]
    
    if (!match[2] %in% all_variables) {
      stop(paste0("The variable \"", var, "\" is not in the list of variables", call. = FALSE))
    }
    
    base_var <- match[2]
    
    # Handle the time component
    if (match[3] != "") {
      time <- as.numeric(sub("f", "", sub("b", "-", match[3]))) # Convert "t2" to 2, "m2" to -2
    } else {
      time <- 0 # Assign 0 for base case
    }
    
    return(list(name = base_var, time = time))
  })
  
  # Group by base variable name
  grouped <- split(parsed, sapply(parsed, function(x) x$name))
  
  # Convert to desired format
  result <- lapply(grouped, function(lst) {
    times <- sapply(lst, function(x) x$time)
    return(times)
  })
  
  return(result)
}




########### example ############

technical_params <- list(n_grid=10 # it might be a vetor of different dimensions
                         )

variables <- list(endogenous=c("cc"), state_endogenous=c("b"), 
                  state_exogenous=c("z"), shocks=c("epsz"))

parameters <- list(bbeta=0.99, ggamma=2, r=1.025, 
                   rho_z=0.6, sigma_z=0.02 #
                   )

grids <- list(b=c(-0.7, 0))

FOCs <- list(f_Euler_cons = cc^(-ggamma) ~ bbeta*r*cc_f1^(-ggamma),
             f_Endo_b = b_f1 ~ exp(z) + b - cc,
             exo_z = z ~ rho_z*z_b1 + epsz)

Ordered_focs <- hierarchical_structure(FOCs, variables)

Timing <- check_and_timing(Ordered_focs$FOCs, variables, parameters)

shocks_grids_and_probabilities <- function(variables, parameters, technical_params){
  exo <- names(variables$state_exogenous)
  Grids <- matrix(NA, nrow = technical_params$n_grid, # it might be a vetor of different dimensions
                  ncol = length(exo))
  Probs <- array(NA, dim = c(technical_params$n_grid, technical_params$n_grid, # it might be a vetor of different dimensions
                             length(exo)))
  for(i in 1:length(exo)){
    Rtauchen::Rtauchen(technical_params$n_grid,# it might be a vetor of different dimensions
                       parameters,
                       parameters,
                       parameters)
    Rtauchen::Tgrid(technical_params$n_grid,# it might be a vetor of different dimensions
                    parameters,
                    parameters,
                    parameters)
  }
} ### modify

System <- function()





# 
# 
# check_and_timing(FOCs, variables)
# 
# variables <- unique(unlist(...))
# 
# varia <- lapply(FOCs, function(foc) setdiff(all.vars(foc), names(parameters)))
# 
# 
# 
# n_state_endo <- lapply(varia, function(v){sum(v %in% c(variables$state_endogenous))})
# 
# n_exo-n_state_endo
