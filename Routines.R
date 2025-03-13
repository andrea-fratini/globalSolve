############### Basis Functions #######

#### chebyshev 

library(pracma)

chebyshev_nodes <- function(n){cos(pi * (2*(1:n)- 1)/(2*n))}

# Function to shrink grid values to Chebyshev space (-1,1)
shrink_grid_chebyshev <- function(grid) {
  return(2 * (grid - min(grid)) / (max(grid) - min(grid)) - 1)
}

# Function to expand back from Chebyshev space to original scale
expand_grid_chebyshev <- function(grid) {
  return((1 + grid) * (max(grid) - min(grid)) / 2 + min(grid))
}

# Function to compute Chebyshev polynomials up to a given degree
chebyshev <- function(grid, degree) {
  x <- grid
  ret_cheb <- matrix(NA, nrow = length(x), ncol = degree)
  
  for (i in 0:(degree - 1)) {
    ret_cheb[, i + 1] <- chebPoly(n = i, x = as.matrix(x))  # Compute Chebyshev polynomial
  }
  
  return(ret_cheb)
}

# Compute the tensor product of Chebyshev polynomials across multiple dimensions
tensor_product_chebyshev <- function(grids, degree) {
  
  # Compute Chebyshev polynomials for each grid separately
  cheb_list <- lapply(grids, function(grid) chebyshev(grid, degree))
  
  # Take the Kronecker product across all dimensions
  tensor_cheb <- Reduce(kronecker, cheb_list)
  
  return(tensor_cheb)
}

wrapper_chebyshev <- function(technical_params, variables){
  
  grids_specification <- technical_params$state_endo_grids
  
  if(length(grids_specification) != length(variables$state_endogenous)){
    
    stop(paste("The number of state endgenous", length(variables$state_endogenous), 
               "does not match the number of grids specifications", length(grids_specification)))
      
  }
  
  nodes <- lapply(grids_specification, function(spec){chebyshev_nodes(spec[3])})
  
  return(tensor_product_chebyshev(nodes, technical_params$poly_type$degree))
  
}


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

#Extract timing for classes of variables

extract_timing <- function(timing, variables){
  return(timing[which(names(timing) %in% variables$state_exogenous)])
}

# Function to generate grids and transition matrices

shocks_grids_and_probabilities <- function(variables, parameters, technical_params){
  
  exo <- variables$state_exogenous
  Grids <- list()
  Probs <- list()
  
  for(i in 1:length(exo)){
    Probs[[exo[i]]] <- Rtauchen::Rtauchen(technical_params$n_grid[i],# it might be a vetor of different dimensions
                                          parameters$sigma_[i],
                                          parameters$rho_[i],
                                          technical_params$width_grid[i])
    
    Grids[[exo[i]]] <- Rtauchen::Tgrid(technical_params$n_grid[i],# it might be a vetor of different dimensions
                                       parameters$sigma_[i],
                                       parameters$rho_[i],
                                       technical_params$width_grid[i]) + parameters$mu_[i]
  }
  
  return(list(Grids=Grids, Probs=Probs))
  
}


########### example ############

technical_params <- list(n_grid=c(10,5), width_grid=c(3,3),
                         state_endo_grids=list(b=c(-0.7, 0, 10)),
                         poly_type=list(type="chebyshev", degree=3))

variables <- list(endogenous=c("cc"), state_endogenous=c("b"), 
                  state_exogenous=c("z", "r"), shocks=c("epsz", "epsr"))

parameters <- list(bbeta=0.99, ggamma=2, mu_=c(z=0,r=1.025), 
                   rho_=c(z=0.6, r=0.8), sigma_=c(z=0.02, r=0.03))

grids <- list(b=c(-0.7, 0))

FOCs <- list(f_Euler_cons = cc^(-ggamma) ~ bbeta*r*cc_f1^(-ggamma),
             f_Endo_b = b_f1 ~ exp(z) + b - cc,
             exo_z = z ~ rho_z*z_b1 + epsz,
             exo_r = r ~ rho_r*r_b1 + epsr)

Ordered_focs <- hierarchical_structure(FOCs, variables)

Timing <- check_and_timing(Ordered_focs$FOCs, variables, parameters)

 ### modify

shocks_grids <- shocks_grids_and_probabilities(variables, parameters, technical_params)


System <- function(FOCs, variables, parameters, technical_params, Ordered_focs, Timing, shocks_grids){
  
  # expand grids for exogenous processes
  
  exo_grid <- expand.grid(shocks_grids$Grids)
  
  exo_grid_idx <- expand.grid(lapply(technical_params$n_grid, seq_len)); names(exo_grid_idx) <- names(exo_grid)
  
  n_exo = dim(exo_grid)[1]

  state_exo_grid <- expand.grid(technical_params$state_endo_grids)
  
  state_exo_grid_idx <- expand.grid(lapply(lapply(technical_params$state_endo_grids, length), seq_len)); names(exo_grid_idx) <- names(exo_grid)
  
  n_state_exo <- dim(exo_grid)[1]
  
  for(exo in 1:n_exo){
    
    exo_b1 <- exo_grid[exo]
    exo_idx_b1 <- exo_grid_idx[exo]
    
    for(state_exo in 1:n_state_exo){
      
      state_exo_b1 <- state_exo_grid[state_exo]
      state_exo_idx_b1 <- state_exo_grid_idx[state_exo]
      
      BASIS
      
    }

  }
  
}




