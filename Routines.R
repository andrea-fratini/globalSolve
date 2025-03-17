rm(list=ls())

############### Basis Functions #######

#### chebyshev 

library(pracma)

chebyshev_nodes <- function(n){cos(pi * (2*(1:n)- 1)/(2*n))}

# Function to shrink grid values to Chebyshev space (-1,1)
shrink_grid_chebyshev <- function(grid) {
  return(2 * (grid - min(grid)) / (max(grid) - min(grid)) - 1)
}

# Function to expand back from Chebyshev space to original scale
expand_grid_chebyshev <- function(grid, grid_spec) {
  return((1 + grid) * (grid_spec[2] - grid_spec[1]) / 2 + grid_spec[1])
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
               "does not match the number of grids specifications", length(grids_specification)), call. = F)
    
  }
  
  nodes <- lapply(grids_specification, function(spec){expand_grid_chebyshev(chebyshev_nodes(spec[3]), spec)})
  
  
  return(list(basis=tensor_product_chebyshev(nodes, technical_params$basis_type$degree), nodes=nodes))
  
}

## wrapper basis

wrapper_basis <- function(technical_params, variables){
  
  if(! technical_params$basis_type$type %in% c("Chebyshev")){
    
    stop(paste(technical_params$basis_type$type, 
               "not in the list of supported basis functions"), call. = F)
    
  }
  
  if(technical_params$basis_type$type == "Chebyshev"){
    
    Basis_nodes <- wrapper_chebyshev(technical_params, variables)
    
    
  }
  
  return(Basis_nodes)
  
}

################# Error handling ##################




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

check_and_timing <- function(FOCs, list_of_vars, params, single=F) {
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
  
  if(single==T){
    
    variables <- setdiff(all.vars(unlist(FOCs)), names(parameters))
    
  }
  
  
  # Function to validate the format of each variable
  validate_var_format <- function(var) {
    # The regular expression that matches the form: "x", "x_ft", "x_bt"
    return(grepl("^([a-z]+)(?:_([fb][0-9]+))?$", var))
  }
  
  parsed <- lapply(variables, function(var) {
    # Check if the variable matches the expected format
    if (!validate_var_format(var)) {
      stop(paste0("The variable \"", var, "\" does not match the expected format"), call. = FALSE)
    }
    
    # Extract the base variable and time component
    match <- regmatches(var, regexec("([a-z]+)(?:_([fb][0-9]+))?", var))[[1]]
    
    if (!match[2] %in% all_variables) {
      stop(paste0("The variable \"", var, "\" is not in the list of variables"), call. = FALSE)
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

f2f = function(input, param){
  #' @param input input of the formula object
  #' @param param values of the parameters associated to the formula
  
  if(is_formula(input) == F | is.list(param) == F){ # small check
    print("Invalid Input")
    return(NULL)
    
  }
  
  dep = f_lhs(input)
  indep = f_rhs(input)
  param_env <- new.env() # create a new environment
  list2env(c(param, list(dep = dep, indep = indep)), envir = param_env)  # Insert everything in the new environment
  
  val = eval(indep, envir = param_env)
  return(val)
}


########### example ############

technical_params <- list(n_grid=c(10), width_grid=c(3), # Add Exo
                         state_endo_grids=list(b=c(-0.7, 0, 10), r=c(0.9, 1, 10)),
                         basis_type=list(type="Chebyshev", degree=3))

variables <- list(endogenous=c("cc"), state_endogenous=c("b", "r"), 
                  state_exogenous=c("z"), shocks=c("epsz", "epsr"))

parameters <- list(bbeta=0.99, ggamma=2, mu_=c(z=0), 
                   rho_=c(z=0.6), sigma_=c(z=0.02))

FOCs <- list(f_Endo_Euler_cons = cc^(-ggamma) ~ bbeta*r*cc_f1^(-ggamma),
             f_State_Endo_b = b_f1 ~ exp(z) + b - cc,
             f_State_Endo_r = r_f1 ~ exp(z) + r - cc,
             exo_z = z_f1 ~ rho_z*z + epsz_f1)


Ordered_focs <- hierarchical_structure(FOCs, variables)

Timing <- check_and_timing(Ordered_focs$FOCs, variables, parameters)
 ### modify

shocks_grids <- shocks_grids_and_probabilities(variables, parameters, technical_params)

Basis <- wrapper_basis(technical_params, variables)

time_span <- unique(unlist(Timing))

coeffs <- rep(1, n_endo*technical_params$basis_type$degree^(n_endo)*n_exo)

# solve for state endo!!!!
# function to get the variables needed to be evaluated

System <- function(FOCs, variables, parameters, technical_params, Ordered_focs, Timing, shocks_grids, basis, coeffs, time_span){
  
  # expand grids for exogenous processes
  
  exo_grid <- expand.grid(shocks_grids$Grids)
  
  exo_grid_idx <- expand.grid(lapply(technical_params$n_grid, seq_len)); names(exo_grid_idx) <- names(exo_grid)
  
  n_exo = dim(exo_grid)[1]

  state_exo_grid <- expand.grid(Basis$nodes) # expand depending on the type of basis, to be corrected!!!
  
  state_exo_grid_idx <- expand.grid(lapply(lapply(Basis$nodes, length), seq_len)); names(exo_grid_idx) <- names(exo_grid)
  
  n_state_exo <- dim(state_exo_grid)[1]
  
  n_endo <- length(variables$endogenous)
  
  # coeffs structure
  
  coeffs <- array(coeffs, dim=c(n_endo, technical_params$basis_type$degree^(n_endo), n_exo)) # da capire le dimensioni !!!!
  
  # Divide the variables in classes
  
  Endo <- FOCs[grep("^f_Endo_", names(FOCs))]
  State_Endo <- FOCs[grep("^f_State_Endo_", names(FOCs))]
  Exo <- FOCs[grep("^exo_", names(FOCs))]
  
  # get the time span of the variables
  
  times <- unique(unlist(Timing))
  times <-  sort(times)
  
  for(t in 1:(length(times)-1)){
    
    to_update <- names(Timing)[sapply(Timing, function(x) any(x == times[t+1]))]
    
    for(exo in 1:n_exo){
      
      exo_b1 <- exo_grid[exo,]
      exo_idx_b1 <- exo_grid_idx[exo,]
      
      temp_coeffs <- coeffs[,,n_exo]
      
      state_exo_ <- state_exo_grid
      
      for(state_exo in 1:n_state_exo){
        
        state_exo_b1 <- state_exo_grid[state_exo,]
        state_exo_idx_b1 <- state_exo_grid_idx[state_exo,]
        
        endo_ <- coeffs[,,state_exo] %*% Basis$basis[state_exo,] # assicurarsi che l'ordine del tensor product sia lo stesso della griglia
        # si potrebbe ottimizzare calcolando il prodotto solo per le endo presenti in t=0
        
        endo_ <- as.numeric(endo_)
        
        if(times[t] != 0){
          names(endo_) <- paste(variables$endogenous, "_", ifelse(times[t] < 0, "b", "f"), times[t], sep = "") # cambia il nome se è != 0
        }
        
        vars_ <- unlist(c(state_exo_b1, endo_))
        
        # calcolo le griglie per il prossimo periodo

        names(Timing)[sapply(Timing, function(x) any(x == times[t+1]))]

        check_and_timing(State_Endo$f_State_Endo_b, variables, parameters, single = T) # Forse non serve, basta dargli le variabili di cui ha bisogno!
        # Cosa fare se nel modello la state endo si updata con variabili avanti nel futuro o molto indietro nel passato? ho bisogno di variabili a un passo
        # per aggiornare le griglie!
        
        States <- lapply(State_Endo, FUN = check_and_timing, list_of_vars=variables, params=parameters, single=T)
        
        # calcolare il valore e generare la nuova griglia
        
        # trasformare le variabili al tempo con il loro nome preciso
        
        
        
      }
      
    }
    
  }
  
}




