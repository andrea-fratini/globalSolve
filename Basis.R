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
               "does not match the number of grids specifications", length(grids_specification)), call. = F)
    
  }
  
  nodes <- lapply(grids_specification, function(spec){chebyshev_nodes(spec[3])})
  
  return(tensor_product_chebyshev(nodes, technical_params$basis_type$degree))
  
}

## wrapper basis

wrapper_basis <- function(technical_params, variables){
  
  if(! technical_params$basis_type$type %in% c("Chebyshev")){
    
    stop(paste(technical_params$basis_type$type, 
               "not in the list of supported basis functions"), call. = F)
    
  }
  
  if(technical_params$basis_type$type == "Chebyshev"){
    
    Basis <- wrapper_chebyshev(technical_params, variables)
    
  }
  
  return(Basis)
  
}


