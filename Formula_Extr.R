# Notation
#
# x = x_{t}, x_m1 = x_{t-1}, x_1 = x_{t+1}
#

#install.packages("rlang")
library(rlang);library(dplyr)

##### Let's try to build a function that preprocess a formula object in R #####

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




param <- list(a = 1, d = 2, x = 2, z = 3)
input = y ~ a * x + d * z

f2f(input, param)


##### Let's move on and try to understand how to work with real equations ####

parameters <- list(bbeta=0.99, ggamma=2, r=1.025, rho_z=0.6, sigma_z=0.02)

grids <- list(b=c(-0.7, 0)) # defines the support of b 

f_Euler_cons <- cc^(-ggamma) ~ bbeta*r*cc_1^(-ggamma)

f_Endo_b <- b_1 ~ exp(z) + b - cc

exo_z <- z ~ rho_z*z_m1 + epsilon_z


compose_function = function(..., param){
 #' @param ... set of equations defined by the user
 #' @param param list of parameters defined
 #' @description
  #' The idea of this function is pretty straightforward. Given a set of functions
  #' given by the user following a certain order, it automatically find the formulas
  #' and the associated parameters to each function. In this way, it is also possible
  #' to update easily the parameters.
 
 stor = list()
 if(!is.list(param)){ # start with a check
  print("Invalid Type od parameters")
  return(NULL)
 }
 
 # Take back the original names of the function given by the user
 fun_names = sapply(substitute(...()), deparse)
 ass_for = list(...)


 for (i in 1:length(fun_names)) {
  param_to_add = all.vars(f_rhs(ass_for[[i]]))[(all.vars(f_rhs(ass_for[[i]])) %in% names(param)) == F] 
  param = append(param, rep(1, length(param_to_add)))
  names(param)[names(param) == ""] = param_to_add
  actual_param = param[names(param) %in% all.vars(f_rhs(ass_for[[i]]))] # forse da rimuovere
  stor[[fun_names[i]]] <- list(input = ass_for[[i]], params = actual_param )
 }
 return(stor)
 
 
}

ass_for = compose_function(f_Endo_b, f_Euler_cons, exo_z,
                             param = parameters)
