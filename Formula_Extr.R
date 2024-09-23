#install.packages("rlang")
library(rlang)

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


##### Let's move on and try to undestand how to work with real equations ####

f_Endo_b = bp ~ r*(y+b-cc)

f_Euler_cons = cc^(-ggamma) ~ bbeta*r*ccp^(-ggamma) + mu

f_Euler_price = q ~ bbeta*(ccp^(-ggamma)*(qp + aalpha*yp))/(cc^(-ggamma)-phi*mu)

f_Val_fun  = v ~ cc^(1-ggamma)/(1-ggamma) + bbeta*vp


compose_function = function(..., param){
 #' @param ... set of equations defined by the user
 #' @param param list of parameters defined
 
 stor = list()
 if(!is.list(param)){ # start with a check
  print("Invalid Type od parameters")
  return(NULL)
 }
 
 # Take back the original names of the function given by the user
 fun_names = sapply(substitute(...()), deparse)
 ass_for = list(...)

 for (i in 1:length(fun_names)) {

  actual_param = param[names(param) %in% all.vars(f_rhs(ass_for[[i]]))]
  stor[[fun_names[i]]] <- list(input = ass_for[[i]], params = actual_param )
 }
 return(stor)
 
 
}


fun_names = compose_function(f_Endo_b, f_Euler_cons, f_Euler_price, f_Val_fun,
                             param = list(r = 1.025, bbeta = 0.99, ggamma = 2, 
                                          aalpha = 0.133, phi = 0.2))



