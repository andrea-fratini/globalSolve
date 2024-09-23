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


