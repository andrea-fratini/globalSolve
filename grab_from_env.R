rm(list=ls())
set.seed(123)


#### Idee ed appunti #####:

# Visto che alcuni termini si ripetono all'interno delle formule
# renderebbe il codice piu efficiente calcolare questi termini comuni solo una volta
# e passarli come argomenti delle funzioni al posto delle variabili?
# Se si modificare il codice per estrarli::::per ora estrae solo PARAMETRI E VARIABILI
rm(list=ls())

library(Rtauchen);library(pracma);library(mpoly);library(MASS);library(rootSolve);library(BB);library(fBasics);library(plotly);library(splines)


f_Endo_b <- bp ~ r*(y+b-cc)

f_Euler_cons <- cc^(-ggamma) ~ bbeta*r*ccp^(-ggamma) + mu

f_Euler_price <- q ~ bbeta*(ccp^(-ggamma)*(qp + aalpha*yp))/(cc^(-ggamma)-phi*mu)

f_Val_fun <- v ~ cc^(1-ggamma)/(1-ggamma) + bbeta*vp

r <- 1.025
bbeta <- 0.99
ggamma <- 2
aalpha <- 0.133
phi <- 0.2

split_vars <- ls()

result <- lapply(split_vars, function(i) {
  var <- strsplit(i, split = "_")[[1]]
  if (var[1] == "f") {
    return(i)  # Return the original variable name if it starts with "f"
  }
})

result

result <- Filter(Negate(is.null), result)

result <- unlist(result)

varz <- ls()
eval(parse(text = result))

param <- which(varz %in% all.vars(eval(parse(text = result[1])))[-1])

eval(parse(text=varz[param]))

bp <- formula2func(eval(parse(text = result[1])))
?parse
eval(parse(text = paste("bp(", varz[param], "=",eval(parse(text=varz[param])), ")")))

# next_b <- function(y,B,cc,params){
#   R <- params$R
#   bb_next <-  R*(y+B-cc)
#   return(bb_next)
# }

utility <- function(cons_guess, ggamma){
  return((cons_guess^(1-ggamma))/(1-ggamma))
}

m.util <- function(cons_guess, ggamma){
  return(cons_guess^(-ggamma))
}

euler.right <- function(m.util, params, y, mu=0){
  return((params$bbeta*params$R*(m.util%*%params$Py[y,]) + mu)^(-1/params$ggamma))
}

euler.assets.right <- function(qp, mutilp, mutil, y, params, mu=0){
  return(params$bbeta*((mutilp*(qp + params$aalpha*params$y_grid))%*%params$Py[y,])/(mutil-params$phi*mu))
}

cont.value <- function(val_guess, params, y){
  return(val_guess%*%params$Py[y,])
}

value_fun.right <- function(cont.value, cons_guess, params){
  return(utility(cons_guess, params$ggamma) + params$bbeta*cont.value)
}

residuals.complete.tot.model.switching_multiplier_splines_baseline <- function(Coeffs, resid.params){
  
  B_basis = resid.params$B_basis
  
  #  regimes <- resid.params$regimes
  
  params = resid.params$params
  
  val_guess_lhs <- numeric(params$B_dim)
  cons_guess_lhs <- numeric(params$B_dim)
  price_guess_lhs <- numeric(params$B_dim)
  cons_guess_lhs_lim <- numeric(params$B_dim)
  mu_guess <- numeric(params$B_dim)
  
  totpar <- params$Y_dim*params$B_dim*4
  
  Coeffs_val <- matrix(Coeffs[1:(totpar/4)], nrow = params$N_coeffs, ncol = params$Y_dim)
  Coeffs_cons <- matrix(Coeffs[((totpar/4) +1):(2*totpar/4)], nrow = params$N_coeffs, ncol = params$Y_dim)
  Coeffs_price <- matrix(Coeffs[((2*totpar/4) +1):(3*totpar/4)],nrow = params$N_coeffs, ncol = params$Y_dim)
  Coeffs_mu <- matrix(Coeffs[((3*totpar/4) +1):totpar], nrow = params$N_coeffs, ncol = params$Y_dim)
  
  res.v <- matrix(0, nrow = params$N_coeffs, ncol = params$Y_dim)
  res.c <- matrix(0, nrow = params$N_coeffs, ncol = params$Y_dim)
  res.q <- matrix(0, nrow = params$N_coeffs, ncol = params$Y_dim)
  res.mu <- matrix(0, nrow = params$N_coeffs, ncol = params$Y_dim)
  
  for(y in 1:params$Y_dim){
    
    # COEFFICIENTS for the VALUE FUNCTION for a fixed Y
    Vcoeffs_given_y <- Coeffs_val[,y] 
    # COEFFICIENTS for CONSUMPTION for a fixed Y
    Ccoeffs_given_y <- Coeffs_cons[,y] 
    # COEFFICIENTS for the ASSET PRICES for a fixed Y
    Pcoeffs_given_y <- Coeffs_price[,y] 
    # COEFFICIENTS for the MULTIPLIER for a fixed Y
    Mcoeffs_given_y <- Coeffs_mu[,y]
    
    # store debt (ignoring if the constraint is binding)
    b_next <- numeric(params$B_dim) 
    # store debt
    b_next_lim <- numeric(params$B_dim) 
    
    for(b in 1:params$B_dim){
      
      # Calculate the value function V(b, s)
      
      val_guess_lhs[b] <- B_basis[b,]%*%Vcoeffs_given_y 
      
      # Calculate consumption c(b, s)
      
      cons_guess_lhs[b] <- B_basis[b,]%*%Ccoeffs_given_y 
      
      # Calculate asset prices q(b, s)
      
      price_guess_lhs[b] <- B_basis[b,]%*%Pcoeffs_given_y 
      
      # Calculate mu mu(b,s)
      
      mu_guess[b] <- B_basis[b,]%*%Mcoeffs_given_y
      
      # store for which values of b'(b, s) the constraint binds
      
      binds <- numeric(params$B_dim) 
      
      if(cons_guess_lhs[b] <0){
        cons_guess_lhs[b] <- 1e-1
      } 
      
      if(price_guess_lhs[b] < 0){
        price_guess_lhs[b] <- 1e-1
      }
      
      y_current <- params$y_grid[y]
      b_current <- params$B_grid[b]
      
      # Calculate b'(b, s) = R*(y+b-c(b, s))
      # at the beginning b'(b, s) = b'_{constrained}(b, s)
      
      b_next[b] <- next_b(y_current, b_current, cons_guess_lhs[b], params) 
      b_next_lim[b] <- next_b(y_current, b_current, cons_guess_lhs[b], params)
      
      # check if - b'(b, s)/R <= φ*q(b, s)
      
      borr_lim <- b_next[b]/params$R + params$phi*price_guess_lhs[b]
      
      if(borr_lim<0){
        
        b_next_lim[b] <- -params$R*params$phi*price_guess_lhs[b] 
        
        cons_guess_lhs_lim[b] <- y_current+b_current-b_next_lim[b]/params$R
        
        binds[b] <- 1
      }
      
      # THE FOLLOWING CHECKS ARE NEEDED TO ENSURE THAT 
      # b'(b, s) DOESN'T MOVE OUTSIDE ITS DOMAIN [-0.7,0]
      
      
      # b'(b, s) IS UNCONSTRAINED
      
      if(b_next[b]< params$R*(y_current+b_current-cons_guess_lhs[b]) | b_next[b] < (-0.7)){
        
        # Since b'(b, s) is determined by c(b, s)
        # fix c(b,s) such that b'(b, s) = b_MAX = max{grid{b}}
        # then by the resource constraint 
        # c(b, s) = y + b - b_MAX/R
        
        
        cons_guess_lhs[b] <- (params$R*(y_current + b_current) - params$B_grid[length(params$B_grid)])/params$R
        b_next[b] <- next_b(y_current, b_current, cons_guess_lhs[b], params)
        
        
      } else if(b_next[b] > 0){
        
        b_next[b] <- -0.02
        
      }
      
      
      
      # b'(b, s) IS CONSTRAINED
      
      if(b_next_lim[b]< params$R*(y_current+b_current-cons_guess_lhs[b]) | b_next_lim[b] < (-0.7)){
        
        # Since b'_{constrained}(b, s) is determined by q(b, s)
        # fix q(b,s) such that b'(b, s) = b_MAX = max{grid{b}}
        # then by the comp. slackness condition 
        # q(b, s) = -b_MAX/(R*φ)
        
        price_guess_lhs[b] <- - params$B_grid[length(params$B_grid)]/(params$R*params$phi)
        b_next_lim[b] <- -price_guess_lhs[b]*params$R*params$phi
        cons_guess_lhs_lim[b] <- y_current+b_current-b_next_lim[b]/params$R
        
      } else if(b_next_lim[b] > 0){
        
        b_next_lim[b] <- -0.02
        
      }
      
    }
    
    # allocate memory for the policy function in t+1
    # both for the unconstrained case: x_guess_rhs
    # and the constrained case: x_guess_rhs_lim
    
    val_guess_rhs <- numeric(params$Y_dim)
    val_guess_rhs_lim <- numeric(params$Y_dim)
    
    cons_guess_rhs <- numeric(params$Y_dim)
    cons_guess_rhs_lim <- numeric(params$Y_dim)
    
    price_guess_rhs <- numeric(params$Y_dim)
    price_guess_rhs_lim <- numeric(params$Y_dim)
    
    
    # Calculate Natural Splines using the grids for
    
    # b'(b, s)
    
    B_next_basis <- bs(b_next, knots = seq(params$b_min, params$b_max, length.out=params$N_coeffs)[-c(1:2,params$N_coeffs-1,params$N_coeffs)],
                       intercept = T, degree=3, Boundary.knots = c(-0.7, 0))
    
    # b'_{constrained}(b, s)
    
    B_next_basis_lim <- bs(b_next_lim, knots = seq(params$b_min, params$b_max, length.out=params$N_coeffs)[-c(1:2,params$N_coeffs-1,params$N_coeffs)],
                           intercept = T, degree=3, Boundary.knots = c(-0.7, 0))
    
    
    for(bp in 1:params$B_dim){
      
      
      for(yp in 1:params$Y_dim){
        
        # COEFFICIENTS for the VALUE FUNCTION for a fixed Y'
        Vcoeffs_given_y_next <- Coeffs_val[,yp]
        # COEFFICIENTS for CONSUMPTION for a fixed Y'
        Ccoeffs_given_y_next <- Coeffs_cons[,yp]
        # COEFFICIENTS for the ASSET PRICES for a fixed Y'
        Pcoeffs_given_y_next <- Coeffs_price[,yp]
        
        
        
        # Calculate the value function V(b', s')
        
        val_guess_rhs[yp] <- B_next_basis[bp,]%*%Vcoeffs_given_y_next
        
        # Calculate the value function V_{constrained}(b', s')
        
        val_guess_rhs_lim[yp] <- B_next_basis_lim[bp,]%*%Vcoeffs_given_y_next
        
        
        
        # Calculate consumption c(b', s')
        
        cons_guess_rhs[yp] <- B_next_basis[bp,]%*%Ccoeffs_given_y_next
        
        # Calculate consumption c_{constrained}(b', s')
        
        cons_guess_rhs_lim[yp] <- B_next_basis_lim[bp,]%*%Ccoeffs_given_y_next
        
        
        
        # Calculate asset prices c(b', s')
        
        price_guess_rhs[yp] <- B_next_basis[bp,]%*%Pcoeffs_given_y_next
        
        # Calculate asset prices q_{constrained}(b', s')
        
        price_guess_rhs_lim[yp] <- B_next_basis_lim[bp,]%*%Pcoeffs_given_y_next
        
        
        
        if(price_guess_rhs[yp] < 0){
          price_guess_rhs[yp] <- 1e-1
        }
        
        if(price_guess_rhs_lim[yp] < 0){
          price_guess_rhs_lim[yp] <- 1e-1
        }
        
        if(cons_guess_rhs[yp] <0){
          cons_guess_rhs[yp] <- 1e-1
        } 
        
        if(cons_guess_rhs_lim[yp] <0){
          cons_guess_rhs_lim[yp] <- 1e-1
        } 
        
        
      }
      
      # Calculate marginal utility in t+1 in the UNCONSTRAINED and CONSTRAINED CASE
      
      m.util.p <- m.util(cons_guess_rhs, params$ggamma)
      
      m.util.p_lim <- m.util(cons_guess_rhs_lim, params$ggamma)
      
      # Calculate marginal utility in t
      
      m.util.current <- m.util(cons_guess_lhs[bp], params$ggamma)
      
      m.util.current_lim <- m.util(cons_guess_lhs_lim[bp], params$ggamma)
      
      # CALCULATE RESIDUALS IN THE UNCONSTRAINED CASE
      
      if(!binds[bp]==1){
        
        # Define multipliers using Garcia-Zangwill trick
        
        alpha_neg <- max(-mu_guess[bp], 0)^(2)
        
        alpha_pos <- max(mu_guess[bp], 0)^(2)
        
        
        euler.rhs <- euler.right(m.util.p, params, y, mu = alpha_pos)
        
        euler.assets.rhs <- euler.assets.right(price_guess_rhs, m.util.p, m.util.current, y, params, mu = alpha_pos)
        
        cont.val.p <- cont.value(val_guess_rhs, params, y)
        
        val.fun.rhs <- value_fun.right(cont.val.p, cons_guess_lhs[bp], params)
        
        res.v[bp,y] <- val.fun.rhs - val_guess_lhs[bp]
        
        res.c[bp,y] <- euler.rhs - cons_guess_lhs[bp]
        
        res.q[bp,y] <- euler.assets.rhs - price_guess_lhs[bp]
        
        res.mu[bp,y] <-  b_next[bp]/params$R + params$phi*price_guess_lhs[bp] - alpha_neg
        
      }
      
      
      
      # CALCULATE RESIDUALS IN THE CONSTRAINED CASE
      
      if(binds[bp]==1){
        
        # Define multipliers using Garcia-Zangwill trick
        
        alpha_neg <- max(-mu_guess[bp], 0)^(2)
        
        alpha_pos <- max(mu_guess[bp], 0)^(2)
        
        
        euler.rhs <- euler.right(m.util.p_lim, params, y, mu = alpha_pos)
        
        euler.assets.rhs <- euler.assets.right(price_guess_rhs_lim, m.util.p_lim, m.util.current, y, params, mu = alpha_pos)
        
        cont.val.p <- cont.value(val_guess_rhs_lim, params, y)
        
        val.fun.rhs <- value_fun.right(cont.val.p, cons_guess_lhs[bp], params)
        
        res.v[bp,y] <- val.fun.rhs - val_guess_lhs[bp]
        
        res.c[bp,y] <- euler.rhs - cons_guess_lhs[bp]
        
        res.q[bp,y] <- euler.assets.rhs - price_guess_lhs[bp]
        
        res.mu[bp,y] <-  b_next_lim[bp]/params$R + params$phi*price_guess_lhs[bp] - alpha_neg
        
      }
      
      
      
    }
    
    
  }
  
  # Vectorize the residuals (to obtain K equations [the residuals] in K uknowns [the coefficients] )
  
  restot <- vec(cbind(res.v, res.c, res.q, res.mu))
  
  return(as.numeric(restot))
  
}
