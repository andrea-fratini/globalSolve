# globalSolve

Structure of the file submitted by the user (as in Dynare):

endogenous variables <- $c("x", "y", "z", "z_{m1}", "x_1", "y_1")$

exogenous variables <- c("\epsilon_z")

parameters <- c(alpha=num, gamma=num, rho_z=num, sigma_z=num)

## Equations

EV_x <- x ~ E[x_1^(alpha) | F] 

EV_y <- y ~ E[x_1^(gamma) | F]

PROC_z <- z ~ rho_z*z_m1 + epsilon_z

# Model objective

Find functions x(z), y(z) that solve the integral equation x = int x_1^(alpha) dfz

The process z is discretized to have K states (a markovchain) on a state space of K numbers (which are reals)

then some techinques among Chebyshev polynomials (code for collocation already implemented, to implement time iteration - to move in c++) , Smolyak polynomials (to implement functions, collocation and time iteration - to move in c++) to obtain solution functions

Methods needed:

Residual function -> to find the zeros of the infinite dimensional problem

Residual error -> to evaluate the model performance

Matching moments -> to find the optimal value of the parameters in a D-dimensional grid

Simulation -> to generate a path for the joint process

Impulse Response -> to generate impulse response functions

Method comparison -> uses previous techinques with same parameters to solve the model and returns policy functions and residual errors

# Computational resources

Be careful to reduce as much as possible the computational cost, i.e. if passing some R functions to c++ try to make less reads from the R environment as possible
