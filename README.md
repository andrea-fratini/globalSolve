# globalSolve

The globalSolve package provides functions to solve DSGE models by means of (global) projection methods. Providing a simple and intuitive approach, the package requires a simple specification of the firs order conditions and the nature of the variables to solve the model, reducing the cost of writing hundreds of lines of code to less then 20 rows.

The solution strategies already implemented are:
  - Tensorial Chebyshev polynomials (for smooth problems with a reduced number of exogenous processes)
  - Finite Elements (for models with non-linearities introduced by inequality constraints)
  - Asintropic grids (for models with a large amount of exogenous processes)

## General problem

The problem usually defined by the optimal conditions of the model can be defined by the following set of functions

$$
\begin{aligned}
z_{t}=h(z_{t-1}, \varepsilon_{z,t}) & \text{Exogenous processes} \\
x_{t+1} = g(x_{t}, z_{t})\\
y_{t+1} = f(x_{t+1}, z_{t+1})
\end{aligned}
$$
  
