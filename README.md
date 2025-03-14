# globalSolve

The **globalSolve** package provides functions to solve DSGE models by means of (global) __projection methods__. Providing a simple and intuitive approach, the package requires a simple specification of the firs order conditions and the nature of the variables to solve the model, reducing the cost of writing hundreds of lines of code to less then 20 rows.

The solution strategies already implemented are:
  - **Tensorial Chebyshev polynomials** (for smooth problems with a reduced number of exogenous processes)
  - **Finite Elements** (for models with non-linearities introduced by inequality constraints)
  - **Asintropic grids** (for models with a large amount of exogenous processes)

## General problem

The problem usually defined by the optimal conditions of the model can be summarized by the following set of functions

$$
\begin{aligned}
z_{t}=h(z_{t-1}, \varepsilon_{z,t}) & \text{Exogenous processes} \\
x_{t+1} = g(x_{t}, z_{t}) & \text{Endogenous state variables} \\
y_{t} = f(x_{t}, z_{t}) & \text{Endogenous variables} 
\end{aligned}
$$

To find a solution to such class of models **globalSolve** approximates the functions $\{g,f\}$ by means of basis functions.


# Example

# Problem:

$$
\begin{aligned}
V(k_t, z_t) = \max_{c_t, l_t} \frac{(c_t^{\tau} (1 - l_t)^{1 - \tau})^{1 - \eta}}{1 - \eta} + \beta \mathbb{E_{t}} V(k_{t+1}, z_{t+1}) \\
\text{s.t. } k_{t+1} = e^{z_t} k_t^{\alpha} l_t^{1 - \alpha} + (1 - \delta) k_t - c_t \\
z_t = \rho z_{t-1} + \epsilon_t \\
\end{aligned}
$$

Assume for $l_{t}$ and $V_{t}$ decision rules $g_{x}$ of the form:

$$
g_{x}(k_t, z_t) = \displaystyle\sum_{i=0}^{K} \theta_{x,i} \varphi_{i}(k_t, z_t)
$$

Then the problem to solve assumes the following structure:

$$
\begin{aligned}
\frac{(c_t^{\tau} (1 - l_t)^{1 - \tau})^{1 - \eta}}{1 - \eta} + \beta \mathbb{E_{t}} g_{V}(k_{t+1}, z_{t+1}) - g_{V}(k_t, z_t) = 0\\
\frac{\tau}{1-\tau} (1-\alpha) e^{z_{t}} k_{t}^{\alpha} g_{l}(k_t, z_t)^{-\alpha} (1-g_{l}(k_t, z_t)) - c_{t} = 0 \\
k_{t+1} = e^{z_{t}} k_{t}^{\alpha} g_{l}(k_t, z_t)^{1-\alpha} + (1-\delta) k_{t} - c_{t}\\
z_{t} = \rho_{z} z_{t-1} + \varepsilon_{z,t}
\end{aligned}
$$

To represent the above group of conditions the package requires the following strcture:
  - Exogenous state variables: ```exo_x```
  - Endogenous state variables: ```f_State_Endo_x```
  - Endogenous variables: ```f_Endo_x```

```{r}
FOCs <- list(f_Endo_Value_fun= (c^(tau)
)
```


