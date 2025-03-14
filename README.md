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
\frac{(c_t^{\tau} (1 - g_{l}(k_t, z_t))^{1 - \tau})^{1 - \eta}}{1 - \eta} + \beta \mathbb{E_{t}} g_{V}(k_{t+1}, z_{t+1}) - g_{V}(k_t, z_t) = 0\\
\frac{\tau}{1-\tau} (1-\alpha) e^{z_{t}} k_{t}^{\alpha} g_{l}(k_t, z_t)^{-\alpha} (1-g_{l}(k_t, z_t)) - c_{t} = 0 \\
k_{t+1} = e^{z_{t}} k_{t}^{\alpha} g_{l}(k_t, z_t)^{1-\alpha} + (1-\delta) k_{t} - c_{t}\\
z_{t} = \rho_{z} z_{t-1} + \varepsilon_{z,t}
\end{aligned}
$$

To represent the above group of conditions the package requires different specifications concerning:
  - **Nature of the variables**
    - Exogenous ```exogenous=c("z")```
    - State Endogenous ```state_endogenous=c("k")```
    - Endogenous ```endogenous=c("V", "l", "c")```
    - Shocks ```shocks=c("epsz")```
  - **Variables timing**
    - Backwards $t-k$:  ```x_bk```
    - In time $t$:  ```x```
    - Forward $t+k$ ```x_fk```
  - **Parameters** ```list(tau=0.4, eta=0.3, beta=0.95, alpha=0.6, delta=0.1, rho_=c(z=0.6), sigma_=c(z=0.02), mu_=c(z=0)) ```
  - **Technical parameters**
      - Grids informations for the Exogenous variables: ```list(z=c(n_points, width_of_the_grid),  ...)```
      - Grids informations for the Endogenous state variables: ```list(k=c(min(grid), max(grid), n_points), ...)```
      - Basis type and required parameters: ``basis_type=list(type="Chebyshev", ..parameters to be defines..)``
  - **Equations**
    - Exogenous state variables: ```exo_x```
    - Endogenous state variables: ```f_State_Endo_x```
    - Endogenous variables: ```f_Endo_x```
  
The complete code should look like this:

```{r}
variables <- list(exogenous=c("z"), state_endogenous=c("k"),
                  endogenous=c("V", "l", "c"), shocks=c("epsz"))

parameters <- list(tau=0.4, eta=0.3, beta=0.95, alpha=0.6, delta=0.1,
                   rho_=c(z=0.6), sigma_=c(z=0.02), mu_=c(z=0)) # mu_ is required if the process is not cenetered in 1 (or its log(*) in 0)

technical_params <- list(exo_grids=list(z=c(n_points, width_of_the_grid)),
                         state_endo_grids=list(b=c(-0.7, 0, 10), r=c(0.9, 1, 10)),
                         basis_type=list(type="Chebyshev", ..parameters to be defines..))

FOCs <- list(f_Endo_Value_fun = V ~ (c^(tau) (1-l)^(1-tau))^(1-eta) / (1-eta) + beta * V_f1
             f_Endo_Euler_cons = c ~ (tau) / (1-tau) exp(z) k^(alpha) l^(-alpha) * (1-l)
             f_State_Endo_k = k_f1 ~ exp(z) k^(alpha) l^(-alpha) + (1-delta) * k - c
             exo_z = z_f1 ~ rhoz * z + epsz
)
```


