## SVmodelRcppSMC

## Some more details:
This is an example of an SV model to be estimated via PMMH using RcppSMC in a
as a package. Also, RcppSMC is used to obtain for log-likelihood estimates 
for a standard toy SV model applied to real data or simulated data.

A pure R implementation of Particle Gibbs with and without ancestor sampling for
the same model is provided too.

The SV model equations (measurements and latent state transition) are given as:

1. $`X_t|(X_{t-1} = x_{t-1})\sim\mathcal{N}\left(x_t|\phi x_{t-1}, \sigma^2\right)\;,`$
2. $`Y_t|(X_t = x_t) \sim \mathcal{N}\left(y_t|0, \beta^2\exp(x_t)\right)\;.`$

Assuming that the variance parameters $\sigma^2$ and $\beta^2$ are unknown, but
$\phi$ is fixed, we consider a Bayesian setting with standard inverse Gamma 
priors on the parameters:

3. $\sigma^2\sim\mathcal{IG}(a=0.01, b=0.01)\;,$
4. $\beta^2\sim\mathcal{IG}(a=0.01, b=0.01)\;.$

This prior setup is conjugate to the model from 1. and 2. such that full 
conditional Gibbs blocks are obtained with distribution for the parameter in 
closed form according to:

5. $p\left(\sigma^2|x_{0:T}, y_{1:T}\right)=p(\sigma^2|x_{0:T})=\mathcal{IG}\left(\sigma^2|a+\frac{T}{2}, b+\frac{1}{2}\sum_{t=1}^T(x_t-\phi x_{t-1})^2\right)\;,$
6. $p\left(\beta^2|x_{0:T}, y_{1:T}\right)=\mathcal{IG}\left(\beta^2|a+\frac{1}{2}, b+\frac{1}{2}\sum_{t=1}^T\exp\left(-x_t\right)y_t^2\right)\;.$
