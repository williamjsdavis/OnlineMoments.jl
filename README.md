[![DOI](https://zenodo.org/badge/555072881.svg)](https://zenodo.org/badge/latestdoi/555072881)

# OnlineMoments.jl
Calculating conditional moments (and variance) from data streams
`OnlineMoments.jl` is a `Julia` package enabling the esimation of conditional moments and conditional variance from empirical streamed time-series data.
This package enables the estimation of drift and diffusion functions using both histogram- and kernel-based regression.

This package accompanies the paper: "Reconstruction of Stochastic Dynamics from Large Datasets" Davis, William (submitted)

## Background

Consider the scalar stochastic differential equation

$`\frac{d}{dt}X(t) = f(X) + g(X)\Gamma(t).`$

The evolution of this process can be described by a Fokker-Planck equation

$`\frac{\partial}{\partial t} p(x,t|x^\prime,t^\prime) = \bigg[-\frac{\partial}{\partial x} D^{(1)}(x) + \frac{\partial^2}{\partial x^2} D^{(2)}(x) \bigg] p(x,t|x^\prime,t^\prime)`$

which contain the Kramers-Moyal (KM) coefficients 

$`D^{(k)}(x) = \lim_{\tau\rightarrow0} \frac{1}{n!\tau} \int_{-\infty}^{\infty} \big[x^\prime - x\big]^k p(x^\prime,t+\tau|x,t)\ dx^\prime.`$

Here $`k=1`$ is the drift coefficient and $`k=2`$ is the diffusion coefficient. These expressions contain moments of the conditional probability density, or "conditional moments"

$`M^{(k)}(\tau,x) = \int_{-\infty}^\infty [x^\prime - x\big]^k p(x^\prime, t + \tau| x,t)\ dx^\prime.`$

These conditional moments are typically estimated in an offline fashion

$`\hat{M}^{(k)}_{ij} = \frac{\sum\limits_{n=1}^{N-i} K_h(\mathcal{X}_j - X_n)\big[X_{n+i} - X_n\big]^k}{\sum\limits_{n=1}^{N-i} K_h(\mathcal{X}_j - X_n)}.`$

In this work, I present *online* formulae for sequential updating

$`\hat{M}^{(k)}_{ij}\big|_{N} = \hat{M}^{(k)}_{ij}\big|_{N-1} + K_h(\mathcal{X}_j - X_{N-i})\times\left(\left[X_N - X_{N-i}\right]^k - \hat{M}^{(k)}_{ij}\big|_{N-1}\right)\Big/W_{ij}\big|_{N},`$

$`W_{ij}\big|_{N} = W_{ij}\big|_{N-1} + K_h\left(\mathcal{X}_j - X_{N-i}\right).`$

The online method, which I call "Online Kernel-Based Regression (OKBR)" scales $`\mathcal{O}(1)`$ in space, compared with Kernel-Based Regression (KBR) which scales as $`\mathcal{O}(N)`$.

<img width="500" alt="Screenshot 2023-07-01 at 3 27 10 PM" src="https://github.com/williamjsdavis/OnlineMoments.jl/assets/38541020/0556d811-17ac-4565-a436-423419371e44">

See the paper for further details.

## Version

- Version 0.1.0 - First release

