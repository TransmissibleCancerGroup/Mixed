require(rstan)

poisson_model <- "
data {
    int n; // number of data points for m (samples)
    int m[n]; // number of mutations per sample in short time period (t)
    int N; // number of data points for M
    int M[N]; // number of mutations since origin (long time period, T)
    
    // prior parameters:
    real<lower=0> a;
    real<lower=0> b;
    real<lower=0> min_permitted_lower_bound;
    real<lower=0> max_permitted_upper_bound;
    real<lower=0> expon_rate;
}
parameters {
    real<lower=0> rate; // rate parameter of the Poisson process (mutations per unit time [year])
    real<lower=0> t; // short time interval in which mutations are generated (years)
    real<lower=min_permitted_lower_bound,
    upper=max_permitted_upper_bound> T; // long time interval (years)
}
model {
    // prior for 't' is Beta value ( [0,1) )
    t ~ beta(a, b);
    
    // prior for 'T' is Uniform ( [min_permitted, max_permitted) )
    T ~ uniform(min_permitted_lower_bound, max_permitted_upper_bound);
    
    // prior on mutation rate is Exponential ( [0, inf) )
    // NB: Expectation(mutation rate) = 1/expon_rate
    rate ~ exponential(expon_rate);
    
    // likelihood of 'm' is Poisson, rate scaled by 't'
    m ~ poisson(rate * t);
    
    // likelihood of 'M' is Poisson, rate scaled by 'T'
    M ~ poisson(rate * T);
}
"

# Compile the model
my_model <- rstan::stan_model(model_code = poisson_model,
                              model_name = "CTVT poisson process")

# List of model parameters + input data
dat <- list(
    n = 2,  # number of post-MRCA mutation estimates
    N = 1,  # number of pre-MRCA mutation estimates
    m = c(23, 27),  # post-MRCA mutation estimates
    M = array(221370),  # pre-MRCA mutation estimates (use array to pass a single value)
    a = 1,  # parameter a of beta prior on 't'
    b = 1,  # parameter b of beta prior on 't'
    min_permitted_lower_bound = 200,  # lower bound on uniform prior on 'T'
    max_permitted_upper_bound = 25000,  # upper bound on uniform prior on 'T'
    expon_rate = 0.025  # parameter of exponential prior on mean mutation rate (corresponds to ~= 40)
)

# Get MCMC samples from the model
trace <- rstan::sampling(my_model, data = dat, iter = 50000, chains = 4)

# List parameter estimates
summary(trace, prob = c(0.025, 0.975))

# Plot some stuff
require(bayesplot)
bayesplot::mcmc_areas(as.array(trace), pars = "t", prob = .95)
bayesplot::mcmc_areas(as.array(trace), pars = "T", prob = .95)
bayesplot::mcmc_areas(as.array(trace), pars = "rate", prob = .95)
