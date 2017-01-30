library(gtools)
library(coda)
library(progress)

likelihood_fn <- function(x, ref.spectra, obs.counts) {
  lc <- linear.combination(ref.spectra, x)
  log.likelihood <- dmultinom(obs.counts, prob = lc, log = TRUE)
  return(log.likelihood)
}

propose <- function(x) {
  r <- runif(1)
  s <- sample(length(x), 2)
  t <- sum(x[s])
  x[s[1]] <- r * t
  x[s[2]] <- (1-r) * t
  x
}

linear.combination <- function(spectra, x) {
  stopifnot(ncol(spectra) == length(x))
  lc <- rowSums(t(t(spectra) * x))
  stopifnot(all(lc > 0))
  return (lc)
}

mcmc.iteration <- function(param, alpha, signatures, data, prior_only = FALSE) {
  proposal <- propose(param)
  if (!prior_only) {
    current <- log(ddirichlet(param, alpha)) + likelihood_fn(param, signatures, data)
    alt <- log(ddirichlet(proposal, alpha)) + likelihood_fn(proposal, signatures, data)
  }
  else {
    current <- log(ddirichlet(param, alpha))
    alt <- log(ddirichlet(proposal, alpha))
  }
  if (alt > current) return(proposal)  # Automatically accept better proposal
  else {
    r <- runif(1)
    if (log(r) < (alt - current)) return(proposal)  # Accept worse proposal with probability (alt/current)
  }
  return(param)  # Reject proposed param, return original param
}

summarise.chain.weights <- function(chain.weights, hpdwidth = 0.95, title = "Signature exposures") {
    chain.weights.summary <- colMeans(chain.weights)
    names(chain.weights.summary) <- 1:length(chain.weights.summary)
    error <- HPDinterval(as.mcmc(chain.weights), prob=hpdwidth)
    bars <- barplot(chain.weights.summary, ylim = c(0, max(error[, 2]*1.05)), main = title,
                    col = ifelse(error[, 1] > 1e-3, "dodgerblue3", "grey90"), cex.names=0.9)
    arrows(bars, error[, 1], bars, colMeans(chain.weights), angle=90, code=1, length=0.05)
    arrows(bars, error[, 2], bars, colMeans(chain.weights), angle=90, code=1, length=0.05)
    legend("topright", legend = paste(prettyNum(100*hpdwidth), "% HPD > 0", sep = ""), fill = "dodgerblue3")
    chain.weights.summary
}

run.mcmc <- function(nsamples, thin, burnin, data, signatures, alpha = NULL, prior_only = FALSE) {
  if(!require(gtools)) stop("Require package gtools")
  if(!require(coda)) stop("Require package coda")
  if(!require(progress)) stop("Require package progress")

  if (is.null(alpha)) alpha <- rep(1, ncol(signatures))

  NSAMPLES <- nsamples  # Record this many samples...
  THIN <- thin          # ... separated by this many iterations ...
  BURNIN <- burnin      # ... after this amount of burn-in

  chain.lik <- rep(0, NSAMPLES)   # store sampled likelihoods,
  chain.prior <- rep(0, NSAMPLES) # sampled priors,
  chain.post <- rep(0, NSAMPLES)  # sampled posterior probabilities
  chain.weights <- matrix(0, nrow=NSAMPLES, ncol=ncol(signatures)) # and sampled weights
  colnames(chain.weights) <- colnames(signatures)

  # PRIOR DISTRIBUTION
  # ==================
  # Our parameter is a vector of 30 weights, which sum to 1.
  # A natural prior to use is the Dirichlet: Prior ~ Dir(alpha).
  # alpha is a vector of length 30. The larger each value of
  # alpha[i], the more the prior favours weight[i].

  # Initial guess is a sample from the prior
  param <- as.vector(rdirichlet(1, alpha))

  burninbar <- progress_bar$new(total = BURNIN, format = "Burn-in: [:bar] :current / :total ")
  for (i in 1:BURNIN) {
    param <- mcmc.iteration(param, alpha, signatures, data, prior_only)
    if (i %% 100 == 0) {
      burninbar$tick(100)
    }
  }

  samplebar <- progress_bar$new(total = NSAMPLES, format = "Sampling: [:bar] :current / :total ")
  for (j in 1:NSAMPLES) {
    for (i in 1:THIN) {
      param <- mcmc.iteration(param, alpha, signatures, data, prior_only)
    }
    chain.lik[j] <- likelihood_fn(param, signatures, data)
    chain.prior[j] <- log(ddirichlet(param, alpha))
    chain.post[j] <-  chain.prior[j] + chain.lik[j]
    chain.weights[j, ] <- param
    samplebar$tick(1)
  }
  probs = as.mcmc(matrix(c(chain.prior, chain.lik, chain.post), nrow=NSAMPLES, ncol=3, byrow = FALSE))
  colnames(probs) <- c("Prior", "Likelihood", "Posterior")
  list(weights = as.mcmc(chain.weights),
       probs = probs)
}

if (!interactive()) {

  library(optparse)
  option_list = list(
    make_option(c("--prior"), type="character", default="w",
                help="Prior strength - choices are w or s, i.e. weak or strong", metavar="character"),
    make_option(c("--plot"), type="character", default=NULL,
                help="Output image file name. Must end .pdf", metavar="character"),
    make_option(c("--data-index"), type="integer", default=NULL,
                help="Which of the 34 data subsets to use", metavar="integer"),
    make_option(c("--burnin"), type="integer", default=500000,
                help="Length of burn-in chain", metavar="integer"),
    make_option(c("--samples"), type="integer", default=2000,
                help="Number of samples to collect after burn-in", metavar="integer"),
    make_option(c("--thin"), type="integer", default=1000,
                help="Number of chain iterations between samples", metavar="integer")
  )
  parser <- OptionParser(option_list=option_list)
  arguments <- parse_args(parser)
  data.idx <- as.integer(arguments$`data-index`)
  if(is.na(data.idx)) stop(parser@usage)
  if (data.idx < 1 | data.idx > ncol(CTVT.COUNT.DATA)) {
    stop(paste("data-index must be between 1 and", ncol(CTVT.COUNT.DATA)))
  }

  # load("../cosmic.signatures.RData")    # COSMIC.SIGNATURES
  # load("../ctvt.count.data.RData")      # CTVT.COUNT.DATA

  cos.sim <- function(x, y) x %*% y / sqrt(x%*%x * y%*%y)

  mydata <- CTVT.COUNT.DATA[, data.idx]
  mydataname <- colnames(CTVT.COUNT.DATA)[data.idx]
  mysignatures <- COSMIC.SIGNATURES

  weak_prior = rep(1, ncol(mysignatures))
  medium_prior = c(2, 1, 1, 1, 2, 1, 2, 1, 1, 1,
                   0.1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 0.1)
  strong_prior = c(100, 2, 2, 2, 100, 2, 100, 2, 2, 2,
                   0.1, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                   2, 2, 2, 2, 2, 2, 2, 2, 2, 0.1)

  if (arguments$`prior` == "w") {
    myprior <- weak_prior
    mypriorname <- "weak prior"
  }

  else if(arguments$`prior` == "m") {
    myprior <- medium_prior
    mypriorname <- "medium prior"
  }
  else if(arguments$`prior` == "s") {
    myprior <- strong_prior
    mypriorname <- "strong prior"
  }
  else stop("Prior must be one of w, m or s")

  posterior <- run.mcmc(arguments$`samples`, arguments$`thin`, arguments$`burnin`, mydata, mysignatures, myprior)
  prior_only <- run.mcmc(arguments$`samples`, arguments$`thin`, arguments$`burnin`, mydata, mysignatures, myprior, prior_only = TRUE)

  pdf(arguments$`plot`, width = 8, height = 10)
  par(mfrow = c(2, 1))
  summarise.chain.weights(posterior$weights, paste(mydataname, "Posterior,", mypriorname))
  summarise.chain.weights(prior_only$weights, paste(mydataname, "Prior only,", mypriorname))

  colours <- c(rep("blue", 16), rep("black", 16), rep("red", 16), rep("grey", 16), rep("green", 16), rep("pink", 16))
  myfit <- linear.combination(mysignatures, colMeans(posterior$weights))
  mycos.sim <- cos.sim(myfit, mydata / sum(mydata))
  barplot(mydata, main = paste(mydataname, "counts"),
          las=2, cex.names=0.5, family="mono", border=NA, col = colours)
  barplot(myfit, main = paste(mypriorname, "fitted spectrum: cosine similarity =", mycos.sim),
          las=2, cex.names=0.5, family="mono", border=NA, col = colours)

  par(mfrow = c(5, 2))
  plot(posterior$weights, auto.layout = FALSE)
  plot(prior_only$weights, auto.layout = FALSE)

  dev.off()
}
