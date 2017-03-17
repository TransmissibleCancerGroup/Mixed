# Simulate 1 sample of `size` mutations from a mixture of
# `proportion` somatic : (1 - `proportion`) germline mutation
# `signature`s, using a multinomial model
simulate_multinom <- function(size, proportion, signatures) {
    exposures <- c(1-proportion, proportion)
    expectation <- signatures %*% exposures
    rmultinom(1, size, expectation)
}

# Simulate 1 sample of `size` mutations from a mixture of
# `proportion` somatic : (1 - `proportion`) germline mutation
# `signature`s, using a Poisson model
simulate_poisson <- function(size, proportion, signatures) {
    exposures <- c(1-proportion, proportion)
    expectation <- size * signatures %*% exposures
    sample <- as.matrix(sapply(expectation, function(x) {rpois(1, x)}),
                        nrow=nrow(signatures))
    rownames(sample) <- rownames(signatures)
    sample
}

# Simulate `n` samples according to a multinomial or Poisson
# model, by passing simulate_multinom or simulate_poisson as
# `samplefun`
multisample <- function(n, samplefun, ...) {
    out <- matrix(nrow = 96, ncol = n)
    for (i in 1:n) out[, i] <- samplefun(...)
    out
}

# Barplot result matrix in `sample_mat` with error bars
# (also works with single samples, but error bars are meaningless)
plot_multisample <- function(sample_mat) {
    barheights <- rowMeans(sample_mat)
    barheight_sd <- apply(sample_mat, 1, sd)
    bars <- barplot(barheights, ylim = c(0, max(barheights+barheight_sd*1.96)))
    segments(bars, barheights - barheight_sd*1.96, bars, barheights + barheight_sd*1.96)
    arrows(bars, barheights - barheight_sd*1.96, bars, barheights + barheight_sd*1.96, 
           angle=90, code=3, length=0.05)
}

noise <- function(size) {
    rpois(96, size)
}
