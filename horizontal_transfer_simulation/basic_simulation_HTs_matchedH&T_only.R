# Horizontal tranfer events - basic simulation study
# Very basic simulation:
#    we observed 17 HT events among 18 haplotypes, that occur in the frequencies given by `totals`
# Simulation, single rep: 
#   Randomly assign 17 HT events to our 18 haplotypes, weighted by their frequency of occurrence


############
# Input data
############
haplotypes <- c("A1a1", "B1", "A1d1a", "C1", "C2", "A1f", "A1h", "A2", "A1e", "A6", "A1b", "A1c", "A1d2", "A1d1", "A4", "A5", "B2", "A1")
totals <- c(136, 96, 94, 39, 24, 18, 14, 11, 9, 8, 6, 6, 5, 3, 3, 2, 2, 1)
cases <- c(5, 1, 10, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
names(totals) <- haplotypes
names(cases) <- haplotypes


################
# Simulated data
################
set.seed(76) # <- Always same results for given random number seed
simulation_reps <- 10000

# Generate simulated data
counts <- matrix(0, nrow=simulation_reps, ncol=length(haplotypes), dimnames=list(1:simulation_reps, haplotypes))
for (i in 1:nrow(counts)) {
    .sample <- sample(haplotypes, size = sum(cases), replace = TRUE, prob = totals)
    .counts <- table(.sample)
    counts[i, names(.counts)] <- .counts
}


######################
# Function definitions
######################
draw_hist <- function(data, nbins, observed = NULL, ...) {
    hist(data, breaks = seq(0, 1, length.out = nbins + 1),
         xlim = c(0, 1), col = "dodgerblue", border = "white",
         xlab = "Proportion of all cases",
         ...)
    if (!is.null(observed)) {
        abline(v = observed, col = "red", lty = 2)
        frac.gt <- sum(data >= observed) / length(data)
        legend("topright", legend = c("Observed proportion of cases",
                                      paste0("(", 100*frac.gt,"% of simulations equal or exceed this)")),
               lty = c(2, 0), col = c("red", "white"))
        text(x = observed, y = -.1, formatC(observed), col = "red")
        cat(paste("Assuming all haplotypes have the same chance of generating a horizontal transfer,\n",
                     "and exactly", sum(cases), "events occur over the cohort,\n",
                     sprintf("the empirical p-value = %.4f\n", frac.gt)))
    }
}

#####################
# Analysis of results
#####################
par(mfrow = c(2, 1))

# Look at the highest proportion of the total number of cases that occur in any single haplotype
cat("Looking at highest proportion in any haplotype...\n")
max.prop <- apply(t(apply(counts, 1, function(r) r/sum(r))), 1, max)
draw_hist(max.prop, nbins = 18, observed = max(cases)/sum(cases),
          main = "Highest proportion of cases in any single haplotype")

# Look at the proportion in our haplotype of interest, A1d1a
cat("Looking at our favourite haplotype...\n")
prop.a1d1a <- t(apply(counts, 1, function(r) r/sum(r)))[, "A1d1a"]
draw_hist(prop.a1d1a, nbins = 18, observed = cases["A1d1a"]/sum(cases),
          main = "Proportion of cases in haplotype A1d1a")


# Uncomment to look at other haplotypes...
# prop.b1 <- t(apply(counts, 1, function(r) r/sum(r)))[, "B1"]
# draw_hist(prop.b1, nbins = 19, observed = cases["B1"]/sum(cases),
#           title = "Proportion of cases in haplotype B1")
# 
# prop.a1 <- t(apply(counts, 1, function(r) r/sum(r)))[, "A1"]
# draw_hist(prop.a1, nbins = 19, observed = cases["A1"]/sum(cases),
#           title = "Proportion of cases in haplotype A1")
