# Custom versions of hdp::hdp_extract_comp_single and hdp::hdp_extract_comp_multi
# from the 'hdp' R package. These custom versions allow the significance level
# of the Bayesian credible intervals to be adjusted through the parameter 'sig.level'

#' function to merge cols of a matrix according to new (numeric) column labels
#' ans will be sorted by numeric column name
merge_cols <- function(mx, labels){
  colnames(mx) <- labels
  num_rows <- nrow(mx)
  # init matrix of 0s with new column names
  ans <- matrix(0, nrow=num_rows, ncol=length(unique(labels)),
                dimnames=list(c(1:num_rows), sort(as.numeric(unique(labels))))
  )

  for (label in colnames(ans)){
    ans[, label] <- rowSums(matrix(mx[, which(colnames(mx) == label)],
                                   nrow=num_rows))
  }
  return(ans)
}

#' function to merge elements of a list according to new (numeric) names
#' ans will be sorted by numeric element names (is this still needed?)
merge_elems <- function(lst, labels){
  names(lst) <- labels
  # init list with new names
  ans <- vector("list", length(unique(labels)))
  names(ans) <- sort(as.numeric(unique(labels)))

  for (label in names(ans)){
    ans[[label]] <- Reduce("+", lst[which(labels == label)])
  }
  return(ans)
}

#' Extracts components from single MCMC chain
#' @export
#' @importFrom "hdp" comp_categ_counts sampling_seed numcluster clust_categ_counts clust_dp_counts numcateg final_hdpState numdp hdp_settings dp pseudoDP comp_dp_counts
#' @importFrom "methods" as validObject
hdp_extract_comp_single <- function(chain, sig.level = 0.95, cos.merge=0.90){

  # input checks
  if (class(chain) != "hdpSampleChain") {
    stop("chain must have class hdpSampleChain")
  }
  if (!validObject(chain)) stop("chain not a valid hdpSampleChain object")
  if (class(cos.merge) != "numeric" | cos.merge >=1 | cos.merge <= 0) {
    stop("cos.merge must be between 0 and 1")
  }
  if (sig.level < 0.0 | sig.level > 1.0) stop("sig.level must be between 0 and 1")

  set.seed(sampling_seed(chain), kind="Mersenne-Twister", normal.kind="Inversion")

  # number of categories, DPs, samples, data items at each DP, and frozen priors
  ncat <- numcateg(final_hdpState(chain))
  ndp <- numdp(final_hdpState(chain))
  nsamp <- hdp_settings(chain)$n
  numdata <- sapply(dp(final_hdpState(chain)), numdata)
  pseudo <- pseudoDP(final_hdpState(chain))
  is_prior <- length(pseudo) > 0
  if (is_prior) {
    priorcc <- 1:length(pseudo)
  }

  # Step (1)
  # Make each ccc (clust_categ_counts) and
  # cdc (clust_dp_counts) matrix have the
  # same number of columns
  maxclust <- max(numcluster(chain))
  clust_label <- 1:maxclust

  ccc_1 <- lapply(clust_categ_counts(chain), function(x){
    ans <- cbind(x, matrix(0, nrow=ncat, ncol=(maxclust-ncol(x)+1)))
    return(ans[, -ncol(ans)])
  })

  cdc_1 <- lapply(clust_dp_counts(chain), function(x){
    ans <- cbind(x, matrix(0, nrow=ndp, ncol=(maxclust-ncol(x)+1)))
    return(ans[, -ncol(ans)])
  })

  # Step (2)
  # Match up raw clusters (matrix columns) across posterior samples (columns not
  # guaranteed to keep same component through all samples)

  # if priors, hold back those columns and just match the new clusters
  if (is_prior){
    ccc_hold <- lapply(ccc_1, function(x) as.matrix(x[,priorcc]))
    ccc_1 <- lapply(ccc_1, function(x) as.matrix(x[,-priorcc]))
  }

  # K-centroids clustering of all raw clusters with cannot-link constraints
  # within each posterior sample, Manhattan distance and median centroid
  mclust <- ncol(ccc_1[[1]])

  if (mclust==1){
    ccc_label <- rep(1, length(ccc_1))

  } else{
    ccc_unlist <- t(do.call(cbind, ccc_1))
    groupfactor <- rep(1:nsamp, each=mclust)
    initial_clust <- rep(1:mclust, times=nsamp)

    ccc_clust <- flexclust::kcca(ccc_unlist, k=initial_clust,
                               group=groupfactor,
                               family=flexclust::kccaFamily("kmedians",
                                                groupFun="differentClusters"))

    # want this plot to be as simple as possible
    # tmp <- matrix(flexclust::clusters(ccc_clust), byrow=T, ncol=mclust)
    # matplot(tmp, type='l', lty=1, main="kmedians")

    ccc_label <- split(flexclust::clusters(ccc_clust), groupfactor)

    remove(ccc_unlist, groupfactor, initial_clust, ccc_clust)

    }


  if(is_prior){
    ccc_label <- lapply(ccc_label, `+`, max(priorcc))
    ccc_label <- lapply(ccc_label, function(x) c(priorcc, x))

    ccc_1 <- mapply(cbind, ccc_hold, ccc_1, SIMPLIFY = FALSE)

    remove(ccc_hold)
  }

  ccc_2 <- mapply(function(ccc, label) {
    colnames(ccc) <- label
    ccc[, order(as.numeric(colnames(ccc)))]
    },
                  ccc_1, ccc_label, SIMPLIFY=FALSE)

  cdc_2 <- mapply(function(cdc, label) {
    colnames(cdc) <- label
    cdc[, order(as.numeric(colnames(cdc)))]
    },
                  cdc_1, ccc_label, SIMPLIFY=FALSE)

  remove(ccc_1, cdc_1, ccc_label)

  # Step (3)
  # Merge the ccc columns with high cosine similarity.
  avgdistn <- matrix(0, nrow=ncat, ncol=maxclust)
  for (i in 1:maxclust){
    distns <- sapply(ccc_2, function(x) x[, i]/sum(x[, i]))
    avgdistn[, i] <- rowMeans(distns, na.rm=T)
  }
  clust_cos <- lsa::cosine(avgdistn)
  clust_same <- (clust_cos > cos.merge & lower.tri(clust_cos))
  same <- which(clust_same, arr.ind=TRUE) # merge these columns

  # if priors, do not merge two prior clusters
  if (is_prior){
    if (length(same)>0) {
      ignore <- which(apply(same, 1, function(x) all(x %in% priorcc)))
      same <- matrix(same[-ignore,], ncol=2)
    }
  }

  # update clust_label vector to reflect the merging of columns.
  if (length(same)>0){
    for (i in 1:nrow(same)){
      clust_label[same[i, 1]] <- clust_label[same[i, 2]]
    }
    remove(i)
  }
  ccc_3 <- lapply(ccc_2, merge_cols, clust_label)
  cdc_3 <- lapply(cdc_2, merge_cols, clust_label)
  clust_label <- colnames(ccc_3[[1]])
  if (any(clust_label != colnames(cdc_3))) stop("problem in step 3!")

  remove(avgdistn, distns, clust_cos, clust_same, same, ccc_2, cdc_2)


  # Step (4)
  # Assign components with no *significantly* non-zero data categories
  # to component '0'
  use_clust <- c()
  for (ii in 1:ncol(ccc_3[[1]])) {
    compii <- sapply(ccc_3, function(x) x[,ii])
    lowerb <- apply(compii, 1, function(y) {
                samp <- coda::as.mcmc(y)
                if (min(sum(!is.na(samp)), sum(!is.nan(samp))) %in% c(0,1)) {
                  NaN
                } else {
                  round(coda::HPDinterval(samp, prob = sig.level)[1], 3)
                }
              })
    if(any(lowerb>0)) use_clust <- c(use_clust, colnames(ccc_3[[1]])[ii])
  }

  # if priors, must include them all
  if(is_prior){
    use_clust <- union(use_clust, as.character(priorcc))
  }

  # update clust_label vector
  clust_label[which(!clust_label %in% use_clust)] <- 0
  ccc_4 <- lapply(ccc_3, merge_cols, clust_label)
  cdc_4 <- lapply(cdc_3, merge_cols, clust_label)
  clust_label <- colnames(ccc_4[[1]])
  if (any(clust_label != colnames(cdc_4))) stop("problem in step 4!")

  remove(compii, ccc_3, cdc_3, ii, lowerb, use_clust)

  # Step (5)
  # Assign components with no *significantly* non-zero sample exposure
  # to component '0' (disregarding DP nodes with no data items (parent nodes))
  use_clust <- c()
  disregard <- which(numdata==0)
  for (ii in 1:ncol(cdc_4[[1]])) {
    compii <- sapply(cdc_4, function(x) x[,ii])
    lowerb <- apply(compii[-disregard,], 1, function(y) {
      samp <- coda::as.mcmc(y)
      if (min(sum(!is.na(samp)), sum(!is.nan(samp))) %in% c(0,1)) {
        NaN
      } else {
        round(coda::HPDinterval(samp, prob = sig.level)[1], 3)
      }
    })
    if(any(lowerb>0)) use_clust <- c(use_clust, colnames(cdc_4[[1]])[ii])
  }

  # if prior sigs, must include them all
  if(is_prior){
    use_clust <- union(use_clust, as.character(priorcc))
  }

  # update clust_label vector
  clust_label[which(!clust_label %in% use_clust)] <- 0
  ccc_5 <- lapply(ccc_4, merge_cols, clust_label)
  cdc_5 <- lapply(cdc_4, merge_cols, clust_label)
  clust_label <- colnames(ccc_5[[1]])
  if (any(clust_label != colnames(cdc_5))) stop("problem in step 5!")

  remove(compii, ccc_4, cdc_4, ii, lowerb, use_clust)

  # Step (6)
  # Rename overall component, order by number of data items (on average)
  # 0th component still goes first

  avg_ndi <- rowMeans(sapply(ccc_5, colSums))
  colorder <- c(1, setdiff(order(avg_ndi, decreasing=T), 1))

  ccc_6 <- lapply(ccc_5, function(x) {
    if (is_prior) {
      colnames(x) <- c("0", paste0("P", priorcc),
                       paste0("N", 1:(ncol(x)-length(priorcc)-1)))[1:ncol(x)]
    }
    x <- x[, colorder]
    if (is_prior){
      update <- which(grepl("N", colnames(x)))
      if (length(update)>0){
        colnames(x)[update] <- paste0("N", 1:length(update))
      }
    } else {
      colnames(x) <- 0:(ncol(x)-1)
    }
    return(x)
  })

  cdc_6 <- lapply(cdc_5, function(x) {
    if (is_prior) {
      colnames(x) <- c("0", paste0("P", priorcc),
                       paste0("N", 1:(ncol(x)-length(priorcc)-1)))[1:ncol(x)]
    }
    x <- x[, colorder]
    if (is_prior){
      update <- which(grepl("N", colnames(x)))
      if (length(update)>0){
        colnames(x)[update] <- paste0("N", 1:length(update))
      }
    } else {
      colnames(x) <- 0:(ncol(x)-1)
    }
    return(x)
  })

  # number of components
  ncomp <- length(colorder)

  remove(ccc_5, cdc_5, avg_ndi, colorder)


  # Step (7)
  # Convert ccc into list of length ncomp, with matrices nsamp*ncat

  # if prior sigs, remove pseudo data from ccc
  if (is_prior){
    pseudodata <- sapply(dp(final_hdpState(chain))[pseudo],
                         function(x) table(factor(x@datass, levels=1:ncat)))
    colnames(pseudodata) <- paste0("P", 1:ncol(pseudodata))

    cn <- colnames(ccc_6[[1]])
    cnp <- which(grepl("P", cn))

    ccc_6 <- lapply(ccc_6, function(x) {
      x[,cnp] <- x[,cnp] - pseudodata[,cn[cnp]]
      return(x)
    })
  }

  ccc_ans <- rep(list(matrix(0, nrow=nsamp, ncol=ncat)), ncomp)
  for (i in 1:ncomp){
    ccc_ans[[i]] <- t(sapply(ccc_6, function(x) x[, i]))
  }
  names(ccc_ans) <- colnames(ccc_6[[1]])

  # Convert cdc into list of length ndp, with matrices nsamp*ncomp
  cdc_ans <- rep(list(matrix(0, nrow=nsamp, ncol=ncomp)), ndp)
  for (i in 1:ndp){
    cdc_ans[[i]] <- t(sapply(cdc_6, function(x) x[i, ]))
  }

  remove(ccc_6, cdc_6)

  # Step (8)
  # Calculate mean and 95% credibility interval for each component's
  # categorical data distribution
  ccc_norm <- lapply(ccc_ans, function(x) x/rowSums(x, na.rm=TRUE))

  ccc_mean <- t(sapply(ccc_norm, colMeans, na.rm=TRUE))

  ccc_credint <- lapply(ccc_norm, function(x) {
    apply(x, 2, function(y) {
      samp <- coda::as.mcmc(y)
      if (min(sum(!is.na(samp)), sum(!is.nan(samp))) %in% c(0,1)) {
          c(NaN, NaN)
      } else {
          round(coda::HPDinterval(samp, prob = sig.level), 3)
      }
    })
  })

  # Step (8)
  # Calculate mean and 95% credibility interval for each DP's
  # distribution over components (counts)
  cdc_norm <- lapply(cdc_ans, function(x) x/rowSums(x, na.rm=TRUE))

  cdc_mean <- t(sapply(cdc_norm, colMeans, na.rm=TRUE))

  cdc_credint <- lapply(cdc_norm, function(x) {
    apply(x, 2, function(y) {
      samp <- coda::as.mcmc(y)
      if (min(sum(!is.na(samp)), sum(!is.nan(samp))) %in% c(0,1)) {
        c(NaN, NaN)
      } else {
        round(coda::HPDinterval(samp, prob = sig.level), 3)
      }
    })
  })

  # add extracted components into chain hdpSampleChain slots
  chain@numcomp <- as.integer(ncomp - 1)

  # proportion of data explained by extracted components?
  avcount <- colMeans(sapply(ccc_ans, rowSums, na.rm=TRUE), na.rm=TRUE)
  chain@prop.ex <- round(1-avcount[1]/sum(avcount), 3)

  chain@comp_cos_merge <- cos.merge

  chain@comp_categ_counts <- ccc_ans
  chain@comp_dp_counts <- lapply(cdc_ans, as, "dgCMatrix")
  chain@comp_categ_distn <- list(mean=ccc_mean,
                                 cred.int=ccc_credint)
  chain@comp_dp_distn <- list(mean=cdc_mean,
                              cred.int=cdc_credint)

  # check validity and return
  if (!validObject(chain)) warning("Not a valid hdpSampleChain object.")
  return(chain)
}

#' Extracts components from multiple MCMC chains
#' @export
#' @importFrom "hdp" comp_categ_counts numcateg final_hdpState numdp hdp_settings dp pseudoDP comp_dp_counts
#' @importFrom "methods" as validObject
hdp_extract_comp_multi <- function(chains, cos.merge=0.90, sig.level=0.95, redo=TRUE){

  # input checks
  if (class(chains) != "hdpSampleMulti") {
    stop("chains must have class hdpSampleMulti")
  }
  if (!validObject(chains)) stop("chains not a valid hdpSampleMulti object")
  if (class(cos.merge) != "numeric" | cos.merge >=1 | cos.merge <= 0) {
    stop("cos.merge must be between 0 and 1")
  }
  if (class(redo) != "logical") stop("redo must be TRUE or FALSE")
  if (sig.level < 0.0 | sig.level > 1.0) stop("sig.level must be in the range (0, 1)")

  # list of hdpSampleChain objects
  chlist <- chains@chains
  nch <- length(chlist)

  # extract components within each chain
  for (i in 1:nch){
    if (!redo & length(comp_categ_counts(chlist[[i]])) > 0) next
    chlist[[i]] <- hdp_extract_comp_single(chlist[[i]], sig.level = sig.level, cos.merge=cos.merge)
  }
  chains@chains <- chlist

  # number of categories, DPs, samples, and data items at each DP
  ncat <- numcateg(final_hdpState(chlist[[1]]))
  ndp <- numdp(final_hdpState(chlist[[1]]))
  nsamp <- sapply(chlist, function(x) hdp_settings(x)$n)
  numdata <- sapply(dp(final_hdpState(chlist[[1]])), numdata)
  pseudo <- pseudoDP(final_hdpState(chlist[[1]]))
  is_prior <- length(pseudo) > 0
  numprior <- length(pseudo)

  remove(i)

  # Step (1) Match up components across chains
  # Get average comp_categ_counts matrix for each chain (average over posterior samples)
  ccc_0 <- lapply(chlist, function(x) apply(simplify2array(comp_categ_counts(x)), 2:3, mean))

  # Make each ccc matrix have the same number of columns
  maxcomp <- max(sapply(ccc_0, ncol))
  ccc_1 <- lapply(ccc_0, function(x){
    ans <- cbind(x, matrix(0, nrow=ncat, ncol=(maxcomp-ncol(x)+1)))
    return(ans[, -ncol(ans)])
  })

  remove(ccc_0)

  # if priors, they will be 1:numprior in order of name
  if (is_prior){
    pr_pos <- lapply(ccc_1, function(mat){
      sapply(paste0("P", 1:numprior), function(cn){
        which(colnames(mat)==cn)})
    })
    comp_mapping <- rep(list(rep(0, maxcomp)), nch)
    comp_mapping <- mapply(function(x, y){x[y] <- 1:length(y); return(x)},
                           comp_mapping, pr_pos, SIMPLIFY=FALSE)
  }

  # K-centroids clustering of components across chains
  # with cannot-link constraints within each chain.
  # Manhattan distance and median centroid.
  # exclude component 0 (and any priors) from each chain before matching
  if (is_prior){
    ccc_1_toclust <- mapply(function(mat, priors) {mat[,-c(1,priors)]},
                            ccc_1, pr_pos, SIMPLIFY=FALSE)
  } else {
    ccc_1_toclust <- lapply(ccc_1, function(mat) mat[,-1])
  }

  mclust <- ncol(ccc_1_toclust[[1]])

  if (mclust!=0) {
    ccc_unlist <- t(do.call(cbind, ccc_1_toclust))
    groupfactor <- rep(1:nch, each=mclust)
    initial_clust <- rep(1:mclust, times=nch)

    ccc_clust <- flexclust::kcca(ccc_unlist, k=initial_clust,
                                 group=groupfactor,
                                 family=flexclust::kccaFamily("kmedians",
                                                              groupFun="differentClusters"))

    # want this plot to be as simple as possible
    # tmp <- matrix(flexclust::clusters(ccc_clust), byrow=T, ncol=mclust)
    # matplot(tmp, type='l', lty=1, main="kmedians")

    clust_map <- split(flexclust::clusters(ccc_clust), groupfactor)

    if (is_prior){
      comp_mapping <- mapply(function(map, otherpos, vals){map[-c(1, otherpos)] <- vals+numprior; return(map)},
                             comp_mapping, pr_pos, clust_map, SIMPLIFY = FALSE)
    } else {
      comp_mapping <- lapply(clust_map, function(x) c(0, x))
    }

    remove(ccc_1, ccc_1_toclust, ccc_unlist, groupfactor, initial_clust,
           ccc_clust, clust_map, mclust)
  }


  # Step (2) Consolidate all ccc and cdc stats from all chains
  # Get comp_categ_counts list for each chain,
  # make the same length, and reorder elements so they match up to comp_mapping
  ccclist <- lapply(chlist, comp_categ_counts)

  cccmerge <- mapply(function(x, y, z){
    ans <- c(x, rep(list(matrix(0, nrow=y, ncol=ncat)), (maxcomp-length(x)+1)))
    ans <- ans[-length(ans)]
    names(ans) <- z
    ans <- ans[order(as.numeric(names(ans)))]
    return(ans)
  }, ccclist, nsamp, comp_mapping, SIMPLIFY = TRUE)

  # rbind matrices from different chains
  ccc_2 <- vector("list", nrow(cccmerge))
  names(ccc_2) <- rownames(cccmerge)
  for (i in 1:length(ccc_2)){
    ccc_2[[i]] <- Reduce(rbind, cccmerge[i, ])
  }

  # Get comp_dp_counts list for each chain,
  # make number of columns the same, and reorder columns to match comp_mapping
  cdclist <- lapply(chlist, comp_dp_counts)
  cdcmerge <- mapply(function(x, y, z) {
    x <- lapply(x, function(w){
            ans <- cbind(w, matrix(0, nrow=y, ncol=(maxcomp-ncol(w)+1)))
            ans <- ans[, -ncol(ans)]
            colnames(ans) <- z
            ans <- ans[,order(as.numeric(colnames(ans)))]
            return(ans)
          })
  }, cdclist, nsamp, comp_mapping, SIMPLIFY=FALSE)

  # rbind matrices from different chains
  cdc_2 <- vector("list", length(cdcmerge[[1]]))
  for (i in 1:length(cdc_2)) {
    cdc_2[[i]] <- Reduce(rbind, lapply(cdcmerge, `[[`, i))
    colnames(cdc_2[[i]]) <- 0:(maxcomp-1)
  }

  remove(i, cdclist, cdcmerge, cccmerge, ccclist, chlist, comp_mapping, maxcomp)

  comp_label <- names(ccc_2)
  if(any(comp_label != colnames(cdc_2[[1]]))) stop("problem in step 2!")

  # Step (3)
  # Merge components with high cosine similarity in ccc

  avgdistn <- sapply(ccc_2, function(x) colMeans(x)/sum(colMeans(x)))
  clust_cos <- lsa::cosine(avgdistn)
  clust_same <- (clust_cos > cos.merge & lower.tri(clust_cos))
  same <- which(clust_same, arr.ind=TRUE) # merge these

  # if priors, do not merge two prior clusters
  if (is_prior){
    if (length(same)>0) {
      ignore <- which(apply(same, 1, function(x) all(x %in% 1:numprior)))
      same <- matrix(same[-ignore,], ncol=2)
    }
  }

  # update comp_label vector to reflect the merging of columns.
  if (length(same)>0){
    for (i in 1:nrow(same)){
      comp_label[same[i, 1]] <- comp_label[same[i, 2]]
    }
    remove(i)
  }
  ccc_3 <- merge_elems(ccc_2, comp_label)
  cdc_3 <- lapply(cdc_2, merge_cols, comp_label)

  remove(ccc_2, cdc_2, avgdistn, clust_cos, clust_same, same)

  comp_label <- names(ccc_3)
  if (any(comp_label != colnames(cdc_3))) stop("problem in step 3!")

  # Step (4)
  # Assign components with no *significantly* non-zero data categories
  # to component '0'

  use_comp <- c()
  for (ii in 1:length(ccc_3)) {
    lowerb <- apply(ccc_3[[ii]], 2, function(y) {
      samp <- coda::as.mcmc(y)
      if (min(sum(!is.na(samp)), sum(!is.nan(samp))) %in% c(0,1)) {
        NaN
      } else {
        round(coda::HPDinterval(samp, prob = sig.level)[1], 3)
      }
    })
    if(any(lowerb>0)) use_comp <- c(use_comp, names(ccc_3)[ii])
  }

  # if priors, must include them all
  if(is_prior){
    use_comp <- union(use_comp, as.character(1:numprior))
  }

  # update comp_label vector
  comp_label[which(!comp_label %in% use_comp)] <- 0
  ccc_4 <- merge_elems(ccc_3, comp_label)
  cdc_4 <- lapply(cdc_3, merge_cols, comp_label)

  remove(ccc_3, cdc_3, ii, lowerb, use_comp)

  comp_label <- names(ccc_4)
  if (any(comp_label != colnames(cdc_4))) stop("problem in step 4!")

  # Step (5)
  # Assign components with no *significantly* non-zero sample exposure
  # to component '0' (disregarding DP nodes with no data items (parent nodes))

  use_comp <- c()
  disregard <- which(numdata==0)
  for (ii in 1:ncol(cdc_4[[1]])) {
    compii <- t(sapply(cdc_4, function(x) x[,ii]))
    lowerb <- apply(compii[-disregard,], 1, function(y) {
      samp <- coda::as.mcmc(y)
      if (min(sum(!is.na(samp)), sum(!is.nan(samp))) %in% c(0,1)) {
        NaN
      } else {
        round(coda::HPDinterval(samp, prob = sig.level)[1], 3)
      }
    })
    if(any(lowerb>0)) use_comp <- c(use_comp, colnames(cdc_4[[1]])[ii])
  }

  # if prior sigs, must include them all
  if(is_prior){
    use_comp <- union(use_comp, as.character(1:numprior))
  }

  # update comp_label vector
  comp_label[which(!comp_label %in% use_comp)] <- 0
  ccc_5 <- merge_elems(ccc_4, comp_label)
  cdc_5 <- lapply(cdc_4, merge_cols, comp_label)

  remove(compii, ccc_4, cdc_4, ii, lowerb, use_comp)

  comp_label <- names(ccc_5)
  if (any(comp_label != colnames(cdc_5))) stop("problem in step 5!")

  # Step (6)
  # Rename overall component, order by number of data items (on average)
  # 0th component still goes first
  avg_ndi <- sapply(ccc_5, function(x) mean(rowSums(x)))
  comporder <- c(1, setdiff(order(avg_ndi, decreasing=T), 1))

  # number of components
  ncomp <- length(comporder)

  if (is_prior) {
    compnames <- c("0", paste0("P", 1:numprior),
                  paste0("N", 1:(ncomp-numprior-1)))[1:ncomp]
    compnames <- compnames[comporder]
    update <- which(grepl("N", compnames))
    if (length(update)>0){
      compnames[update] <- paste0("N", 1:length(update))
    }
  } else {
    compnames <- 0:(ncomp-1)
  }

  ccc_ans <- ccc_5[comporder]
  names(ccc_ans) <- compnames

  cdc_ans <- lapply(cdc_5, function(x) {
    x <- x[, comporder]
    colnames(x) <- compnames
    return(x)
  })

  remove(ccc_5, cdc_5, avg_ndi, comporder, compnames)


  # Step (7)
  # Calculate mean and 95% credibility interval for each component's
  # categorical data distribution
  ccc_norm <- lapply(ccc_ans, function(x) x/rowSums(x, na.rm=TRUE))

  ccc_mean <- t(sapply(ccc_norm, colMeans, na.rm=TRUE))

  ccc_credint <- lapply(ccc_norm, function(x) {
    apply(x, 2, function(y) {
      samp <- coda::as.mcmc(y)
      if (min(sum(!is.na(samp)), sum(!is.nan(samp))) %in% c(0,1)) {
        c(NaN, NaN)
      } else {
        round(coda::HPDinterval(samp, prob = sig.level), 3)
      }
    })
  })

  # Calculate mean and 95% credibility interval for each DP's
  # distribution over components (counts)
  cdc_norm <- lapply(cdc_ans, function(x) x/rowSums(x, na.rm=TRUE))

  cdc_mean <- t(sapply(cdc_norm, colMeans, na.rm=TRUE))

  cdc_credint <- lapply(cdc_norm, function(x) {
    apply(x, 2, function(y) {
      samp <- coda::as.mcmc(y)
      if (min(sum(!is.na(samp)), sum(!is.nan(samp))) %in% c(0,1)) {
        c(NaN, NaN)
      } else {
        round(coda::HPDinterval(samp, prob = sig.level), 3)
      }
    })
  })

  # add extracted components into hdpSampleMulti slots
  chains@numcomp <- as.integer(ncomp - 1)

  # proportion of data explained by extracted components?
  avcount <- colMeans(sapply(ccc_ans, rowSums, na.rm=TRUE), na.rm=TRUE)
  chains@prop.ex <- round(1-avcount[1]/sum(avcount), 3)

  chains@comp_cos_merge <- cos.merge

  chains@comp_categ_counts <- ccc_ans
  chains@comp_dp_counts <- lapply(cdc_ans, as, "dgCMatrix")
  chains@comp_categ_distn <- list(mean=ccc_mean,
                                 cred.int=ccc_credint)
  chains@comp_dp_distn <- list(mean=cdc_mean,
                              cred.int=cdc_credint)

# check validity and return
  if (!validObject(chains)) warning("Not a valid hdpSampleMulti object.")
  return(chains)
}
