# get data set
data(kidney, package = 'SimSeq')
# extract data set parameters
.get_params <- function(counts, group) {
  dge <- edgeR::DGEList(counts = counts, group = group)
  dge <- edgeR::calcNormFactors(dge)
  dge <- edgeR::estimateCommonDisp(dge)
  dge <- edgeR::estimateTagwiseDisp(dge)
  disp <- dge$tagwise.dispersion # Dispersion
  mean <- apply(counts,1,mean)
  list(disp = disp, mean = mean)
}
#' counts: a matrix of counts
#' groups: a vector giving treatment groups
get_params <- function(counts, groups, mean_count_filter = 10) {
  mean_total = apply(counts,1,mean)
  index_filter = which(mean_total > mean_count_filter)
  counts <- counts[index_filter, ]
  params <- list()
  for (group_name in unique(groups)) {
    params[[group_name]] <- .get_params(counts = counts[, groups == group_name], 
                                        group = factor(groups[groups == group_name]))
  }
  params[["total"]] <- .get_params(counts = counts, group = factor(groups))
  params
}


# simulate a data set
sim_data_set <- function(params, G, N, p_diff_exp = 0.2, p_up = 0.5, basefc = 1.5, 
                         sample_effects = TRUE) {
  idx <- sample(1:length(params$total$mean), size = G)
  mu1 <- params$total$mean[idx]
  mu2 <- mu1
  n_diff_exp <- round(G * p_diff_exp)
  upidx <- 1:round(n_diff_exp * p_up)
  dnidx <- (max(upidx) + 1):n_diff_exp
  fc <- basefc + rexp(n = n_diff_exp, rate = 1)
  mu2[upidx] <- mu2[upidx]*fc[upidx]
  mu2[dnidx] <- mu2[dnidx]/fc[dnidx]
  disp <- params$total$disp
  phi1 <- disp[idx]
  phi2 <- phi1
  
  counts <- matrix(nrow = G, ncol = N)
  if (sample_effects == TRUE) {
    S1 <- runif(N / 2, min = 0.7, max = 1.3)
    S2 <- runif(N / 2, min = 0.7, max = 1.3)
  } else {
    S1 <- S2 <- rep(1, N / 2)
  }
  
  for (i in 1:G) {
    counts[i, 1:(N / 2)] <- sapply(S1, FUN = function(x) rnbinom(1, 1 / phi1[i], mu = mu1[i] * x))
    counts[i, (N / 2 + 1):N] <- sapply(S2, FUN = function(x) rnbinom(1, 1 / phi2[i], mu = mu2[i] * x))
  }
  list(counts = counts, mu1 = mu1, mu2 = mu2, phi1 = phi1, phi2 = phi2, S1 = S1, S2 = S2, fc = fc)
}

kparams <- get_params(counts = kidney$counts, groups = kidney$treatment, mean_count_filter = 10)
ds <- sim_data_set(params = kparams, G = 1000, N = 10, basefc = 1.5)

sim_list <- list(
  
)