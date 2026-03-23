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

kparams <- get_params(counts = kidney$counts, groups = kidney$treatment, mean_count_filter = 10)
saveRDS(kparams, "w_wout_test1/data/kparams.rds")