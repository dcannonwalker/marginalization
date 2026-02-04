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
  index_filter = which(mean_total < mean_count_filter)
  counts <- counts[index_filter, ]
  params <- list()
  for (group_name in unique(groups)) {
    params[[group_name]] <- .get_params(counts = counts[, groups == group_name], 
                                        group = factor(groups[groups == group_name]))
  }
  params[["total"]] <- .get_params(counts = counts, group = factor(groups))
  params
}

kparams <- get_params(counts = kidney$counts, groups = kidney$treatment)

# simulate a data set
n.var <- 100
random.index = sample(1:length(kparams$total$mean), size=n.var)
sample.mean1 = kparams$total$mean[random.index]
sample.mean2 = sample.mean1
n.diffexp <- n.var / 2
fraction.upregulated <- 0.5
upindex = 1:round(n.diffexp*fraction.upregulated)
dnindex = round(n.diffexp*fraction.upregulated+1):n.diffexp
basefc <- 1.5
fc <- basefc + rexp(n = n.diffexp, rate = 1)
sample.mean2[upindex] = sample.mean2[upindex]*fc[upindex]
sample.mean2[dnindex] = sample.mean2[dnindex]/fc[dnindex]
disp.total=kparams$total$disp
sample.disp1 = disp.total[random.index]
sample.disp2 = sample.disp1

counts = matrix(nrow=n.var, ncol = 2*s)
if(random_sampling==TRUE){
  rand1=runif(s,min=0.7,max=1.3)
  rand2=runif(s,min=0.7,max=1.3)
}else{
  rand1=rep(1,s)
  rand2=rep(1,s)
}
for(i in 1:n.var)
{
  counts[i,1:s] = sapply(rand1, FUN = function(x) rnbinom(1, 1/sample.disp2[i], mu=sample.mean2[i]*x))
  counts[i,(s+1):(2*s)] = sapply(rand2, FUN = function(x) rnbinom(1, 1/sample.disp1[i], mu=sample.mean1[i]*x))
}
