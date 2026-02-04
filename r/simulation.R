# get data set
data(kidney, package = 'SimSeq')
# extract data set parameters
counts=kidney$counts # RNA-seq read count data
index.cancer=(1:72)*2 # cancer sample index
index.normal=index.cancer-1 # normal sample index
counts= counts[,c(index.cancer, index.normal)] # Arrange samples for convinience
mean.total = apply(counts,1,mean)
index.filter = which(mean.total < 10)
counts <- counts[!index.filter, ]
.get_params <- function(counts, group) {
  dge=DGEList(counts=counts, group = group)
  dge=calcNormFactors(dge)
  dge=estimateCommonDisp(dge)
  dge=estimateTagwiseDisp(dge)
  disp = dge$tagwise.dispersion # Dispersion
  mean = apply(counts,1,mean)
  list(disp = disp, mean = mean)
}
# Get count mean and dispersion using edgeR package

# Mean and dispersion values are obtained separately from cancer and normal samples when different dispersion is assummed between two sample types.
# Mean and dispersion values from normal samples
params.normal <- .get_params(counts = counts[, 73:144], group = factor(rep(2, 72)))
# Mean and dispersion values from cancer samples
params.cancer <- .get_params(counts = counts[, 1:72], group = factor(rep(1, 72)))

# Gene filtering: Genes having small read count (<10) are filtered
#### edgeR would recommend doing this filtering BEFORE fitting the models - 
#### also much nicer code to do it that way???
mean.total = apply(counts,1,mean)
index.filter = which(mean.total < 10)
mean.total = mean.total[-index.filter]
disp.normal = disp.normal[-index.filter]
disp.cancer = disp.cancer[-index.filter]
mean.normal = mean.normal[-index.filter]
mean.cancer = mean.cancer[-index.filter]


# Mean and dispersion values obtained using all samples when same dispersion is assummed between two sample types..
dge.total = DGEList(counts = counts, group = factor(c(rep(1,72),rep(2,72))))
dge.total = calcNormFactors(dge.total)
dge.total = estimateCommonDisp(dge.total)
dge.total = estimateTagwiseDisp(dge.total)
disp.total = dge.total$tagwise.dispersion
disp.total = disp.total[-index.filter]