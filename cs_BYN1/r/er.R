library(edgeR)
set.seed(5)
sim_list <- readRDS("cs_BYN1/sim_list.rds")
G <- sim_list$G
# edgeR::edgeRUsersGuide()
dgel <- DGEList(counts = sim_list$y_g, group = sim_list$x_g)
colnames(dgel) <- paste0(rep(1:5, 2), "_", rep(c("ctl", "trt"), each = 5))

# normalization
dgel <- normLibSizes(dgel)

# design
trt <- factor(sim_list$x_g)
design <- model.matrix(~trt)
rownames(design) <- colnames(dgel)

# dispersion
dgel <- estimateDisp(dgel, design = design)

# glm
fit <- glmFit(dgel, design)

# tests 
lrt <- glmLRT(fit, coef = 2)
tt <- topTags(lrt, n = G, sort.by = "none")

saveRDS(list(fit = fit, lrt = lrt, tt = tt), file = "cs_BYN1/erfit.rds")

# extras
# plotMDS(dgel)
# dgel$samples
