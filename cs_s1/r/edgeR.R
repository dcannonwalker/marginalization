library(edgeR)
set.seed(5)
sim_list <- readRDS("cs_s1/data/sim_list.rds")
G <- sim_list$G
N_g <- sim_list$N_g
# edgeR::edgeRUsersGuide()
dgel <- DGEList(counts = sim_list$y_g, group = sim_list$x_g)
colnames(dgel) <- paste0(rep(1:(N_g / 2), 2), "_", rep(c("ctl", "trt"), each = (N_g / 2)))

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

saveRDS(list(fit = fit, lrt = lrt, tt = tt), file = "cs_s1/data/eR_fit.rds")

# extras
# plotMDS(dgel)
# dgel$samples
