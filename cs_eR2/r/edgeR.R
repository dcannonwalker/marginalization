library(edgeR)
set.seed(5)
sim_list <- readRDS("cs_eR2/data/sim_list.rds")
sim_list$y_g
# edgeR::edgeRUsersGuide()
dgel <- DGEList(counts = sim_list$y_g, group = sim_list$x_g)
colnames(dgel) <- paste0(rep(1:5, 2), "_", rep(c("ctl", "trt"), each = 5))

# normalization
dgel <- normLibSizes(dgel)

# design
sim_list$z_g
pair <- factor(rep(1:5, 2))
trt <- factor(sim_list$x_g)
design <- model.matrix(~trt + pair)
rownames(design) <- colnames(dgel)

# dispersion
dgel <- estimateDisp(dgel, design = design)

# glm
fit <- glmFit(dgel, design)

# tests 
lrt <- glmLRT(fit, coef = 2)
tt <- topTags(lrt, n = 50, sort.by = "none")

saveRDS(list(fit = fit, lrt = lrt, tt = tt), file = "cs_eR2/data/eR_fit.rds")

# extras
# plotMDS(dgel)
# dgel$samples
