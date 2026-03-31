library(edgeR)
set.seed(5)
sim_list <- readRDS("w_wout_test2/data/sim_list.rds")
G <- sim_list$G
# edgeR::edgeRUsersGuide()
dgel <- DGEList(counts = sim_list$y_g, group = sim_list$x_g)
colnames(dgel) <- paste0(rep(1:(sim_list$N_g / 2), 2), "_", rep(c("ctl", "trt"), each = sim_list$N_g / 2))

# normalization
dgel <- normLibSizes(dgel)

# design
trt <- factor(sim_list$x_g)
pairs <- factor(rep(1:(sim_list$N_g / 2), 2))
design <- model.matrix(~trt + pairs)
rownames(design) <- colnames(dgel)

# dispersion
dgel <- estimateDisp(dgel, design = design)

# glm
fit <- glmFit(dgel, design)

# tests 
lrt <- glmLRT(fit, coef = 2)
tt <- topTags(lrt, n = G, sort.by = "none")

saveRDS(list(fit = fit, lrt = lrt, tt = tt), file = "w_wout_test2/data/erfit.rds")

lrt$table
hist(lrt$unshrunk.coefficients[, 2], breaks = 200)
hist(lrt$table$logFC, breaks = 200)
lrt$coefficients
