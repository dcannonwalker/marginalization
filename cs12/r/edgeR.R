library(edgeR)
set.seed(5)
sim_list <- readRDS("cs9/data/sim_list.rds")
sim_list$y_g
?edgeR
?edgeR::`edgeR-package`
edgeR::edgeRUsersGuide()
dgel <- DGEList(counts = sim_list$y_g)
dgel <- normLibSizes(dgel)
dgel$samples
