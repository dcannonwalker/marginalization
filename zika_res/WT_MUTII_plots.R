library(DESeq2)
library(ggplot2)
fboth <- readRDS("zika_res/deseq_WT_MUTII_RNA_RPF.rds")
frna <- readRDS("zika_res/deseq_WT_MUTII_RNA.rds")
frpf <- readRDS("zika_res/deseq_WT_MUTII_RPF.rds")

xx <- mcols(fboth)
yy <- mcols(frna)
zz <- mcols(frpf)
disp <- tibble(group_id = rownames(yy), disp_rna = yy$dispersion, disp_rpf = zz$dispersion)
disp_plot <- ggplot(disp, aes(disp_rna, disp_rpf)) + 
  geom_point() +
  theme_minimal() + 
  geom_abline(intercept = 0, slope = 1, color = "red") + 
  ggtitle("Dispersion estimates for RNA samples vs. RPF samples", subtitle = "WT vs. MUTII") +
  xlab("RNA dispersion") +
  ylab("RPF dispersion")
saveRDS(disp_plot, "zika_res/WT_MUTII_disp_plot.rds")

# trt effect in RNA
type_lfc_rna <- tibble(group_id = rownames(xx), lfc_rna_both = xx$type_WT_vs_MUTII, lfc_rna_rna = yy$type_WT_vs_MUTII)
type_lfc_rna_plot <- ggplot(type_lfc_rna, aes(lfc_rna_both, lfc_rna_rna)) + 
  geom_point() + 
  theme_minimal() + 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  ggtitle("Estimated treatment logFC in RNA samples", subtitle = "WT vs. MUTII") +
  xlab("estimate from full model") +
  ylab("estimate from RNA only model")
saveRDS(type_lfc_rna_plot, "zika_res/WT_MUTII_type_lfc_rna_plot.rds")
# trt effect in RPF
type_lfc_rpf <- tibble(group_id = rownames(xx), lfc_rpf_both = xx$type_WT_vs_MUTII + xx$typeWT.protocolRPF, lfc_rpf_rpf = zz$type_WT_vs_MUTII)
type_lfc_rpf_plot <- ggplot(type_lfc_rpf, aes(lfc_rpf_both, lfc_rpf_rpf)) + 
  geom_point() + 
  theme_minimal() + 
  geom_abline(slope = 1, intercept = 0, color = "red") +
  ggtitle("Estimated treatment logFC in RPF samples", subtitle = "WT vs. MUTII") +
  xlab("estimate from full model") +
  ylab("estimate from RPF only model")
saveRDS(type_lfc_rpf_plot, "zika_res/WT_MUTII_type_lfc_rpf_plot.rds")
