# ======== de.helpers =========
# functions from de.helpers -- consider installing the package to the project
# later
calc_perf_metrics <- function(p, tn, model, presorted = FALSE, nfdr = NULL, fdr_mtd = c("bfdr", "fdr")) {
  if (!presorted) {
    ord <- order(p)
    p <- p[ord]
    tn <- tn[ord]
  }
  fpr <- sapply(p, calc_rate, p = p, s = tn)
  tpr <- sapply(p, calc_rate, p = p, s = 1 - tn)
  tfdr <- sapply(p, calc_fdr, p = p, tn = tn)
  
  if (is.null(nfdr)) {
    fdr_mtd <- match.arg(fdr_mtd)
    fdr_fn <- switch (fdr_mtd,
                      bfdr = calc_bfdr,
                      fdr = function(p) p.adjust(p, method = "fdr")
    )
    nfdr <- fdr_fn(p)
  }
  tibble::tibble(model = model, p = p, tn = tn, fpr = fpr, tpr = tpr, tfdr = tfdr, nfdr = nfdr)
}
calc_fdr <- function(p0, p, tn) {
  sum(tn[p <= p0]) / sum(p <= p0)
}
calc_rate <- function(p0, p, s) {
  x <- sum(s[p <= p0])
  x / sum(s)
}
#' Calculate the estimated false discovery rate
#' @param p A vector of posterior probabilities
#' @export
calc_bfdr <- function(p, presorted = FALSE) {
  if (!presorted) p <- p[order(p)]
  sapply(p, function(p0) {
    sum(p[p <= p0]) / sum(p <= p0)
  })
}