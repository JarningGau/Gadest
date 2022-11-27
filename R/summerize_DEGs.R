#' Calculating connection specificity index (CSI) matrix.
#' @param pccMat pearson correlation coefficient (PCC) matrix of genes.
#' @param cores the number of threads used. Default: 1
#' @return a CSI matrix
#' @export
CalculateCSIMatrix <- function(pccMat, cores=1) {
  CSI <- function(r1, r2) {
    delta <- pccMat[r1,r2]
    r.others <- setdiff(colnames(pccMat), c(r1,r2))
    N <- sum(pccMat[r1, r.others] < delta) + sum(pccMat[r2, r.others] < delta)
    M <- length(r.others) * 2
    return(N/M)
  }
  csiMat <- parallel::mclapply(rownames(pccMat), function(i) sapply(colnames(pccMat), function(j) CSI(i, j)), mc.cores = cores)
  csiMat <- do.call(rbind, csiMat)
  rownames(csiMat) <- rownames(pccMat)
  return(csiMat)
}

#' Binarizing the connection specificity index (CSI) matrix.
#' @param csiMat a CSI matrix
#' @param cutoff threshold for binarization. Default: 0.8
#' @return a binarized CSI matrix
#' @export
BinarizeCSIMatrix <- function(csiMat, cutoff=0.8) {
  csiMat.binary <- matrix(as.numeric(csiMat >= cutoff), nrow = nrow(csiMat))
  colnames(csiMat.binary) <- colnames(csiMat)
  rownames(csiMat.binary) <- rownames(csiMat)
  return(csiMat.binary)
}
