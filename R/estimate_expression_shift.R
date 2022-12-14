#' @importFrom magrittr "%>%" "%<>%"
NULL

## calculate NMI for all cell pairs: it takes around 30s for 4960 cell pairs
nmi <- function(x,y){
  pxy <- cbind(x, y) %>% as.data.frame() %>% table %>% entropy::freqs.empirical()
  px <- rowSums(pxy)
  py <- colSums(pxy)
  Hx <- entropy::entropy.empirical(px, unit = "log2")
  Hy <- entropy::entropy.empirical(py, unit = "log2")
  mi <- entropy::mi.empirical(pxy)
  return(mi/sqrt(Hx*Hy))
}

## zip two vector into list
## examples: zip(letters, LETTERS)
zip <- function(a, b) {
  lapply(seq_along(a), function(i) c(a[i], b[i]))
}

#' @title Calculating median NMI for a given cells.
#' @param X matrix, discrete gene expression matrix. Rows, genes; columns, cells.
#' @param cells.1 vector, selected cells for group1. Match to colnames(X)
#' @param cells.2 vector, selected cells for group2. Match to colnames(X)
#' @param N.cells int, sample size of selected cells.
#' @param seed int, seed for sampling.
#' @param sample.on.cell.pairs bool, whether sample on cell pairs for NMI calculation.
#' @param N.cell.pairs int, sample size of cell pairs for NMI calculation.
#' @return numeric, median NMI
sample.nmi <- function(X, cells.1, cells.2=NULL, N.cells=100, seed=1, sample.on.cell.pairs=TRUE, N.cell.pairs=100) {
  set.seed(seed)
  if (is.null(cells.2)) {
    N.cells <- ifelse(length(cells.1) >= N.cells, N.cells, length(cells.1))
    cell.pairs <- utils::combn(
      x = sample(x = cells.1, size = N.cells, replace = FALSE),
      m = 2,
      simplify = F)
  } else {
    a <- sample(cells.1, size = (N.cells/2)^2, replace = T)
    b <- sample(cells.2, size = (N.cells/2)^2, replace = T)
    cell.pairs <- zip(a, b)
  }
  if(sample.on.cell.pairs) {
    N.cell.pairs <- ifelse(length(cell.pairs) >= N.cell.pairs, N.cell.pairs, length(cell.pairs))
    cell.pairs <- sample(cell.pairs, N.cell.pairs, replace = FALSE)
  }
  sapply(cell.pairs, function(xx) {
    x <- X[, xx[1]]
    y <- X[, xx[2]]
    nmi(x, y)
  }) %>% stats::median()
}

#' Overall expression shift of the provided gene expression matrix.
#' @param counts gene expression matrix
#' @param cell.types A named vector contains cell types for each cell in 'counts'
#' @param groups A named vector contains groups for each cell in 'counts'. Only 2 levels are surpported now.
#' @param n.bins number of bins for calculation NMI. Default: 20
#' @param min.cells minimal number of cells for each group. Default: 10
#' @param n.samples number of samplings for calculating expression shift of each cell type. Default: 100
#' @param n.cells number of cells sampled for calculation NMI in each group for each cell type. Default: min.cells
#' @param n.cells.pairs number of cell pairs sampled for calculation NMI on 'n.cells' sampled cells. Default: C(10,2)=45
#' @param cores the number of threads used. Default: 1
#'
#' @return A data.frame contains expression change between groups
#' @export
EstimateExpressionShift <- function(counts, cell.types, groups, n.bins=20, min.cells=10, n.samples=100, n.cells = min.cells, n.cells.pairs = 45, cores=1){
  ## check params
  if (missing(counts) || missing(cell.types) || missing(groups)) {
    stop("Must provide 'counts', 'cell.types', and 'groups'.")
  }
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("No rownames and (or) colnames in 'counts'.")
  }
  if (is.null(names(cell.types)) || is.null(names(groups))) {
    stop("No names in 'cell.types' and(or) 'groups'")
  }
  if (!all(colnames(counts) == names(cell.types))) {
    stop("The cell names between 'counts' (columns) and 'cell.types' are not matched.")
  }
  if (!all(colnames(counts) == names(groups))) {
    stop("The cell names between 'counts' (columns) and 'groups' are not matched.")
  }
  cores <- ifelse(cores > 0, cores, parallel::detectCores())
  ## divided profile into bins
  cat("\nstep1. Discretizing the gene expression matrix ...\n")
  disc.expr <- apply(counts, 2, function(xx) {
    infotheo::discretize(xx, disc = "equalwidth", n.bins) %>% unlist()
  })
  rownames(disc.expr) <- rownames(counts)
  ## get pairs of cells for comparation
  cell.counts <- table(cell.types, groups) %>% as.data.frame()
  cell.counts <- cell.counts[cell.counts$Freq >= min.cells, ]
  # cell.counts <- subset(cell.counts, Freq >= min.cells)
  cell.type.group.pairs <- table(cell.counts$cell.types)

  ## only 2 levels in groups are surpported now. [2022-05-23]
  cell.type.used <- names(cell.type.group.pairs)[cell.type.group.pairs == 2]
  if (is.null(cell.type.used)) {
    stop("No cell type were selected. Try a smaller 'min.cells'.")
  }
  cat("\nstep2. Calculating the global expression shift ...\n")
  ## get expression change
  group.levels <- unique(groups)
  ## es.df: rows = samples, columns = cell types, values = es.index
  es.df <- do.call(rbind, parallel::mclapply(1:n.samples, function(xx) {
    sapply(cell.type.used, function(ct) {
      ## within groups
      cells.1 <- names(groups)[cell.types == ct & groups == group.levels[1]]
      cells.2 <- names(groups)[cell.types == ct & groups == group.levels[2]]
      a <- sample.nmi(X = disc.expr, cells.1 = cells.1, N.cells = n.cells,
                      seed = xx, sample.on.cell.pairs = T, N.cell.pairs = n.cells.pairs)
      b <- sample.nmi(X = disc.expr, cells.1 = cells.2, N.cells = n.cells,
                      seed = xx, sample.on.cell.pairs = T, N.cell.pairs = n.cells.pairs)
      ## between groups
      c <- sample.nmi(X = disc.expr, cells.1 = cells.1, cells.2 = cells.2, N.cells = n.cells,
                      seed = xx, sample.on.cell.pairs = T, N.cell.pairs = n.cells.pairs)
      es.index <- (a+b) / (2*c) # a,b = NMI within groups; c = NMI between groups.
      return(es.index)
    })
  }, mc.cores = cores)
  )
  es.df %<>% as.data.frame()
  es.df <- es.df %>% dplyr::mutate(rep = 1:nrow(es.df), .before=1)
  # es.df <- es.df %>%
  #   as.data.frame() %>%
  #   dplyr::mutate(rep = 1:nrow(.), .before=1)
  return(es.df)
}
