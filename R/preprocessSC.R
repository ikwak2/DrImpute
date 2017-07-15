#' A function for preprocessing gene expression matrix.
#'
#' Preprocess gene expression data 
#' 
#' @param X Gene expression matrix (Gene by Cell). 
#'
#' @param min.expressed.genes Cell level filtering criteria. For a given cell, if the number of expressed genes are less than min.expressed.gene, we filter it out.  
#'
#' @param min.expressed.cell Gene level filtering criteria. For a given gene, if the number of expressed cells are less than min.expressed.cell, we filter it out.  
#'
#' @param max.expressed.ratio Gene level filtering criteria. For a given gene, if the ratio of expressed cells are larger than max.expressed.ratio, we filter it out.
#' 
#' @param normalize.by.size.effect Normaize using size factor.
#'
#' @return Filtered gene expression matrix
#'
#' @author Wuming Gong
#'
#' @references
#' Il-Youp Kwak, Wuming Gong, Kaoko Koyano-Nakagawa and Daniel J. Garry (2017+)
#' DrImpute: Imputing dropout eveents in single cell RNA sequencing data
#' 
#' @examples
#'
#' data(usoskin_ex)
#' exdata <- preprocessSC(exdata)
#'
#' @seealso \code{\link{SCimpute}}

preprocessSC <- function(X, min.expressed.gene = 0, min.expressed.cell = 2, max.expressed.ratio = 1, normalize.by.size.effect = FALSE){

	M0 <- ncol(X)
	N0 <- nrow(X)

	cat('----------------------------------------------------------------\n')
	cat('Preprocess single cell RNA-seq expression matrix\n')
	cat('----------------------------------------------------------------\n')
	cat(sprintf('number of input genes(nrow(X))=%.d\n', N0))
	cat(sprintf('number of input cells(ncol(X))=%.d\n', M0))
	X <- X[, colSums(X > 1) >= min.expressed.gene, drop = FALSE]
	M <- ncol(X)
	cat(sprintf('number of input cells that express at least %d genes=%.d\n', min.expressed.gene, M))
	X <- X[rowSums(X > 1) <= max.expressed.ratio * M & rowSums(X > 1) >= min.expressed.cell, , drop = FALSE]
	N <- nrow(X)
	cat(sprintf('number of input genes that are expressed in at least %d cells and at most %.0f%% cells=%.d\n', min.expressed.cell, max.expressed.ratio * 100, N))
	cat(sprintf('sparsity of expression matrix=%.1f%%\n', 100 * (N * M - sum(X > 0)) / (N * M)))

	if (normalize.by.size.effect){
		cat(sprintf('scaling raw read counts by size factor\n'))
	  sf <- apply((X + 1) / exp(rowMeans(log(X + 1))), 2, median)
		X <- t(t(X) / sf)
	}

	as.matrix(X)

} # end of preprocess
