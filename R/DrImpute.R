#' Imputing dropout events in single-cell RNA-sequencing data.
#'
#' Imputing dropout events in single-cell RNA-sequencing data.
#'
#' @param X Gene expression matrix (gene by cell). 
#'
#' @param ks Number of cell clustering groups. Default set to ks = 10:15. 
#'
#' @param dists Distribution matrices to use. Default is set to c("spearman", "pearson"). "eucleadian" can be added as well. 
#'
#' @param method Use "mean" for mean imputation, "med" for median imputation.
#'
#' @param cls User can manually provide clustering information. Using different base clusterings. each row represent different clusterings. each column represent each cell.
#'
#' @param seed User can provide a seed.
#'
#' @param zerop zero percentage of resulting imputation is at least zerop. 
#'
#' @return Imputed Gene expression matrix (gene by cell).
#'
#' @author Il-Youp Kwak
#'
#' @references
#' Il-Youp Kwak, Wuming Gong, Kaoko Koyano-Nakagawa and Daniel J. Garry (2017+)
#' DrImpute: Imputing dropout eveents in single cell RNA sequencing data
#'
#' @examples
#'
#' data(usoskin_ex)
#' exdata <- preprocessSC(exdata)
#' logdat <- log(exdata+1)
#' cls <- getCls(logdat)
#' 
#' logdat_imp <- DrImpute(logdat, cls = cls)
#' 

DrImpute <- function(X, ks=10:15, dists=c("spearman", "pearson"), method = "mean", cls=NULL, seed = 1, zerop = 0) {

    set.seed(seed)
    if(is.null(cls))
        cls <- getCls(X=X, ks = ks, dists = dists)

    nc <- dim(X)[2]
    ng <- dim(X)[1]
    Tot <- nc * ng
    czp <- sum(X==0) / Tot

    if( method == "mean") {
        nX <- imp0clC(X, cls)
        imps <- nX[X == 0 & nX != 0]
        nimp =  max(0,sum(X==0) - Tot * zerop)

        if( nimp > 0 ) {
            threshold = sort(imps, decreasing=T)[min(nimp, length(imps))]
            nX[ X == 0 & nX < threshold ] = 0
        }

    } else if (method == "median") {
        nX <- imp0clC2(X, cls)
    } 

    rownames(nX) <- rownames(X)
    colnames(nX) <- colnames(nX)

    cat("\n Zero percentage : \n")
    cat(sprintf("Before impute : %.f percent. \n", sum(X==0) / Tot *100) )
    cat(sprintf("After impute : %.f percent. \n", sum(nX==0) / Tot *100) )
    cat(sprintf("%.f percent of zeros are imputed. \n", sum(X==0) / Tot *100 - sum(nX==0) / Tot *100) )
    nX
}

