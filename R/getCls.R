#' get base clustering results using SC3 based clustering methods.
#'
#' Similarity matrix constructed using "pearson", "spearman" or "euclidean". K-means clustering is performed on first few number of principal components of similarity matrix.
#' 
#' @param X Log transformed gene expression matrix (Gene by Cell). 
#'
#' @param ks Number of cell clustering groups. Default set to ks = 10:15. 
#'
#' @param dists Distribution matrices to use. Default is set to c("spearman", "pearson"). "euclidean" can be added as well. 
#'
#' @param dim.reduc.prop Proportion of principal components to use for K-means clustering.
#'
#' @return A matrix object, Each row represent different clustering results.
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
#' @seealso  \code{\link{DrImpute}}  \code{\link{preprocessSC}}

getCls <- function(X, ks=10:15, dists=c("spearman", "pearson"), dim.reduc.prop = 0.05) {

    N = dim(X)[2]
    ddim = round(N * dim.reduc.prop)

    Ds <- NULL
    i = 1;
    for( d in dists) {
        if(d == "spearman") {
            cat("Calculating Spearman distance. \n")
            Ds[[i]] <- as.matrix(1 - cor(X, method = "spearman"))
            i = i + 1;
        } else if (d == "pearson") {
            cat("Calculating Pearson distance. \n")
            Ds[[i]] <- as.matrix(1 - cor(X, method = "pearson"))
            i = i + 1;
        } else if (d == "euclidean") {
            cat("Calculating Euclidean distance. \n")
            Ds[[i]] <- as.matrix(dist(t(X), method = "euclidean"))
            i = i + 1;
        }
    }

    cls <- NULL
    for(k in ks) {
        cat(sprintf(" Clustering for k : %d\n", k) )
        for(d in 1:length(Ds)) {
            cls <- rbind(cls, kmeans(prcomp(Ds[[d]], center = TRUE, scale. = TRUE)$rotation[,1:ddim], centers=k, iter.max = 1e+09, nstart = 1000)$cluster )
        }
    }
    cat(sprintf("cls object have %d number of clustering sets.\n\n", length(ks) * length(Ds) ))
    cls
}
