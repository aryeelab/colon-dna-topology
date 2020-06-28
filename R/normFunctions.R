
copyNumberNormalize <- function( object, cnvData, returnNormMat=FALSE ){
    over <- findOverlaps( rowRanges(object), cnvData )
    cnvsDf <- as.data.frame( mcols( cnvData ) )
    cnvCor <- split( subjectHits( over ), rownames(object)[queryHits( over )] )
    cnvCor <- lapply( cnvCor, function(x){
        colMeans(cnvsDf[x,])
    })
    cnvCor <- do.call(rbind, cnvCor)
    cnvCor <- cnvCor[,match(colnames(object), colnames(cnvCor))]
    colnames(cnvCor) <- colnames(object)
    cnvCor[is.na( cnvCor )] <- 1
    getmode <- function(v) {
        uniqv <- unique(v)
        uniqv[which.max(tabulate(match(v, uniqv)))]
    }
    cnvCor <- t(t(cnvCor)/apply( cnvCor, 2, getmode ))
    cnvCor <- cnvCor/exp(rowMeans(log(cnvCor)))
    if( returnNormMat ){
        return(cnvCor)
    }
    object <- estimateSizeFactors( object, normMatrix=cnvCor )
    object
}

normalizeTrend <- function( object, prevNorm=TRUE ){
    logCounts <- log( counts( object, normalized=TRUE ) + 0.5 )
    ab <- log( rowMeans( counts( object, normalized=TRUE) + 0.5 ) )
    offs <- matrix(0, nrow(logCounts), ncol(logCounts), byrow = TRUE)
    for (x in seq_len(ncol(logCounts))) {
        fit <- loessFit(logCounts[, x], ab )
        offs[, x] <- fit$fitted
    }
    offs <- offs - rowMeans(offs)
    if( prevNorm ){
        normalizationFactors(object) <- normalizationFactors(object)*exp(offs)
    }else{
        normMatrix <- exp( offs )
        sf <- sizeFactors( object )
        nf <- t(t( normMatrix ) * sf)
        nf <- nf/exp( rowMeans( log( nf ) ) )
        normalizationFactors( object ) <- nf
    }
    object
}
