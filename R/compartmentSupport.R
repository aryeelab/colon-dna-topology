
getObservedVsExpected <- function( prueba ){
    numCols <- ncol(prueba)
    prueba2 <- as( prueba, "sparseMatrix" )
    indx1 <- 1:ncol(prueba)
    for( i in 1:( numCols - 1 )){
        bnd <- band( prueba2, i, i )
        sumBand <- sum( bnd, na.rm=TRUE )
        indx2 <- indx1 + i
        inRange <- indx2 <= numCols
        indx1 <- indx1[inRange]
        indx2 <- indx2[inRange]
        meanBand <- sumBand / sum( inRange )
        if( sumBand > 0 ){
            for( k in seq_len( length(indx1) ) ){
                mr <- prueba[indx1[k],indx2[k]]
                if( is.na( mr ) ){
                    next()
                }else if( mr == 0 ){
                    next()
                }else{
                    prueba[indx1[k],indx2[k]] <- mr / meanBand
                    prueba[indx2[k],indx1[k]] <- mr / meanBand
                }
            }
        }
    }
    diag(prueba) <- diag(prueba) / mean( diag(prueba) )
    as.matrix(prueba)
}

getObservedVsExpected2 <- function( prueba ){
    numCols <- ncol(prueba)
    prueba2 <- as( prueba, "sparseMatrix" )
    indx1 <- 1:ncol(prueba)
    for( i in 0:( numCols - 1 )){
        bnd <- as( band( prueba2, i, i ), "sparseMatrix" )
        sumBand <- sum( bnd )
        hasValues <- as.matrix(bnd > 0)
        num <- numCols - i
        meanBand <- sumBand / num
        prueba[hasValues] <- prueba[hasValues] / meanBand
    }
    as.matrix(prueba)
}

#' @export minfi
estimateCompartmentForCorMat <- function( prueba, what="eigen" ){
    keep <- colSums( prueba ) > 0
    corMat <- cor( prueba[keep,keep] )
    a <- rep(NA, nrow(prueba))
    b <- minfi:::.getFirstPC(corMat, method="exact")
    a[keep] <- b
    a <- a * sqrt(length(a))
    if( what == "eigen" ){
        a
    }else{
        ifelse( a > 0, "A", "B" )
    }
}

estimateCompartmentsForSample <- function( hic, res="150000", chr="chr14"){
    Oij <- hic@resolutionNamedList[[res]][[chr]]
    stopifnot( nrow(Oij) == ncol(Oij) )
    Oij <- as.matrix(Oij)
    prueba <-  getObservedVsExpected( Oij )
    rt <- as.vector( estimateCompartmentForCorMat( prueba ) )
    rt
}

#' @importFrom Matrix band
#' @export
estimateCompartmentsForSamples <- function( sampleNames, chrs, res="100000" ){
    allEigens <- bplapply( sampleNames, function(y){
        cat( sprintf( "Reading sample %s...\n", y ) )
        data(sprintf("hic_%s", y))
        hic <- get(sprintf("hic_%s", y))
        allComps <- lapply( chrs, function(x){
            cat( sprintf( "processing in res %s chromosome %s in sample %s\n", res, x, y) )
            estimateCompartmentsForSample( hic, res=res, chr=x )
        } )
        names( allComps ) <- chrs
        allComps
    } )
    names( allEigens ) <- sampleNames
    allEigens <- lapply( chrs, function(x){
        sapply( allEigens, "[[", x )
    })
    names( allEigens ) <- chrs
    allEigens
}
