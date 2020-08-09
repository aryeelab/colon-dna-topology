
#library(HDF5Array)
#bs <- loadHDF5SummarizedExperiment("/data/aryee/bernstein/seqcapepi/SummarizedData/Robjects/se_hdf5")
#loops <- readRDS("/data/aryee/areyes/GB/glioma_topology/output/rds/loopsAll.rds")
#regions <- loops@anchors
#names( regions ) <- sprintf("anchors%0.6d", seq_len( length( regions ) ) )

getMeanMeth <- function( bs, regions, numCpg=0, subRegions=NULL, minCov=3 ){
    if( is.null( names(regions) ) ){ stop( "Regions need to be named") }
    ##    cov <- as.matrix( getCoverage( bs, type="Cov" ) )
    ##    M <- as.matrix( getCoverage( bs, type="M" ) )
    cat( "reading cov data\n" )
    cov <- as.matrix(getCoverage( bs, type="Cov" ))
    cat( "reading M data\n" )
    M <- as.matrix(getCoverage( bs, type="M" ))
    cat( "diviging\n" )
    cov[cov < minCov] <- NA
    meth <- M/cov
    rm(cov, M)
    gc()
    ##    meth <- M/cov
    ##    meth[cov < 3] <- NA
    cat( "granges data\n" )
    bsRanges <- rowRanges( bs )
    if( !is.null( subRegions ) ){
        keep <- findOverlaps( bsRanges, subRegions )
        keep <- sort( unique( queryHits( keep ) ))
        bsRanges <- bsRanges[keep]
        meth <- meth[keep,]
    }
    gc()
    ovl <- findOverlaps( bsRanges, regions )
    ovlSp <- split( queryHits(ovl), names(regions)[subjectHits(ovl)] )
    ovlSp <- ovlSp[lengths( ovlSp ) > numCpg]
    cat( "tmp matrix\n" )
    methMat <- matrix( NA, ncol=ncol(bs), nrow=length(regions) )
    colnames(methMat) <- colnames(bs)
    rownames(methMat) <- names(regions)
    methMatTmp <- t(vapply( ovlSp, function(x){
        colMeans(meth[x,,drop=FALSE], na.rm=TRUE)
    }, numeric(ncol(meth)) ))
    methMat[rownames(methMatTmp),colnames(methMatTmp)] <- methMatTmp
    methMat
}

getSumCov <- function( bs, regions ){
    if( is.null( names(regions) ) ){ stop( "Regions need to be named") }
    ##    cov <- as.matrix( getCoverage( bs, type="Cov" ) )
    ##    M <- as.matrix( getCoverage( bs, type="M" ) )
    cat( "reading cov data\n" )
    cov <- as.matrix(getCoverage( bs, type="Cov" ))
    bsRanges <- rowRanges( bs )
    ovl <- findOverlaps( bsRanges, regions )
    ovlSp <- split( queryHits(ovl), names(regions)[subjectHits(ovl)] )
    cat( "tmp matrix\n" )
    methMat <- matrix( NA, ncol=ncol(bs), nrow=length(regions) )
    colnames(methMat) <- colnames(bs)
    rownames(methMat) <- names(regions)
    methMatTmp <- t(vapply( ovlSp, function(x){
        colSums(cov[x,,drop=FALSE], na.rm=TRUE)
    }, numeric(ncol(cov)) ))
    methMat[rownames(methMatTmp),colnames(methMatTmp)] <- methMatTmp
    methMat
}

