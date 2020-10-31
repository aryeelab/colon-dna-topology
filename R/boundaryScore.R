
calculateBoundaryScores <- function( mats, res, steps ){
    boundaryScores <- dplyr::bind_rows( bplapply( mats, function(smp){
##        mt <- readRDS( mats[smp] )
        mt <- get(sprintf("hic_%s", smp))
        allChrs <- names( mt@resolutionNamedList[[res]] )
        res <- dplyr::bind_rows( lapply( allChrs, function(chr){
            mat <- as.matrix( mt@resolutionNamedList[[res]][[chr]] )
            mat <- getObservedVsExpected( mat )
            mat[lower.tri(mat)] <- 0
            res <- dplyr::bind_rows(
                              lapply( steps,
                                     function(binNum){
                                         cat(sprintf("Sample %s, chromosome %s, bin %s\n", smp, chr, binNum))
                                         b <- sparseHiC:::.boundary_score_matrix( mat, num_bins=binNum )
                                         data.frame( sampleID=smp, chr=chr, pos=b$pos,
                                                    within_vs_cross=b$within_vs_cross, bin_num=binNum )
                                     } ) )
            res
        } ) )
    } ))
    boundaryScores
}

getTiledScores <- function( boundary, any=TRUE, refSamp=NULL, distanceThr=300000, numBins=11,
                           na.rm=TRUE, meanCenter=TRUE ){
    boundMat <- as.matrix(boundary[,!colnames( boundary ) %in% c("chr", "pos"),])
    boundMat <- (t((t(boundMat) - colMedians( boundMat, na.rm=TRUE ))/ colMads( boundMat, na.rm=TRUE) ) )
    boundary <- cbind( chr=boundary$chr, pos=boundary$pos, as.data.frame( boundMat ) )
    grs <- GRanges( boundary$chr, IRanges( boundary$pos-20000 + 1, boundary$pos+20000 ) )
    seqlengths(grs) <- seqlengths(blocks)[names(seqlengths(grs))]
    grs <- trim(grs)
    boundary_cutoff <- qnorm(0.9)
    if( any ){
        boundPass <- which(rowSums( boundMat > boundary_cutoff ) > 0)
        avgScore <- rowMeans( boundMat, na.rm=TRUE )
    }else{
        boundPass <- which( boundMat[,refSamp] > boundary_cutoff )
        avgScore <- boundMat[,refSamp]
    }
    boundaryGrp <- reduce(grs[boundPass,])
    ovl <- findOverlaps(grs[boundPass,], boundaryGrp)
    ovl <- split( boundPass[queryHits( ovl )], subjectHits( ovl ) )
    boundPass2 <- vapply( ovl, function(x){ x[which.max(avgScore[x])] }, numeric(1))
    names(boundPass2) <- NULL
    boundaryUnion <- grs[boundPass2,]
    nrs <- distanceToNearest(boundaryUnion)
    nrs <- nrs[mcols(nrs)$distance > distanceThr,]
    boundaryUnion <- boundaryUnion[unique(c(queryHits(nrs), subjectHits(nrs))),]
    boundaryUnion <- shift(resize( resize( boundaryUnion, 1, fix="center" ), 40000*numBins, fix="center"), 1)
    names( boundaryUnion ) <- sprintf( "boundary%0.9d", seq_len(length(boundaryUnion)) )
    boundaryUnion <- unlist(tile( boundaryUnion, numBins))
    mcols(boundaryUnion)$id <- sapply(strsplit(names(boundaryUnion), "\\."), "[[", 1)
    mcols(boundaryUnion)$binNum <- seq_len(numBins)
    names(boundaryUnion) <- NULL
    ovl <- findOverlaps( boundaryUnion, grs )
    ovl <- split( subjectHits(ovl), queryHits( ovl ) )
    stopifnot(all(lengths( ovl ) == 1))
    tiledBoundaryScores <- lapply( ovl, function(x){
        colMeans(boundMat[x,,drop=FALSE], na.rm=TRUE)
    } )
    tiledBoundaryScores <- do.call(rbind, tiledBoundaryScores)
    tiledBoundaryScores <- cbind(as.data.frame( mcols(boundaryUnion[as.numeric(names(ovl))])), tiledBoundaryScores)
#    tiledBoundaryScores <- tiledBoundaryScores[rowSums( is.na( tiledBoundaryScores) ) == 0,]
    tiledBoundaryScores <- reshape2::melt( tiledBoundaryScores,
                                          id.vars=c("id", "binNum"),
                                          variable.name="sample",
                                          value.name="score")
    if( meanCenter ){
        tiledBoundaryScores <- tiledBoundaryScores %>%
            dplyr::group_by( id, sample ) %>%
            dplyr::mutate( score=score - mean(score, na.rm=TRUE) ) %>%
            dplyr::ungroup( )
    }
    if( na.rm ){
        tiledBoundaryScores <- na.omit( tiledBoundaryScores )
        tiledBoundaryScores <- tiledBoundaryScores %>%
            dplyr::group_by( id, sample ) %>%
            dplyr::filter( dplyr::n() == numBins )
    }
##    sum(is.na( tiledBoundaryScores$score ))
    tiledBoundaryScores
}
