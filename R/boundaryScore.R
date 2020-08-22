
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
