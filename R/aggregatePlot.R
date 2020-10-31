plotLoopAggregate <- function( loopData, keepLoops, disruptedLoops, bias=1, levs=NULL, furtherAgg=NULL, maxVal=NULL, furtherLevs=NULL, zlim ){
    if( !is.null( maxVal ) ){
        loopData <- loopData %>%
            dplyr::filter( first %in% seq(-maxVal, maxVal, 1),
                          second %in% seq(-maxVal, maxVal, 1) )
    }
    loopData <- loopData %>%
        dplyr::filter( loopID %in% keepLoops ) %>%
        dplyr::group_by( sample, loopID ) %>%
        dplyr::mutate( nMedian=mean(n), nNorm=(n / mean(n)) )
    loopData <- dplyr::ungroup( loopData )
    datSum <- loopData %>%
        dplyr::group_by( sample, loopID ) %>%
        dplyr::summarize( med=unique(nMedian) ) %>%
        dplyr::ungroup() %>%
        dplyr::group_by( loopID ) %>%
        dplyr::summarize( sm=sum( med > 0 ) )
    keepCov <- datSum$loopID[datSum$sm == length( unique(loopData$sample ))]
    distNorm <- loopData %>%
        dplyr::filter( !loopID %in% disruptedLoops ) %>%
        dplyr::mutate( dist=sqrt( first^2 + second^2 ) ) %>%
        dplyr::group_by( sample, dist ) %>%
        dplyr::summarize( normFac=median( nNorm, na.rm=TRUE ) ) %>%
        dplyr::ungroup()
    distNormNorm <- distNorm %>%
        dplyr::group_by( dist ) %>%
        dplyr::summarize( distNorm=median(normFac) )
    distNorm <- dplyr::full_join( distNorm, distNormNorm, by="dist" ) %>%
        dplyr::mutate( normFac=normFac-distNorm )
    loopData <- loopData %>%
        dplyr::mutate( dist=sqrt( first^2 + second^2 ) ) %>%
        dplyr::full_join( distNorm, by=c("sample", "dist") )
    loopData$nNorm <- loopData$nNorm - loopData$normFac
    loopData <- loopData %>%
        dplyr::filter( loopID %in% keepCov ) %>%
        dplyr::mutate(
                   group=factor( ifelse( loopID %in% disruptedLoops, "Lost", "Stable" ),
                                levels=c("Stable", "Lost") ) ) %>%
        dplyr::group_by( first, second, sample, group ) %>%
        dplyr::summarize( aggScore=mean(nNorm, na.rm=TRUE) ) %>%
        as.data.frame
    if( !is.null( levs ) ){
        loopData$sample <- factor( as.character( loopData$sample ), levels=levs )
    }
    if( !is.null( furtherAgg ) ){
        loopData$sample <- furtherAgg[loopData$sample]
        loopData <- loopData %>%
            dplyr::group_by( second, first, group, sample ) %>%
            dplyr::summarize( aggScore=mean( aggScore ) )
        if( !is.null( furtherLevs ) ){
            loopData$sample <- factor( loopData$sample, levels=furtherLevs )
        }
    }
    if( !is.null( zlim ) ){
        loopData$aggScore <- pmin( loopData$aggScore, zlim )
    }
    loopData %>%
        ggplot( aes( rev(second), rev(first), fill=aggScore ) ) +
        geom_tile() + facet_grid( group~sample ) +
        coord_fixed() +
        scale_fill_gradientn(
            colors=( colorRampPalette( brewer.pal(9, "Blues"), bias=bias )(100)),
            na.value="grey50" )
}


getDESeqRes <- function( pairCountList=NULL, anchorIndexes=NULL, hicDsd=NULL,
                        thr=5, thr2=2, normalizeTrend=FALSE ){
    if( is.null( hicDsd ) ){
        hicSub <- list()
        for( i in names(pairCountList)){
            kp1 <- pairCountList[[i]]$first %in% anchorIndexes
            kp2 <- pairCountList[[i]]$second %in% anchorIndexes
            hicSub[[i]] <- pairCountList[[i]][kp1&kp2,]
        }
        hicSub <- do.call(rbind, hicSub)
        hicSub <- reshape2::dcast( data=hicSub, first + second ~ sample, value.var="n", fill=0 )
        cols <- !grepl("first|second", colnames(hicSub))
        kp <- rowSums( hicSub[,cols] > thr2 ) > thr2 & hicSub$first != hicSub$second
        hicSub <- hicSub[kp,]
        rownames(hicSub) <- NULL
        hicDsd <- DESeqDataSetFromMatrix(
            countData = hicSub[,cols],
            colData = data.frame(
                sample = colnames(hicSub)[cols],
                groups = ifelse( colnames(hicSub)[cols] %in% c("BRD3162", "BRD3179"), "CIMP", "NonCIMP" ) ),
            design = ~groups )
        rowData(hicDsd)$first <- hicSub$first
        rowData(hicDsd)$second <- hicSub$second
    }
    hicDsd <- estimateSizeFactors( hicDsd )
    hicDsd <- hicDsd[rowMeans( counts( hicDsd, normalized=TRUE ) ) > thr,]
    hicDsd <- estimateSizeFactors( hicDsd )
    if( normalizeTrend ){
        hicDsd <- normalizeTrend( hicDsd, prevNorm=FALSE )
    }
    hicDsd <- DESeq( hicDsd )
    meanExpr <- vapply( split( seq_len( ncol( hicDsd ) ), colData( hicDsd )$groups ),
           function(x){
               rowMeans( counts( hicDsd, normalized=TRUE )[,x] )
           }, FUN.VALUE=numeric(nrow(hicDsd)) )
    sumExpr <- vapply( split( seq_len( ncol( hicDsd ) ), colData( hicDsd )$groups ),
                       function(x){
                           rowSums( counts( hicDsd, normalized=FALSE )[,x] )
                       }, FUN.VALUE=numeric(nrow(hicDsd) ) )
    sf2 <- estimateSizeFactorsForMatrix( sumExpr )
    sumExpr <- t(t( sumExpr ) / sf2)
    rawLfc <- log2( meanExpr[,"CIMP"] / meanExpr[,"NonCIMP"] )
    rawLfc2 <- log2( sumExpr[,"CIMP"] / sumExpr[,"NonCIMP"] )
    hicRes <- results( hicDsd, contrast=c("groups", "CIMP", "NonCIMP") )
    rs <- fdrtool(results(hicDsd, contrast=c("groups", "CIMP", "NonCIMP"))$stat,
                  statistic= "normal", plot = FALSE)
    hicRes$log2FoldChangeRaw <- rawLfc
    hicRes$log2FoldChangeSum <- rawLfc2
    hicRes$pvalue <- rs$pval
    hicRes$padj <- p.adjust( rs$pval, method="BH" )
    hicRes$first <- rowData(hicDsd)$first
    hicRes$second <- rowData(hicDsd)$second
    hicRes$sampNum <- rowSums( counts(hicDsd) > 0  )
    rowData(hicDsd) <- hicRes
    hicDsd
}
