countHicInRegions <- function( interactionFile, regions, chunkSize=10000000, keepChr,
                              #tmpFile="tmp_2343allValidPairs",
                              summarizedFile=TRUE,
                              minDistance=NULL, keepRegions=NULL, intra=FALSE ){
    i=0
    container=list()
    tmpFile <- tempfile()
    if( !summarizedFile ){
        system( sprintf( "gunzip -c %s | awk '$2 == $5 {print $2,$3,$6}' > %s",
                       interactionFile, tmpFile ) )
    }else{
        tmpFile <- interactionFile
    }
    eof=FALSE
    while(!eof){
        i <- i + 1
        cat( sprintf("Processing chunk %s\n", i) )
        allValidPairsFile <- fread( tmpFile,
                                   nrows=chunkSize,
                                   colClasses=c("character", "integer", "integer"),
                                   skip=chunkSize*(i-1) )
        rowNumb <- nrow( allValidPairsFile )
        allValidPairsFile <- allValidPairsFile[allValidPairsFile$V1 %in% keepChr,]
        if( !is.null( minDistance ) ){
            keep <- abs( allValidPairsFile$V3 - allValidPairsFile$V2 ) > minDistance
            allValidPairsFile <- allValidPairsFile[keep,]
        }
        r1 <- GRanges( allValidPairsFile$V1,
                      IRanges(start=allValidPairsFile$V2,
                              end=allValidPairsFile$V2) )
        r2 <- GRanges( allValidPairsFile$V1,
                      IRanges(start=allValidPairsFile$V3,
                              end=allValidPairsFile$V3) )
        if( !is.null( keepRegions ) ){
            over1 <- findOverlaps( r1, keepRegions )
            over2 <- findOverlaps( r2, keepRegions )
            over1 <- seq_len( length( r1 ) ) %in% queryHits( over1 )
            over2 <- seq_len( length( r2 ) ) %in% queryHits( over2 )
            over1 <- over1 & over2
            r1 <- r1[which(over1)]
            r2 <- r2[which(over1)]
        }
#        keep <- seqnames(r1) %in% keepChr #&
#            as.character( seqnames(r2) ) %in% keepChr &
#            as.character( seqnames(r1) ) == as.character( seqnames(r2) )
#        r1 <- r1[keep,]
#        r2 <- r2[keep,]
        over1 <- findOverlaps( r1, regions )
        over2 <- findOverlaps( r2, regions )
        inBoth <- intersect( queryHits( over1 ), queryHits( over2 ) )
        over1 <- over1[queryHits( over1 ) %in% inBoth,]
        over2 <- over2[queryHits( over2 ) %in% inBoth,]        
        stopifnot( all( queryHits( over2 ) == queryHits( over1 ) ) )
        mai <- subjectHits( over2 ) >= subjectHits( over1 )
        first <- ifelse( mai, subjectHits( over1 ), subjectHits( over2 ) )
        second <- ifelse( !mai, subjectHits( over1 ), subjectHits( over2 ) )
        container[[i]] <- data.frame( first=first, second=second ) %>%
            dplyr:::group_by(first, second) %>%
            dplyr:::summarize( n=dplyr::n() )
        eof <- ifelse( rowNumb==chunkSize, FALSE, TRUE )
    }
    if( !summarizedFile ){
        file.remove( tmpFile )
    }
    container <- do.call( rbind, container )
    container <- container %>%
        dplyr:::group_by( first, second ) %>%
        dplyr:::summarize( n=sum(n) )
    container
}

countHicIntraRegions <- function( interactionFile, regions, chunkSize=10000000, keepChr, verbose=TRUE ){
    i=0
    container=list()
    tmpFile <- interactionFile
    eof=FALSE
    while(!eof){
        i <- i + 1
        if( verbose ){
            cat( sprintf("Processing chunk %s\n", i) )
        }
        allValidPairsFile <- fread( tmpFile,
                                   nrows=chunkSize,
                                   drop=c(1, 4, 7, 8, 9, 10, 11, 12, 13),
                                   skip=chunkSize*(i-1), fill=TRUE)
        stopifnot( ncol(allValidPairsFile) == 4 )
        lastLine <- allValidPairsFile[nrow(allValidPairsFile),]
##        print(lastLine)
##                                   colClasses=c("character", "integer",
##                                                "character", "integer") )
        rowNumb <- nrow( allValidPairsFile )
        allValidPairsFile <- allValidPairsFile[allValidPairsFile$V2 %in% keepChr,]
        allValidPairsFile <- allValidPairsFile[allValidPairsFile$V5 %in% keepChr,]
        r1 <- GRanges( allValidPairsFile$V2,
                      IRanges( start=allValidPairsFile$V3,
                              end=allValidPairsFile$V3 ) )
        r2 <- GRanges( allValidPairsFile$V5,
                      IRanges( start=allValidPairsFile$V6,
                              end=allValidPairsFile$V6 ) )
        over1 <- findOverlaps( r1, regions )
        over2 <- findOverlaps( r2, regions )
        inBoth <- intersect( queryHits( over1 ), queryHits( over2 ) )
        over1 <- over1[queryHits( over1 ) %in% inBoth,]
        over2 <- over2[queryHits( over2 ) %in% inBoth,]        
        stopifnot( all( queryHits( over2 ) == queryHits( over1 ) ) )
        mai <- subjectHits( over2 ) >= subjectHits( over1 )
        first <- ifelse( mai, subjectHits( over1 ), subjectHits( over2 ) )
        second <- ifelse( !mai, subjectHits( over1 ), subjectHits( over2 ) )
        container[[i]] <- data.frame( first=first, second=second ) %>%
            dplyr:::group_by(first, second) %>%
            dplyr:::summarize(n=dplyr::n())
        if( verbose ){
            cat(sprintf("Read %s lines\n", rowNumb) )
        }
        if( !rowNumb == chunkSize ){
            xx <- system(sprintf("tail -1 %s", tmpFile), intern=TRUE)
            xx <- as.data.frame(matrix(strsplit(xx, "\t")[[1]][c(2, 3, 5, 6)],
                                       nrow=1))
            cat( sprintf("Processed %s lines for sample %s\n",
                         chunkSize*i + rowNumb,
                         basename(tmpFile) ))
            if( !all(as.data.frame(lastLine) == xx) ){
                stop( sprintf("Last line does not match for file %s", tmpFile ) )
            }
        }
        eof <- ifelse( rowNumb==chunkSize, FALSE, TRUE )
    }
    container <- do.call( rbind, container )
    container <- container %>%
        dplyr:::group_by( first, second ) %>%
        dplyr:::summarize( n=sum(n) )
    container
}

#' 
#' @import GenomicRanges
getMatrixFromCont <- function( chr, merged, container, keepChr, symmetric=TRUE ){
    mergedSub <- merged[as.character(seqnames(merged)) %in% keepChr,]
    mergedSub <- mergedSub[as.character(seqnames(mergedSub)) == chr]
    matSize <- sum(as.character(seqnames( merged )) == chr)
    matBlock <- matrix( 0, nrow=matSize, ncol=matSize )
    chrIndexes <- which( as.vector( seqnames( merged ) == chr ) )
    containerSub <- container[container$first %in% chrIndexes,]
    containerSub <- as.data.frame(containerSub)
    minIdx <- min(chrIndexes)
    containerSub$first <- containerSub$first - minIdx + 1
    containerSub$second <- containerSub$second - minIdx + 1
    allFirsts <- sort( unique( containerSub$first ) )
    for( i in allFirsts ){
        if( i %% 20 == 0 ){
            cat( sprintf("processed %s rows\n", i ) )
        }
        toFill <- containerSub$first == i
        matBlock[i,containerSub$second[toFill]] <- containerSub$n[toFill]
    }
    if( symmetric ){
        matBlock[lower.tri( matBlock )] <- t( matBlock )[lower.tri( matBlock )]
    }
    matBlock
}

getAggregateMat <- function( interactFile, query, nBins, offset, aggregateFun=mean, asMatrix=TRUE, chunkSize=10000000, widthScaling=FALSE, distNormFun=mean ){
    if( is( query, "list" ) ){
        nBins2 <- (nBins/2) + offset
        queryLeft <- query[[1]]
        queryRight <- query[[2]]
        queryO <- queryLeft
        off <- round( ( width( queryLeft ) / (nBins/2) ) * offset )
        start( queryLeft ) <- start( queryLeft ) - off
        off <- round( ( width( queryRight ) / (nBins/2) ) * offset )
        end( queryRight ) <- end( queryRight ) + off
        queryLeft <- trim( queryLeft )
        queryRight <- trim( queryRight )
        queryLeft <- tile( queryLeft, n=nBins2 )
        queryRight <- tile( queryRight, n=nBins2 )
        queryLeft <- unlist( queryLeft )
        mcols( queryLeft )$blockIndex <- rep( seq_len( length( queryO ) ), each=nBins2 )
        mcols( queryLeft )$binNum <- rep( seq_len( nBins2 ), length( queryO ) )
        queryRight <- unlist( queryRight )
        mcols( queryRight )$blockIndex <- rep( seq_len( length( queryO ) ), each=nBins2 )
        mcols( queryRight )$binNum <- rep( nBins2+seq_len( nBins2 ), length( queryO ) )
        query <- c( queryLeft, queryRight )
        nBins2 <- nBins2*2
#        return(query)
    }else{
        off <- round( ( width( query ) / nBins ) * offset )
        nBins2 <- nBins + (offset*2)
        queryO <- query
        start( query ) <- start( query ) - off
        end( query ) <- end( query ) + off
        query <- trim( query )
        query <- tile( query, n=nBins2 )
        query <- unlist( query )
        mcols( query )$blockIndex <- rep( seq_len( length( queryO ) ), each=nBins2 )
        mcols( query )$binNum <- rep( seq_len( nBins2 ), length( queryO ) )
    }
    selfOver <- findOverlaps( query )
    selfOver <- selfOver[queryHits( selfOver ) != subjectHits( selfOver ),]
    query <- query[!seq_along(query) %in% unique(queryHits( selfOver ))]
    pr <- findOverlaps( query )
    stopifnot(all( queryHits(pr) == subjectHits(pr) ))
    keepChr <- unique( as.character( seqnames( query ) ) )
    counts <- countHicInRegions( interactFile, query, keepChr=keepChr, chunkSize=chunkSize )
    counts$binLeft <- query$binNum[counts$first]
    counts$binRight <- query$binNum[counts$second]
    counts$indexLeft <- query$blockIndex[counts$first]
    counts$indexRight <- query$blockIndex[counts$second]
    counts <- counts[counts$indexLeft == counts$indexRight,]
    counts$binWidth <- width(query)[counts$first]
    counts$group <- queryO$group[counts$indexLeft]
    if( widthScaling ){
        counts$interactScore <- counts$n / counts$binWidth
        counts$interactScore <- ( counts$interactScore / sum( counts$interactScore ) ) * 10^6
    }else{
        counts$interactScore <- counts$n
    }
    score <- counts %>%
        dplyr:::group_by( binLeft, binRight, group ) %>%
        dplyr:::summarize( intScore=aggregateFun(interactScore) ) %>%
        as.data.frame
#    score$intScore <- log10( score$intScore )
    score$distance <- score$binRight - score$binLeft
    if( !is.null(distNormFun) ){
        normFact <- score %>%
            dplyr:::group_by( distance, group ) %>%
            dplyr:::summarize( normFactor=distNormFun(intScore) )
        score <- dplyr:::left_join( score, normFact, by = c( "distance", "group" ) )
        score <- dplyr::mutate( score, normIntScore = intScore / normFactor )
    }else{
        score <- dplyr::mutate( score, normIntScore = intScore )
    }
    if( asMatrix ){
        allMats <- lapply( unique(score$group), function(xx){
            scoreSub <- dplyr::filter( score, group == xx )
            mat <- matrix( NA, ncol=nBins2, nrow=nBins2 )
            for( i in seq_len( nrow( scoreSub ) ) ){
                x1 <- scoreSub$binLeft[i]
                x2 <- scoreSub$binRight[i]
                mat[x1,x2] <- scoreSub$normIntScore[i]
                mat[x2,x1] <- mat[x1,x2]
            }
            mat
        } )
        names( allMats ) <- unique( score$group )
        return( allMats )
    }else{
        return( dplyr::select( score, binLeft, binRight, group, intScore, normIntScore ) )
    }
}

countHiCRectangles <- function( interactionFile, regionRow, regionCol,
                               keepChr=NULL, chunkSize=50000000, #tmpFile="tmpFile20129.txt",
                               summarizedFile=TRUE ){
    stopifnot( all( as.character(seqnames( regionRow )) %in% keepChr ) )
    stopifnot( all( as.character(seqnames( regionCol )) %in% keepChr ) )
    i=0
    container=list()
    if( !summarizedFile ){
        cat("summarizing\n")
        tmpFile <- tempfile()
        system(sprintf("gunzip -c %s | awk '$2 == $5 {print $2,$3,$6}' > %s",
                       interactionFile, tmpFile) )
    }else{
        tmpFile <- interactionFile
    }
    eof=FALSE
    while(!eof){
        i <- i + 1
        cat( sprintf("Processing chunk %s\n", i) )
        allValidPairsFile <- fread( tmpFile,
                                   nrows=chunkSize,
                                   colClasses=c("character", "integer", "integer"),
                                   skip=chunkSize*(i-1) )
        rowNumb <- nrow( allValidPairsFile )
        allValidPairsFile <- allValidPairsFile[allValidPairsFile$V1 %in% keepChr,]
        r1 <- GRanges( allValidPairsFile$V1,
                      IRanges(start=allValidPairsFile$V2,
                              end=allValidPairsFile$V2) )
        r2 <- GRanges( allValidPairsFile$V1,
                      IRanges(start=allValidPairsFile$V3,
                              end=allValidPairsFile$V3) )
        over1 <- findOverlaps( r1, regionRow )
        over1 <- over1[!duplicated( queryHits( over1 ) )]
        over2 <- findOverlaps( r2, regionCol )
        over2 <- over2[!duplicated( queryHits( over2 ) )]
        overP1 <- findOverlaps( r2, regionRow )
        overP1 <- overP1[!duplicated( queryHits( overP1 ) )]
        overP2 <- findOverlaps( r1, regionCol )
        overP2 <- overP2[!duplicated( queryHits( overP2 ) )]
        first <- c( subjectHits( over1 ), subjectHits( overP1 ) )
        second <- c( subjectHits( over2 ), subjectHits( overP2 ) )
        container[[i]] <- data.frame( first=first, second=second ) %>%
            dplyr:::group_by( first, second ) %>%
            dplyr:::summarize( n=dplyr::n() )
        eof <- ifelse( rowNumb==chunkSize, FALSE, TRUE )
    }
    if( !summarizedFile ){
        file.remove( tmpFile )
    }
    container <- do.call( rbind, container )
    container <- container %>%
        dplyr:::group_by( first, second ) %>%
        dplyr:::summarize( n=sum(n) )
    container
}

abRatio <- function( interactionsFile, regionRow, regionCol, compVector, keepChr ){
    rects <- countHiCRectangles( interactionsFile, regionRow=regionRow,
                                regionCol=regionCol, keepChr=keepChr )
    rects$chr <- as.character( seqnames(regionCol) )[rects$second]
    rectsSp <- split( rects, rects$chr )
    thresholds <- tapply( compVector, as.character( seqnames( regionCol ) ),
                         function(x){
                             quantile( x, c(0.2, 0.8), na.rm=TRUE )[c("20%", "80%")]
                         } )
    overMask <- findOverlaps( regionRow, regionCol )
    for( i in keepChr ){
        veryB <- thresholds[[i]][["20%"]]
        veryA <- thresholds[[i]][["80%"]]
        vecVal <- compVector[rectsSp[[i]]$second]
        rectsSp[[i]] <- rectsSp[[i]][which(vecVal < veryB | vecVal > veryA),]
        rectsSp[[i]]$compartment <- factor( ifelse( compVector[rectsSp[[i]]$second] > 0, "A", "B" ) )
    }
    rectsSp <- do.call(rbind, rectsSp)
    abCounts <- t( sapply( unique( rectsSp$first ), function( indChr ){
        toMask <- subjectHits( overMask )[queryHits( overMask ) == indChr]
        contRSub <- rectsSp[rectsSp$first == indChr,]
        contRSub <- contRSub[!contRSub$second %in% seq( toMask[1] - 20, toMask[1] + 20 ),]
        if( nrow( contRSub ) > 0 ){
            return( tapply( contRSub$n, contRSub$compartment, sum, default=0 ) )
        }else{
            return( as.array(c("A"=0, "B"=0 ) ))
        }
    } ) )
    rownames(abCounts) <- names(regionRow)[unique( rectsSp$first )]
    abCountsEmpty <- matrix( 0, ncol=2, nrow=length( regionRow ) )
    rownames( abCountsEmpty ) <- names( regionRow )
    colnames( abCountsEmpty ) <- colnames( abCounts )
    abCountsEmpty[rownames(abCounts),] <- abCounts
    abCountsEmpty <- as.data.frame( log2( abCountsEmpty + 1 ) )
    abCountsEmpty$abRatio <- abCountsEmpty$A - abCountsEmpty$B
    abCountsEmpty[names(regionRow),]
}

abRatio2 <- function( interactionsFile, regionRow, regionCol, compVector, keepChr ){
    rects <- countHiCRectangles( interactionsFile, regionRow=regionRow,
                                regionCol=regionCol, keepChr=keepChr, summarizedFile=FALSE )
    rects$chr <- as.character( seqnames(regionCol) )[rects$second]
    allChrData <- lapply( keepChr, function( chr ){
        cat( sprintf( "Processing %s\n", chr ) )
        colIndexes <- as.character(seqnames( regionCol )) %in% chr
        rowIndexes <- as.character(seqnames( regionRow )) %in% chr
        rectsSub <- rects[rects$chr %in% chr,]
        minColIndex <- min( which(colIndexes) )
        minRowIndex <- min( which(rowIndexes) )
        mat <- matrix( 0, ncol=sum( colIndexes ), nrow=sum( rowIndexes ) )
        rectsSub$first <- rectsSub$first - minRowIndex + 1
        rectsSub$second <- rectsSub$second - minColIndex + 1
        mat <- matrix( 0, ncol=sum( colIndexes ), nrow=sum( rowIndexes ) )
        firstVec <- rectsSub$first
        secondVec <- rectsSub$second
        countsVec <- rectsSub$n
        cat( sprintf("Building matrix for chromosome %s\n", chr ) )
        for( i in sort( unique( secondVec ) ) ){
            auxInd <- secondVec == i
            mat[firstVec[auxInd],i] <- countsVec[auxInd]
        }
        compVectorSub <- compVector[colIndexes]
        thresholds <- quantile( compVectorSub, c(0.2, 0.8), na.rm=TRUE )[c("20%", "80%")]
        veryB <- compVectorSub < thresholds[["20%"]]
        veryA <- compVectorSub > thresholds[["80%"]]
        ovl <- findOverlaps( regionRow[rowIndexes,], regionCol[colIndexes,] )
        ovl <- as.data.frame( ovl )
        cat( sprintf( "Calculating ratio for %s\n", chr ) )
        abRatio <- t( sapply( seq_len( nrow( mat ) ), function(i){
            mask <- ovl$subjectHits[ovl$queryHits == i][1]
            mask <- !seq_len( ncol( mat ) ) %in% seq( mask - 25, mask + 25 )
            c( A=sum( mat[i,which(veryA & mask)] ), B=sum( mat[i,which(veryB & mask)] ) )
        } ) )
        abRatio <- as.data.frame( log2( abRatio + 1 ) )
        abRatio$abRatio <- abRatio$A - abRatio$B
        abRatio
    } )
    allChrData <- do.call( rbind, allChrData )
    rownames(allChrData) <- names(regionRow)
    allChrData
}

## require(data.table)
## require(GenomicRanges)
## require(magrittr)
## require(fastmatch)

## @importFrom fastmatch fmatch
## @importFrom data.table fread
##`%fin%` <- function(x, table) {
##   fmatch(x, table, nomatch = 0L) > 0L
##}

getBoundaryScores <- function( countsLongSub, coords, howmuch=200000, resolution=1000, chromosome="chr1" ){
    chrIndexes <- which(as.vector( seqnames( coords ) == chromosome ))
    nBins <- howmuch/resolution
    countsLongSub <- countsLongSub[countsLongSub$first %in% chrIndexes,]
    countsLongSub <- as.data.frame(countsLongSub)
    countsLongSub <- countsLongSub[(countsLongSub$second - countsLongSub$first) < (nBins*2 + 3),]
    maxi <- length(chrIndexes)
    iter <- seq_len( maxi- nBins*2 + 1 )
    prueba2 <- lapply( iter, function(i){
        if( i %% 200 == 0 ){
            cat(sprintf("processed %s bins\n", i))
        }
        left <- chrIndexes[i:(i+nBins-1)]
        mm <- i+nBins*2-1
        right <- chrIndexes[(i+nBins):mm]
        allInd <- c( left, right )
        rest <- min( allInd )
        countsLongSub <- countsLongSub[countsLongSub$first %in% allInd & countsLongSub$second %in% allInd,]
        lvls <- as.character(seq_len( nBins*2 ))
        countsLongSub$first <- factor(countsLongSub$first - rest +1, levels=lvls )
        countsLongSub$second <- factor(countsLongSub$second - rest + 1, levels=lvls )
        mt <- with( countsLongSub, {
            mt <- matrix( 0, nrow=nlevels(first), ncol=nlevels(second),
                         dimnames=list(levels(first), levels(second) ) )
            mt[cbind(first, second)] <- n
            mt
        } )
        right <- as.character(right - rest + 1)
        left <- as.character(left - rest + 1)
        c(a=sum(mt[left,left]), b=sum(mt[right,right]), c=sum(mt[left,right]))
    } )
    compBins <- getCompBinning( resolution, chromosome )
    compBins <- compBins[start(compBins) > howmuch]
    compBins <- compBins[iter]
    compBins <- shift( compBins, -1*resolution/2 )
    mcols(compBins) <- do.call(rbind, prueba2)
    compBins
}

getLoopData <- function( samp, bins, keepChr, allInteractFiles, ctcfLoopsLeft, ctcfLoopsRight, offset=5 ){
    bins <- bins[as.character(seqnames(bins)) %in% keepChr]
    countsLong <- countHicInRegions( allInteractFiles[samp], bins, chunkSize=50000000, summarizedFile=FALSE, keepChr=keepChr )
    keepRun <- keepChr
#    nJobs <- length(keepRun)
#    resourcesList <- list("queue" = "big", "memory"="25000", "ncpus"="1", "memorylimit"="25000" )
#    source("~/cluster/useCluster.R")
    allLoops <-
        lapply(
            keepRun,
            function( chr, bins, countsLong, keepChr, offset, ctcfLoopsLeft, ctcfLoopsRight ){
#                source("/data/aryee/areyes/GB/glioma_topology/code/support/countHiCReads.R")
                cat(sprintf("processing %s...\n", chr) )
                mat1 <- getMatrixFromCont( chr, bins, countsLong, keepChr )
                cat(sprintf("finshed building matrix %s...\n", chr) )
                binsSub <- bins[as.character( seqnames( bins )) %in% chr]
                maxIdx <- length( binsSub )
                ctcfLoopsLeftSub <- ctcfLoopsLeft[as.character(seqnames( ctcfLoopsLeft )) %in% chr]
                ctcfLoopsRightSub <- ctcfLoopsRight[as.character(seqnames( ctcfLoopsRight)) %in% chr]
                leftBin <- subjectHits( findOverlaps( ctcfLoopsLeftSub, binsSub ) )
                rightBin <- subjectHits( findOverlaps( ctcfLoopsRightSub, binsSub ) )
                names( leftBin ) <- names(ctcfLoopsLeftSub)
                names( rightBin ) <- names(ctcfLoopsRightSub)
                stopifnot( all( names(leftBin) == names(rightBin) ) )
                countsChr <- lapply( names( leftBin ), function(x){
                    leftBinAll <- seq( leftBin[x] - offset, leftBin[x] + offset )
                    rightBinAll <- seq( rightBin[x] - offset, rightBin[x] + offset )
                    if( !all(leftBinAll > 0 & leftBinAll <= maxIdx) ){
                        return(NULL)
                    }
                    if( !all(rightBinAll > 0 & rightBinAll <= maxIdx) ){
                        return(NULL)
                    }
                    countsSub <-
                        reshape2::melt( mat1[leftBinAll,rightBinAll], value.name="n",
                                       varnames=c("first", "second") )
                    countsSub$first <- countsSub$first - offset - 1
                    countsSub$second <- countsSub$second - offset - 1
                    countsSub$loopID <- x
                    countsSub
                })
                countsChr <- do.call(rbind, countsChr)
                countsChr
            }, bins=bins, countsLong=countsLong,
            keepChr=keepChr, offset=offset, ctcfLoopsLeft, ctcfLoopsRight )
    allLoops <- do.call( rbind, allLoops )
    allLoops$sample <- samp
    allLoops
}
