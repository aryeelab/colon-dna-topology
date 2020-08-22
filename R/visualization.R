plotCoverage <- function( bamDir, plot_gr, sizeFactors=NULL, fun=function(x){x},
                         nm="", nms="", smooth=0, norm=FALSE, facet=FALSE,
                         highlight= NULL, sampleGroups=NULL, meanSmooth=FALSE,
                         ksmooth=1, itersmooth=1, binSize=100, lwd=0.2){
    names( sampleGroups ) <- nms
    bfl <- BamFileList( bamDir )
    names( bfl ) <- nms
    if( !is.null( sizeFactors ) ){
        names( sizeFactors ) <- nms
    }
    seqlevelsStyle(plot_gr) <- seqlevelsStyle( bfl )
    start( plot_gr ) <- start( plot_gr ) - smooth - 1
    end( plot_gr ) <- end( plot_gr ) + smooth
    sbp <- ScanBamParam( which=plot_gr )
    chr <- as.character(seqnames(plot_gr))
    df <- do.call( rbind, lapply( names(bfl), function(w){
        cov <- coverage( readGAlignments( bfl[[w]], param=sbp ) )[chr]
        print(names(cov))
        bins <- unlist(tile(plot_gr, n=width(plot_gr)/binSize))
        bins <- keepSeqlevels(bins, chr)
        bins <- binnedSum( bins, cov, "cov" )
        cov <- bins$cov/200
#        cov <- as.numeric(cov)[start(plot_gr):end(plot_gr)]
        if( !is.null( sizeFactors ) ){
            cov <- cov/sizeFactors[w]
        }else if( norm ){
            cov <- cov / mean(cov, na.rm=TRUE)
        }
        cov <- pmax( fun(cov), 0 )
        if( smooth > 0){
            cov <- stats::filter(cov, rep(1/smooth, smooth))
        }
        if( meanSmooth ){
            cov <- minfi:::.meanSmoother( cov, ksmooth, itersmooth )
        }
        df <- data.frame(
            genomicPos=rowMeans( data.frame( start( bins ), end( bins ) ) ),
            cov=cov,
            sample=w )
        df <- na.omit(df)
        df
    } ) )
    if( is.null( sampleGroups ) ){
        df$sampleGroups <- df$sample
    }else{
        df$sampleGroups <- sampleGroups[df$sample]
    }
    p4 <- ggplot( df, aes( genomicPos, cov, col=sampleGroups, group=sample ) ) +
        ylab( nm ) +
        geom_line(lwd=lwd)
    if( !is.null( highlight ) ){
        highlight <- subsetByOverlaps( highlight, plot_gr )
        p4 <- p4 +
            geom_rect( data = as.data.frame( highlight ), inherit.aes=FALSE,
                      aes( xmin=start, xmax=end, ymin=-Inf, ymax=Inf ),
                      fill="gray", alpha=0.35 )
    }
    if( facet ){
        p4 <- p4 + facet_grid( sample ~ . )
    }
    p4 <- p4 + xlim( start( plot_gr ), end( plot_gr ) )
    p4
}


plotCoverage2 <- function( bamDir, plot_gr, sizeFactors=NULL, fun=function(x){x},
                          nm="", nms="", smooth=0, norm=FALSE, facet=FALSE,
                          highlight= NULL, sampleGroups=NULL, meanSmooth=FALSE,
                          ksmooth=1, itersmooth=1, binSize=100, lwd=0.2){
    names( sampleGroups ) <- nms
    bfl <- bamDir
    names( bfl ) <- nms
    if( !is.null( sizeFactors ) ){
        names( sizeFactors ) <- nms
    }
    start( plot_gr ) <- start( plot_gr ) - smooth - 1
    end( plot_gr ) <- end( plot_gr ) + smooth
    chr <- as.character(seqnames(plot_gr))
    df <- do.call( rbind, lapply( names(bfl), function(w){
        cov <- import( bfl[w], which=plot_gr )
        gappedRegions <- gaps(cov)
        gappedRegions <- gappedRegions[seqnames(gappedRegions) == chr]
        gappedRegions <- gappedRegions[strand(gappedRegions) == "*"]
        mcols(gappedRegions)$score <- 0
        cov <- sort(c(gappedRegions, cov))
        cov <- rep( cov$score, width(cov) )
        covRle <- Rle(cov)
        bins <- unlist(tile(plot_gr, n=width(plot_gr)/binSize))
        bins <- keepSeqlevels(bins, chr)
        mcols(bins)$cov <- viewSums(Views(covRle, IRanges( start=start(bins), end=end(bins) )))
        cov <- bins$cov/binSize
        if( !is.null( sizeFactors ) ){
            cov <- cov/sizeFactors[w]
        }else if( norm ){
            cov <- cov / mean(cov, na.rm=TRUE)
        }
        cov <- pmax( fun(cov), 0 )
        if( smooth > 0){
            cov <- stats::filter(cov, rep(1/smooth, smooth))
        }
        if( meanSmooth ){
            cov <- minfi:::.meanSmoother( cov, ksmooth, itersmooth )
        }
        df <- data.frame(
            genomicPos=rowMeans( data.frame( start( bins ), end( bins ) ) ),
            cov=cov,
            sample=w )
        df <- na.omit(df)
        df
    } ) )
    if( is.null( sampleGroups ) ){
        df$sampleGroups <- df$sample
    }else{
        df$sampleGroups <- sampleGroups[df$sample]
    }
    p4 <- ggplot( df, aes( genomicPos, cov, col=sampleGroups, group=sample ) ) +
        geom_line(lwd=lwd) +
        ylab( nm )
    if( !is.null( highlight ) ){
        highlight <- subsetByOverlaps( highlight, plot_gr )
        p4 <- p4 +
            geom_rect( data = as.data.frame( highlight ), inherit.aes=FALSE,
                       aes( xmin=start, xmax=end, ymin=-Inf, ymax=Inf ),
                       fill="gray", alpha=0.35 )
    }
    if( facet ){
        p4 <- p4 + facet_grid( sample ~ . )
    }
    p4 <- p4 + xlim( start( plot_gr ), end( plot_gr ) )
    p4
}

plotCovProfileForSample <- function(
    sampleName, protein,
    plotGr,
    bamPaths=system.file("extdata/geodata", package="ColonDNATopology"),
    suffix=".bw", ... ){
    samples <- file.path( bamPaths,
                          paste0( paste( sampleName, protein, sep="-"), suffix ) )
    stopifnot(file.exists( samples ))
    pp <- plotCoverage2( samples,
                        plot_gr=plotGr,
                        nm=protein,
                        nms=sampleName, ... )
}

plotGenesFromGRangesList <- function( grangesList, plot_gr, gpSelf=FALSE ){
    models <- subsetByOverlaps( grangesList, plot_gr, ignore.strand=TRUE )
    newModels <- lapply( seq_along(models), function(i){
        xx <- models[[i]]
        xx <- keepSeqlevels(xx, unique(as.character(seqnames(xx))) )
        gps <- gaps( xx, start=start(range(xx)), end=end(range(xx)) )
        gps <- gps[strand(gps) == unique(as.character(strand(xx)))]
        strand(gps) <- "*"
        newModel <- c(xx, gps)
        mcols(newModel)$gene_id <- names(models)[i]
        mcols(newModel)$tx_name <- mcols(newModel)$gene_id
        mcols(newModel)$exonic_part <- NULL
        mcols(newModel)$tx_id <- mcols(newModel)$gene_id
        mcols(newModel)$model <- ifelse( as.character(strand(newModel)) == "*", "gap", "exon")
        newModel
    })
    newModels <- GRangesList(newModels)
    names(newModels) <- names(models)
    bp <- autoplot( newModels,
                   plot_gr, aes(type=model), group.selfish=gpSelf )@ggplot + xlim(start(plot_gr), end(plot_gr))
    bp
}


plotHiCProfileForSample <- function( sample, res="50000", plotGR, fun=function(x){ log10(x+1)},
                                    protocol="hic", matType="raw", boundaries=NULL, ctcfSites=NULL, autoZlim=NULL, ... ){
    #data(sprintf("hic_%s", sample))
    mat <- get(sprintf("hic_%s", sample))
    chr <- as.character( seqnames(plotGR) )
    mat <- mat@resolutionNamedList[[res]][[chr]]
    mat <- as.matrix( mat )
    start <- as.numeric(rownames(mat)) + 1
    granges <- GRanges( chr, IRanges( start, start + as.numeric(res) - 1 ) )
    rownames(mat) <- seq_len(ncol(mat))
    colnames(mat) <- seq_len(ncol(mat))
    if( !is.null( autoZlim ) ){
        zlim <- quantile( fun(mat), autoZlim )
    }
    p <- matrixPlotter( fun(mat), granges, plotGR, zlim=zlim, ... )
    p# + xlim(start(plotGR), end(plotGR))
}


plotGenomicVector <- function( coordinates, vector, plotGR, highlight, iter=3, title=""){
#    start( plotGR ) <- start( plotGR ) - 200000
#    end( plotGR ) <- end( plotGR ) + 200000
    keep <- subjectHits( findOverlaps( plotGR, coordinates ) )
    dfP <- data.frame(
        genomicPos=rowMeans( data.frame( start( coordinates ), end( coordinates ) )[keep,] ),
        score=vector[keep] )
    if( iter > 0 ){
        dfP$score <- minfi:::.meanSmoother( dfP$score, iter=iter )
    }
    pp <- ggplot(dfP, aes( genomicPos, score ) ) +
        geom_line() + xlim( start(plotGR), end(plotGR) )
    if( !is.null( highlight ) ){
        highlight <- subsetByOverlaps( highlight, plotGR )
        pp <- pp +
            geom_rect( data = as.data.frame( highlight ), inherit.aes=FALSE,
                         aes( xmin=start, xmax=end, ymin=-Inf, ymax=Inf ),
                         fill="gray", alpha=0.35 )
    }
    pp + xlim( start(plotGR), end(plotGR) ) + ylab(title)
}

plotLoopsFromRanges <- function( loopsGR, plot_gr, rangeAlpha=c(0.1, 0.7), rangeSize=c(0, 2), valueVars=NULL, widthLims=NULL, sd=100 ){
    if( is.null( valueVars ) ){
        mcols( loopsGR )$value <- 1
    }else{
        mcols( loopsGR )$value <- valueVars
    }
    set.seed(sd)
    mcols( loopsGR )$hts <- runif( length( loopsGR ), 10, 500 )
    loopsGR <- subsetByOverlaps( loopsGR, plot_gr )
    if( length(loopsGR) == 0 ){
        return(NULL)
    }
    loopsGR$loopType <- "same"
    pl <- ggplot( loopsGR ) + geom_arch( aes( height = hts, size = value, alpha=value ) )
    pl <- pl@ggplot + xlim(start(plot_gr), end(plot_gr)) +
        ylab("") + theme(axis.text.y=element_blank()) #+
#        scale_colour_manual( values = brewer.pal( length( unique( mcols( archs )$loopType ) ), "Dark2") )
    pl <- pl + guides( size=FALSE, alpha=FALSE )
    if( is.null( widthLims ) ){
        widthLims <- c( 0, max( mcols( loopsGR )$value ) )
    }
    pl <- pl +
        scale_alpha_continuous( range=rangeAlpha, guide="none", limits=widthLims ) +
        scale_size_continuous( range=rangeSize, guide="none", limits=widthLims )
    pl
}

plotAllHics <- function( allSamples, plotGR, heightProp=0.38, res="40000", withLabels=TRUE,
                        highlight, colorBias=0.9, labs=NULL, ... ){
    rmtheme <- theme(legend.pos="none",
                     axis.text.x=element_blank(),
                     axis.title.x=element_blank(),
                     axis.ticks.y=element_blank(),
                     plot.margin = margin(0, .1, .1, 0, "cm"))
    allP <- bplapply( seq_len(length(allSamples)), function(x){
        y <- allSamples[x]
        p <- plotHiCProfileForSample( y, res=res, plotGR,
                                     heightProp=heightProp, #zlim=hicZlim,
                                     colorBias=colorBias,
                                     highlight=highlight, ... )
        if( x != length(allSamples) ){
            p <- p + rmtheme #+ ylab(x)
        }else{
            p <- p + xlab("Genome (Mb)") +
                scale_x_continuous(
                    labels=function( x ){
                        m <- round( x/1000000, 2 )
                        ifelse( (x %% 1000000) %in% c(0, 5e5), m, "")
                    } ) +
                theme( axis.ticks.y=element_blank(), plot.margin = margin(0, .1, .1, 0, "cm") )
        }
        if( !is.null( labs ) ){
            p <- p + ylab(labs[x])
        }
        p
    } )
    if( !withLabels ){
        allSamples <- rep( "", length( allSamples ) )
    }
    hts <- rep( 1, length( allSamples ) )
    hts[length(allSamples)] <- 1.6
    ggarrange(
        plotlist=allP, nrow=length( allSamples ),
        labels=allSamples, heights=hts, align="v")
}

plotCompartmentHeatmap <- function( coordinates, vec, plotGR, ylim=3.5, bias=1 ){
    plotGR2 <- GenomicRanges::resize( plotGR, width( plotGR ) * 1.3, fix="center" )
    keep <- queryHits( findOverlaps( coordinates, plotGR2 ) )
    dfP <- as.data.frame( coordinates[keep,] )[,c("start", "end")]
    dfP$vec <- vec[keep]
    dfP$vec <- pmin( pmax( dfP$vec, -ylim ), ylim)
    cols1 <- colorRampPalette( c( "#ffffff", "#0002aa" ) )(10)
    cols2 <- colorRampPalette( c( "#ffffff", "#f7cd46" ) )(10)
    cols1 <- rev( colorRampPalette( cols1, bias=bias )(35) )
    cols2 <- colorRampPalette( cols2, bias=bias )(35)
    cols <- rev(c(cols1, cols2))
    brks <- seq( -ylim, ylim, length.out=71 )
    dfP %>%
        ggplot( aes( xmin=start, xmax=end, ymin=0, ymax=1, fill=vec ) ) +
        geom_rect() +
        scale_fill_gradientn(
            colours=cols, limits=c(-ylim, ylim),
            breaks=brks, labels=NULL ) +
        theme( axis.text.y=element_blank(), axis.ticks.y=element_blank() ) +
        xlim( start(plotGR), end(plotGR) )
}

plotDataTracks <- function( plotGR, verbose=TRUE, heightProp=1/5,
                           hicHeight=.6, hicZlim=TRUE, samps=NULL, grp=NULL,
                           blocks, ctcfLoops, loops, res="50000", smooth=100, ... ){
    if( verbose ) cat("Plotting methylation\n")
    p3 <- plotMethBins(
        hansen,
        plotGR, extend=0, binSize=12000,
        grp=colData(hansen)$type,
        iter=1, typeMeth="raw", highlight=blocks )
    if( verbose ) cat("Plotting HiC\n")
    p2 <- plotHiCProfileForSample( samps[1], res=res, plotGR,
                                  heightProp=heightProp, autoZlim=hicZlim,#zlim=hicZlim,
                                  colorBias=0.9,
                                  highlight=blocks, ... )
    if( verbose ) cat("Plotting CTCF\n")
    p4 <- plotCovProfileForSample( sampleName=samps[2], protein="CTCF",
                                  plotGr=plotGR, smooth=smooth, sizeFactors=1,
                                  highlight=blocks, sampleGroups=grp )
    if( verbose ) cat("Plotting K27Ac\n")
    p6 <- plotCovProfileForSample( sampleName=samps[3], protein="K27Ac",
                                  plotGr=plotGR, smooth=smooth, sizeFactors=1,
                                  highlight=blocks, sampleGroups=grp )
    if( verbose ) cat("Plotting K9me3\n")
    p7 <- plotCovProfileForSample( sampleName=samps[4], protein="K9me3",
                                  plotGr=plotGR, smooth=smooth, sizeFactors=1,
                                  highlight=blocks, sampleGroups=grp )
    if( verbose ) cat("Plotting K27me3\n")
    p8 <- plotCovProfileForSample( sampleName=samps[5], protein="K27me3",
                                  plotGr=plotGR, smooth=smooth, sizeFactors=1,
                                  highlight=blocks, sampleGroups=grp )
    if( verbose ) cat("Plotting A/B ratio\n")
    p9 <- plotGenomicVector( wholeGenome, scoreMatNorm[,samps[6]], plotGR,
                            blocks, iter=1, title="A/B" )
    if( verbose ) cat("Plotting genes\n")
    pGenes <- plotGenes( plotGR )
    if( verbose ) cat("Plotting loops\n")
    if( grp == "Normal" ){
        kp <- loops@colData$groups %in% "normal"
    }else{
        kp <- !loops@colData$groups %in% "normal"
    }
    pLoops <- plotArchLoops( loops[,which(kp)], plotGR, grp="groups", facet=FALSE )
    pLoops2 <- plotLoopsFromRanges( ctcfLoops, plotGR )
    print("paso\n")
    rmtheme <- theme(legend.pos="none",
                     axis.text.x=element_blank(),
                     axis.title.x=element_blank(),
                     plot.margin = margin(0, 0, .1, 0, "cm"))
    scaleCols <- c("#762a83", "#1b7837")
    names( scaleCols ) <- c("Cancer", "Normal")
    scaleCol <- scale_color_manual(values=scaleCols)
    ggarrange( p2 + rmtheme,
              p9 + rmtheme,
              pLoops + rmtheme,
              pLoops2 + rmtheme,
              p3 + rmtheme + scaleCol,
              p4 + rmtheme + scaleCol,
              p6 + rmtheme + scaleCol,
              p7 + rmtheme + scaleCol,
              p8 + rmtheme + scaleCol,
              pGenes + theme(legend.pos="none") +
              scale_x_continuous(labels=function(x){round(x/1000000, 2)}) +
              xlab("Genome (Mb)"),
              nrow = 10, align="v",
              heights=c(hicHeight, 0.75, 0.75, 0.75, 1, 1, 1, 1, 1, 1.5) )
}

plotDataTracks2 <- function( plotGR, verbose=TRUE, heightProp=1/5,
                           hicHeight=.6, hicZlim=TRUE, samps=NULL, grp=NULL,
                           blocks, ctcfLoops, loops, res="50000", smooth=100, ... ){
    ## if( verbose ) cat("Plotting methylation\n")
    ## p3 <- plotMethBins(
    ##     hansen,
    ##     plotGR, extend=0, binSize=12000,
    ##     grp=colData(hansen)$type,
    ##     iter=1, typeMeth="raw", highlight=blocks )
    if( verbose ) cat("Plotting HiC\n")
    p2 <- plotHiCProfileForSample( samps[1], res=res, plotGR,
                                  heightProp=heightProp, autoZlim=hicZlim,#zlim=hicZlim,
                                  colorBias=0.9,
                                  highlight=blocks, ... )
    ## if( verbose ) cat("Plotting CTCF\n")
    ## p4 <- plotCovProfileForSample( sampleName=samps[2], protein="CTCF",
    ##                               plotGr=plotGR, smooth=smooth, sizeFactors=1,
    ##                               highlight=blocks, sampleGroups=grp )
    if( verbose ) cat("Plotting K27Ac\n")
    p6 <- plotCovProfileForSample( sampleName=samps[3], protein="K27Ac",
                                  plotGr=plotGR, smooth=smooth, sizeFactors=1,
                                  highlight=blocks, sampleGroups=grp )
    ## if( verbose ) cat("Plotting K9me3\n")
    ## p7 <- plotCovProfileForSample( sampleName=samps[4], protein="K9me3",
    ##                               plotGr=plotGR, smooth=smooth, sizeFactors=1,
    ##                               highlight=blocks, sampleGroups=grp )
    ## if( verbose ) cat("Plotting K27me3\n")
    ## p8 <- plotCovProfileForSample( sampleName=samps[5], protein="K27me3",
    ##                               plotGr=plotGR, smooth=smooth, sizeFactors=1,
    ##                               highlight=blocks, sampleGroups=grp )
    ## if( verbose ) cat("Plotting A/B ratio\n")
    ## p9 <- plotGenomicVector( wholeGenome, scoreMatNorm[,samps[6]], plotGR,
    ##                         blocks, iter=1, title="A/B" )
    if( verbose ) cat("Plotting genes\n")
    pGenes <- plotGenes( plotGR )
    if( verbose ) cat("Plotting loops\n")
    if( grp == "Normal" ){
        kp <- loops@colData$groups %in% "normal"
    }else{
        kp <- !loops@colData$groups %in% "normal"
    }
#    pLoops <- plotArchLoops( loops[,which(kp)], plotGR, grp="groups", facet=FALSE )
#    pLoops2 <- plotLoopsFromRanges( ctcfLoops, plotGR )
    print("paso\n")
    rmtheme <- theme(legend.pos="none",
                     axis.text.x=element_blank(),
                     axis.title.x=element_blank(),
                     plot.margin = margin(0, 0, .1, 0, "cm"))
    scaleCols <- c("#762a83", "#1b7837")
    names( scaleCols ) <- c("Cancer", "Normal")
    scaleCol <- scale_color_manual(values=scaleCols)
    ggarrange( p2 + rmtheme +
               coord_cartesian(
                   xlim=c(start(plotGR),end(plotGR)),
                   ylim=c(0, width(plotGR)*heightProp)),
##              p9 + rmtheme,
##              pLoops + rmtheme,
##              pLoops2 + rmtheme,
##              p3 + rmtheme + scaleCol,
##              p4 + rmtheme + scaleCol,
              p6 + rmtheme + scaleCol,
##              p7 + rmtheme + scaleCol,
##              p8 + rmtheme + scaleCol,
              pGenes + theme(legend.pos="none") +
              scale_x_continuous(labels=function(x){round(x/1000000, 2)}) +
              xlab("Genome (Mb)"),
              nrow = 3, align="v",
              heights=c(hicHeight, 0.75, 1.5) )
}

plotSubcompartmentHeatmap <- function( coordinates, vec, plotGR ){
    plotGR2 <- GenomicRanges::resize( plotGR, width( plotGR ) * 1.3, fix="center" )
    keep <- queryHits( findOverlaps( coordinates, plotGR2 ) )
    dfP <- as.data.frame( coordinates[keep,] )[,c("start", "end")]
    dfP$vec <- vec[keep]
    dfP %>%
        ggplot( aes( xmin=start, xmax=end, ymin=0, ymax=1, fill=vec ) ) +
        geom_rect() +
        scale_fill_manual( values=c(`A`="#0002aa", `B`="#f7cd46", `I`="#66BFE3") ) +
        theme( axis.text.y=element_blank(), axis.ticks=element_blank(), axis.line=element_blank() ) +
        xlim( start(plotGR), end(plotGR) )
}

plotChipsWithMeth <- function( samps, prot, coors, hl, normFacs, brks=c(0, 10, 20), ylims=c(0, 25) ){
    names(normFacs) <- samps
    allP <- mapply(
        function(x, y, l){
            p <- plotCovProfileForSample( x, protein=prot, coors,
                                         smooth=0, meanSmooth=TRUE,
                                         binSize=2500,
                                         ksmooth=1,
                                         itersmooth=1,
                                         sizeFactors=y,
                                         highlight=hl,
                                         bamPath="/data/aryee/bernstein/chip/SummarizedData/bam",
                                         suffix=".bam",
                                         sampleGroups=l) + #ylim(0, 25) +
                geom_area(aes(fill=sampleGroups)) +
                theme(legend.pos="none")
            p <- p +
                scale_fill_manual(values=scaleCols) +
                scale_colour_manual(values=scaleCols) +
                scale_y_continuous( breaks=brks, limits=ylims ) +
                ylab("")
            p
        }, samps, normFacs, l=ifelse(grepl("N$", samps), "Normal", "Cancer"),
    SIMPLIFY=FALSE)
    meth <- plotMethBins(
        hansen,
        coors, extend=0, binSize=5000,
        grp=colData(hansen)$type,
        iter=2, typeMeth="raw", highlight=cn2 ) +
        scale_color_manual(values=scaleCols) +
        rmtheme +
        ylab("Meth.") +
        theme(legend.pos="none") +
        scale_y_continuous(labels=function(x){round(x, 1)}, breaks=c(0.5, 1), limits=c(0.2, 1))
    pl <- ggarrange(
        meth + theme(plot.margin=unit(c(0.5, 0, 0, 0), "cm") ),
        allP[[1]] + rmtheme,
        allP[[2]] + rmtheme,
        allP[[3]] + rmtheme,
        allP[[4]] + theme(legend.pos="none") +
        scale_x_continuous(labels=function(x){round(x/1000000, 2)}) +
        theme(plot.margin = unit(c(0.2, 0, 0.5, 0), "cm")) +
        xlab("Genome (Mb)"),
        align="v", nrow=5, heights=c(1, .8, .8, .8, 1.8) )
    pl <- annotate_figure(
        pl,
        top = text_grob(sprintf("%0.1f Mb", round( width( hl )/1000000, 1)),
                        hjust = -0.2,
                        vjust=2.9,
                        size = 10),
        left=text_grob(paste0("H3", prot), size=13, rot=90, vjust=2.3 ))
    pl
}

plotMultChipsWithMeth <- function( samps, prot, coors, hl, normFacs, brks=c(0, 10, 20), ylims=c(0, 25) ){
    names(normFacs) <- samps
    allP <- mapply(
        function(x, y, l, protZ, brksZ, ylimsZ){
            p <- plotCovProfileForSample( x, protein=protZ, coors,
                                         smooth=0, meanSmooth=TRUE,
                                         binSize=2500,
                                         ksmooth=1,
                                         itersmooth=1,
                                         sizeFactors=y,
                                         highlight=hl,
                                         bamPath="/data/aryee/bernstein/chip/SummarizedData/bam",
                                         suffix=".bam",
                                         sampleGroups=l) + #ylim(0, 25) +
                geom_area(aes(fill=sampleGroups)) +
                theme(legend.pos="none")
            p <- p +
                scale_fill_manual(values=scaleCols) +
                scale_colour_manual(values=scaleCols) +
                scale_y_continuous( breaks=brksZ, limits=ylimsZ ) +
                ylab(sprintf("%7s", paste0("H3", protZ))) +
                theme( axis.title.y=element_text(angle=0, size=11, hjust=0) )
            p
        }, samps, normFacs, l=ifelse(grepl("N$", samps), "Normal", "Cancer"),
        prot, brks, ylims, SIMPLIFY=FALSE)
    meth <- plotMethBins(
        hansen,
        coors, extend=0, binSize=5000,
        grp=colData(hansen)$type,
        iter=2, typeMeth="raw", highlight=cn2 ) +
        scale_color_manual(values=scaleCols) +
        rmtheme +
        ylab("  DNAme") +
        theme(legend.pos="none") +
        scale_y_continuous(labels=function(x){round(x, 1)}, breaks=c(0.5, 1), limits=c(0.2, 1)) +
        theme( axis.title.y=element_text(angle=0, size=11, hjust=0) )
    pl <- ggarrange(
        meth + theme(plot.margin=unit(c(0.5, 0, 0, 0), "cm") ),
        allP[[1]] + rmtheme,
        allP[[2]] + rmtheme,
        allP[[3]] + rmtheme,
        allP[[4]] + rmtheme,
        allP[[5]] + rmtheme,
        allP[[6]] + theme(legend.pos="none") +
        scale_x_continuous(labels=function(x){round(x/1000000, 2)}) +
        theme(plot.margin = unit(c(0.2, 0, 0.5, 0), "cm") ) +
        xlab("Genome (Mb)"),
        align="v", nrow=7, heights=c(1, .8, .8, .8, .8, .8, 2.1) )
#    pl <- annotate_figure(
#        pl,
#        top = text_grob(sprintf("%0.1f Mb", round( width( hl )/1000000, 1)),
#                        hjust = -1.2,
#                        vjust=2.9,
#                        size = 9)#,
#        left = text_grob( "H3K27Ac   H3K27me3    H3K9m3", size=10, rot=90, vjust=2.3 )
#    )
    pl
}


plotMultChipsWithMeth2 <- function( samps, prot, coors, hl, normFacs, brks=c(0, 10, 20), ylims=c(0, 25) ){
    names(normFacs) <- samps
    allP <- mapply(
        function(x, y, l, protZ, brksZ, ylimsZ){
            p <- plotCovProfileForSample( x, protein=protZ, coors,
                                         smooth=0, meanSmooth=TRUE,
                                         binSize=2500,
                                         ksmooth=1,
                                         itersmooth=1,
                                         sizeFactors=y,
                                         highlight=hl,
                                         bamPath="/data/aryee/bernstein/chip/SummarizedData/bam",
                                         suffix=".bam",
                                         sampleGroups=l) + #ylim(0, 25) +
                geom_area(aes(fill=sampleGroups)) +
                theme(legend.pos="none")
            p <- p +
                scale_fill_manual(values=scaleCols) +
                scale_colour_manual(values=scaleCols) +
                scale_y_continuous( breaks=brksZ, limits=ylimsZ ) +
                ylab(sprintf("%7s", paste0("H3", protZ))) +
                theme( axis.title.y=element_text(angle=0, size=11, hjust=0) )
            p
        }, samps, normFacs, l=ifelse(grepl("N$", samps), "Normal", "Cancer"),
        prot, brks, ylims, SIMPLIFY=FALSE)
    meth <- plotMethBins(
        hansen,
        coors, extend=0, binSize=5000,
        grp=colData(hansen)$type,
        iter=2, typeMeth="raw", highlight=cn2 ) +
        scale_color_manual(values=scaleCols) +
        rmtheme +
        ylab("  DNAme") +
        theme(legend.pos="none") +
        scale_y_continuous(labels=function(x){round(x, 1)}, breaks=c(0.5, 1), limits=c(0.2, 1)) +
        theme( axis.title.y=element_text(angle=0, size=11, hjust=0) )
    psubcomp <- plotSubcompartmentHeatmap( wholeGenome, subCompVec, coors ) +
        rmtheme
    pl <- ggarrange(
        psubcomp,
        meth + theme(plot.margin=unit(c(0.5, 0, 0, 0), "cm") ),
        allP[[1]] + rmtheme,
        allP[[2]] + rmtheme,
        allP[[3]] + rmtheme,
        allP[[4]] + rmtheme,
        allP[[5]] + rmtheme,
        allP[[6]] + theme(legend.pos="none") +
        scale_x_continuous(labels=function(x){round(x/1000000, 2)}) +
        theme(plot.margin = unit(c(0.2, 0, 0.5, 0), "cm") ) +
        xlab("Genome (Mb)"),
        align="v", nrow=8, heights=c(.5, 1, .8, .8, .8, .8, .8, 2.1) )
#    pl <- annotate_figure(
#        pl,
#        top = text_grob(sprintf("%0.1f Mb", round( width( hl )/1000000, 1)),
#                        hjust = -1.2,
#                        vjust=2.9,
#                        size = 9)#,
#        left = text_grob( "H3K27Ac   H3K27me3    H3K9m3", size=10, rot=90, vjust=2.3 )
#    )
    pl
}
