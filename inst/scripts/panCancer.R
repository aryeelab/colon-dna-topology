
library(minfi)
source("motifInitialize.R")

require(doParallel)
registerDoParallel(cores = 5)

#####################################
### Using existing colon B blocks ###
#####################################

blocks <- readRDS("../output/rds/mergedBlocks.rds")
names(blocks) <- sprintf("block%.5d", seq_len( length( blocks ) ) )
comp <- readRDS("../output/rds/compDiffObject.rds")
blocksB <- subsetByOverlaps( blocks, reduce(comp[which(mcols( comp )$subComp == "B"),] ) )


######################
### Calling blocks ###
######################
definedBlocksForProj <- function( proj ){
    print(proj)
    meth <- readRDS( sprintf("/data/aryee/areyes/tcga/scripts/meth/SE/meth_TCGA-%s.rds", proj) )
    sampleType <- gsub( "\\S$", "", sapply( strsplit( colData( meth )$cases, "-" ), "[[", 4 ) )
    colData( meth )$type <-
                      dplyr::case_when(
                                 sampleType %in% c("01", "02", "03", "04", "09", "40") ~ "Tumor",
                                 sampleType %in% c("10", "11", "12", "14") ~ "Normal",
                                 !is.na(sampleType) ~ "Other" )
    meth <- meth[,colData( meth )$type %in% c("Tumor", "Normal")]
    if(!all( c("Normal", "Tumor") %in% colData(meth)$type )){
        return(NULL)
    }
    spBarcode <- strsplit( colData( meth )$cases, "-" )
    colData(meth)$patient <- paste( sapply( spBarcode, "[[", 1 ),
                                   sapply( spBarcode, "[[", 2 ),
                                   sapply( spBarcode, "[[", 3 ), sep="-" )
    meth <- meth[,!duplicated( colData( meth )[,c("patient", "type")] )]
    cl <- cpgCollapse( meth, what="Beta", verbose=TRUE )
    mod <- model.matrix( ~type, colData( meth ) )
    blocks <- blockFinder( cl$object, mod, what="Beta",
                          cluster=cl$blockInfo$pns, pickCutoff=TRUE, B=10 )
##    kp <- p.adjust( blocks$table$p.value, method="BH" ) < 0.1
    blocks <- blocks$table#[which(kp),]
#    blocks <- blocks[blocks$value < 0,]
    blocks_gr <- with( blocks,
                      GRanges(chr, IRanges(start,end), value=value, p.value=p.value) )
    saveRDS( blocks_gr,
            file=sprintf( "/data/aryee/areyes/GB/glioma_topology/output/blocksPan/blocks_%s.rds", proj) )
    blocks_gr
}

tcgaProjects <- gsub("meth_TCGA-|.rds", "", list.files("/data/aryee/areyes/tcga/scripts/meth/SE/"))
tcgaProjects

fl <- "/data/aryee/areyes/GB/glioma_topology/output/blocksPerCancer.rds"
if(!file.exists(fl)){
    allBlocks <- lapply( tcgaProjects, definedBlocksForProj )
    names( allBlocks ) <- tcgaProjects
    saveRDS( allBlocks, file=fl )
}
allBlocks <- readRDS( fl )

lengths(lapply( allBlocks, function(x){ if (is.null(x))return(NULL); x[width( x ) > 25000] } ))

############################
### Create block summary ###
############################
blockSum <- t(as.data.frame(lapply( allBlocks, function(bl){
    if(is.null(bl)){
        return(c(hypo=NA, hyper=NA, hypoPerc=NA, hypoMean=NA, hypoSum=NA))
    }
    bl <- bl[which(p.adjust( bl$p.value, method="BH" ) < 0.25)]
    bl <- bl[width(bl) > 25000]
    c(hypo=sum(bl$value < 0),
      hyper=sum(bl$value > 0),
      hypoPerc=100*sum(bl$value < 0)/length(bl),
      hypoMean=mean(width(bl[bl$value < 0]) ),
      hypoTot=sum(width(bl[bl$value < 0]) ) )
} )))
blockSum <- as.data.frame(blockSum)
blockSum$project <- rownames(blockSum)

hasBlocks <- blockSum$proj[which(blockSum$hypo > 50 &
                                 blockSum$hypoSum > 1000000)]


blockScoreDf <- bplapply( tcgaProjects, function(proj){
    print(proj)
    methAll <-
        readRDS( sprintf("/data/aryee/areyes/tcga/scripts/meth/SE/meth_TCGA-%s.rds",
                         proj) )
    islands <- minfi:::.getIslandAnnotation( object = methAll )
    probes <- rownames( islands )[which(islands$Relation_to_Island == "OpenSea")]
    meth <- methAll[probes,]
    ## blocks <- allBlocks[[proj]]
    ## blocks <- blocks[which(p.adjust( blocks$p.value, method="BH" ) < 0.25)]
    ## blocks <- blocks[which(blocks$value < 0)]
    ## blocks <- blocks[width(blocks) > 25000]
    blocks <- blocksB
    meth <- subsetByOverlaps( meth, blocks )
    colData(methAll)$blockScore <- colMeans( assays( meth )[["Beta"]], na.rm=TRUE )
    sampleType <- gsub( "\\S$", "",
                       sapply( strsplit( colData( methAll )$cases, "-" ), "[[", 4 ) )
    colData( methAll )$type <-
                      dplyr::case_when(
                                 sampleType %in% c("01", "02", "03", "04", "09", "40") ~ "Tumor",
                                 sampleType %in% c("10", "11", "12", "14") ~ "Normal",
                                 !is.na(sampleType) ~ "Other" )
#    methAll <- methAll[,sampleType %in% c("01", "02", "03", "04", "09", "40")]
    spName <- strsplit( colData( methAll )$cases, "-" )
    colData(methAll)$patient <-
                       paste(
                           sapply( spName, "[[", 1 ),
                           sapply( spName, "[[", 2 ),
                           sapply( spName, "[[", 3), sep="-")
    blockScoreDf <- colData(methAll)[,c("patient", "blockScore", "type")]
    blockScoreDf <- blockScoreDf[!duplicated( blockScoreDf[,c("patient", "type")] ),]
    blockScoreDf <- as.data.frame( blockScoreDf )
    blockScoreDf$proj <- proj
    blockScoreDf
}, BPPARAM=MulticoreParam(5))
names( blockScoreDf ) <- hasBlocks
blockScoreDf <- dplyr::bind_rows( blockScoreDf )

library(magrittr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

png("../output/plots/jitter_blocksB_pancancer.png", width=10, height=10, unit="in", res=300)
blockScoreDf %>%
    dplyr::filter( type %in% c("Tumor", "Normal") ) %>%
    ggplot( aes( type, blockScore, col=type ) ) +
    geom_jitter(width=0.25, size=0.4) +
    facet_wrap( ~proj ) +
    panel_border( color="black", size=1 ) +
    scale_color_manual( values=c(`Normal`="#4daf4a90", `Tumor`="#984ea390") ) +
    theme(axis.line=element_blank(), legend.pos="none") +
    labs(y="Block methylation", x="")
dev.off()

hypoTumors <- blockScoreDf %>%
    dplyr::filter( type == "Tumor" ) %>%
    dplyr::group_by( proj ) %>%
    dplyr::summarize( hypo=sum(blockScore < 0.6), tot=dplyr::n() ) %>%
    dplyr::mutate( perc=100*hypo/tot ) %>%
    dplyr::filter( perc > 20 ) %>%
    dplyr::pull( proj )

###########################################
### Adding stromal fraction information ###
###########################################

immune <- readRDS("../input/mmc2.rds")
immune <- immune[,c("TCGA.Participant.Barcode", "TCGA.Study", "Stromal.Fraction")]
colnames( immune ) <- c("patient", "proj", "stromalFraction")
blockScoreDf <- dplyr::left_join( blockScoreDf, immune )

#################################################################
### Survival analysis based only on global methylation levels ###
#################################################################
library(magrittr)
library(survival)
library(survminer)

blockScoreDf %>%
    dplyr::group_by( proj ) %>%
    dplyr::summarize( md=median( stromalFraction, na.rm=TRUE ) )


panCox <- dplyr::bind_rows( lapply( tcgaProjects, function( project ){
    print(project)
    blockDataSub <- blockScoreDf %>%
        dplyr::filter( proj == project )
    blockDataSub <- na.omit( blockDataSub )
    survDat <- read.delim(sprintf("/data/aryee/areyes/tcga/scripts/survival/%s_survival.txt.gz", project),
                          header=TRUE)
    survDat <- data.frame(
        patient=survDat$`X_PATIENT`,
        fustat=as.numeric(survDat$`OS` == 1),
        futime=survDat$`OS.time` ) 
    survDat <- survDat[!duplicated( survDat ),]    
    survDat <- dplyr::left_join( blockDataSub, survDat )
    survDat <- na.omit( survDat )
    if( nrow(survDat) == 0 ){
        return(NULL)
    }
#    survDat <- survDat[survDat$stromalFraction < 0.6,]
    survDat$blockScoreBin <- ifelse( survDat$blockScore > quantile( survDat$blockScore, 0.9 ),
                                    "high", "low" )
    survDat$blockScoreBin <- factor( survDat$blockScoreBin, levels=c("low", "high") )
    survDat$stromalFractionBin <-  ggplot2::cut_interval( survDat$stromalFraction, 4 )
    surv_object <- Surv(time = survDat$futime, event = survDat$fustat)
    fit.coxph <- coxph(
        surv_object ~ blockScoreBin,
        data = survDat )
    summ <- broom::tidy(fit.coxph)
    rs <- as.data.frame( (summ[summ$term == "blockScoreBinhigh",c("estimate", "conf.low", "conf.high")]) )
    rs$proj <- project
    rs
} ) )

panCox <- panCox[is.finite(panCox$`conf.low`),]
panCox$proj <- factor( panCox$proj )
panCox <- panCox %>%
    dplyr::mutate(
               proj=forcats::fct_reorder( proj, estimate ),
               between=ifelse( sign(conf.low)==sign(conf.high), "no", "yes" ) )

library(cowplot)
theme_set(theme_cowplot())

png("../output/plots/coxp_meth_pancancer.png",
    height=5, width=3.5, unit="in", res=300)
panCox %>%
    dplyr::filter( proj %in% hypoTumors ) %>%
    ggplot( aes( proj, estimate, col=between ) ) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
    geom_point() +
    geom_hline(yintercept=0, linetype="dashed") +
    labs(x="Tumor type", y="High methylation hazard ratio") +
    theme( legend.pos="none") +
    scale_color_manual(values=c(no="red", yes="black")) +
    coord_flip() 
dev.off()

#################################################################
### Survival analysis based only on global methylation levels ###
#################################################################

coorData <- rowRanges(readRDS("../output/expressionData_anno_stromFiltered.rds"))
hypomethGenes <- readRDS("../output/rds/candidate_IB_list.rds")

library(DESeq2)
library(edgeR)
library(limma)

resultsDEFile <- "../output/rds/panCancerDE.rds"
if(!file.exists(resultsDEFile)){
    resultsDE <- lapply( tcgaProjects, function( project ){
        print( project )
        path <- "/data/aryee/areyes/tcga/scripts/rnaseq/SE/"
        dsd <- readRDS( sprintf("%s/dsd_TCGA-%s.rds", path, project))
        dsd <- dsd[rowRanges( dsd )$gene_type == "protein_coding",]
        sampleType <- gsub( "\\S$", "",
                           sapply( strsplit( colData( dsd )$cases, "-" ), "[[", 4 ) )
        colData(dsd)$sample_type <-
                       ifelse( sampleType %in% c("01", "02", "03", "04", "09", "40"),
                              "Tumor", "Normal" )
        colData(dsd)$sample_type <-
                       factor( colData(dsd)$sample_type, levels=c("Normal", "Tumor") )
        spBarcode <- strsplit( colData( dsd )$cases, "-" )
        colData( dsd )$patient <- paste( sapply( spBarcode, "[[", 1 ),
                                sapply( spBarcode, "[[", 2 ),
                                sapply( spBarcode, "[[", 3 ), sep="-" )
        dsd <- dsd[,!duplicated( colData(dsd)[,c("patient", "sample_type")] )]
        design(dsd) <- ~ sample_type
##        blGenes <- names(coorData)[queryHits(findOverlaps( coorData, allBlocks[[project]] ))]
        blGenes <- names(coorData)[queryHits(findOverlaps( coorData, blocksB, type="within" ))]
        if( all(colData(dsd)$sample_type == "Tumor") ){
            return(NULL)
        }
        dge <- DGEList( counts( dsd ) )
        keep2 <- !rownames(dsd) %in% blGenes
        design <- model.matrix( ~ sample_type,
                               as.data.frame( colData( dsd ) ) )
        dge$samples <- calcNormFactors( dge[which(keep2),] )$samples
        normCounts <- t(t(counts(dsd))/dge$samples$norm.factors)
        v <- voom(dge, design)
        fit <- lmFit(v, design)
        fit <- eBayes(fit)
        dfRes <- data.frame(
            geneID=rownames(dsd),
            meanBaseDe=rowMeans(normCounts),
            meanNormal=rowMeans(normCounts[,colData( dsd )$sample_type == "Normal"]),
            meanTumor=rowMeans(normCounts[,colData( dsd )$sample_type == "Tumor"]),
            pvalueDe=fit$p.value[,"sample_typeTumor"],
            qvalueDe=p.adjust( fit$p.value[,"sample_typeTumor"], method="BH" ),
            project=project,
            inBlock=rownames(dsd) %in% blGenes,
            inHypometh=rownames(dsd) %in% hypomethGenes )
        dfRes
    } )
    resultsDE <- dplyr::bind_rows( resultsDE )
    saveRDS( resultsDE, file=resultsDEFile )
}
resultsDE <- readRDS( resultsDEFile )

resultsBLFile <- "../output/rds/panCancerBL.rds"

if(!file.exists( resultsBLFile )){
    resultsBL <- bplapply( tcgaProjects, function(project){
        print(project)
        path <- "/data/aryee/areyes/tcga/scripts/rnaseq/SE/"
        dsd <- readRDS( sprintf("%s/dsd_TCGA-%s.rds", path, project))
        dsd <- dsd[rowRanges( dsd )$gene_type == "protein_coding",]
        sampleType <- gsub( "\\S$", "",
                           sapply( strsplit( colData( dsd )$cases, "-" ), "[[", 4 ) )
        colData(dsd)$sample_type <-
                       ifelse( sampleType %in% c("01", "02", "03", "04", "09", "40"),
                              "Tumor", "Normal" )
        dsd <- dsd[,colData(dsd)$sample_type == "Tumor"]
        blockScoreSub <- blockScoreDf %>%
            dplyr::filter( proj == project, type == "Tumor" )
        spBarcode <- strsplit( colData( dsd )$cases, "-" )
        colData( dsd )$patient <- paste( sapply( spBarcode, "[[", 1 ),
                                        sapply( spBarcode, "[[", 2 ),
                                        sapply( spBarcode, "[[", 3 ), sep="-" )
        dsd <- dsd[,!duplicated(colData( dsd )$patient)]
        colData(dsd) <- DataFrame(
            dplyr::left_join(
                       as.data.frame(colData( dsd )[,c("patient", "sample_type")]),
                       blockScoreSub[,c("patient", "blockScore", "stromalFraction")] ) )
        dsd <- dsd[,rowSums(is.na(colData(dsd))) == 0]
        dsd <- dsd[,which(colData(dsd)$stromalFraction < 0.6)]
        colData(dsd)$stromalFractionBin <-
                       ggplot2::cut_interval( colData(dsd)$stromalFraction, 5 )
        blGenes <- names(coorData)[queryHits(findOverlaps( coorData, blocksB ))]
        if( ncol(dsd) < 30 ){
            return(NULL)
        }
        dge <- DGEList( counts( dsd ) )
        keep2 <- !rownames(dsd) %in% blGenes
        dge$samples <- calcNormFactors( dge[which(keep2),] )$samples
        normCounts <- t(t(counts(dsd))/dge$samples$norm.factors)
        design <- model.matrix( ~ blockScore,
                               as.data.frame( colData( dsd ) ) )
        v <- voom(dge, design)
        fit <- lmFit(v, design)
        fit <- eBayes(fit)
        design <- model.matrix( ~ stromalFractionBin + blockScore,
                               as.data.frame( colData( dsd ) ) )
        v <- voom(dge, design)
        fit2 <- lmFit(v, design)
        fit2 <- eBayes(fit2)
        dfRes <- data.frame(
            geneID=rownames(dsd),
            meanBaseBl=rowMeans(normCounts),
            coefBl=fit$coefficient[,"blockScore"],
            pvalueBl=fit$p.value[,"blockScore"],
            qvalueBl=p.adjust( fit$p.value[,"blockScore"], method="BH" ),
            coefBl2=fit2$coefficient[,"blockScore"],
            pvalueBl2=fit2$p.value[,"blockScore"],
            qvalueBl2=p.adjust( fit2$p.value[,"blockScore"], method="BH" ),
            project=project,
            inBlock=rownames(dsd) %in% blGenes,
            inHypometh=rownames(dsd) %in% hypomethGenes )
        dfRes
    }, BPPARAM=SerialParam())
    resultsBL <- dplyr::bind_rows( resultsBL )
    saveRDS( resultsBL, file=resultsBLFile )
 }
resultsBL <- readRDS(resultsBLFile)

resultsBL %>%
    dplyr::filter( inHypometh ) %>%
    dplyr::group_by( project, inBlock ) %>%
    dplyr::summarize( num = sum(qvalueBl2 < 0.25 & coefBl2 > 0) ) %>%
    tidyr::pivot_wider( values_from="num", names_from="inBlock" ) %>%
    as.data.frame

png("../output/plots/volcano_hypometh_pancancer.png", width=12, height=12,
    unit="in", res=300)
resultsBL %>%
    dplyr::filter( project %in% hypoTumors ) %>%
    dplyr::arrange( inHypometh ) %>%
    ggplot( aes( coefBl*-1, -log10(pvalueBl), col=inHypometh ) ) +
    geom_point( size=0.4 ) +
    xlim(-17, 17) +
    scale_color_manual(values=c(`TRUE`="red", `FALSE`="#00000050")) +
    facet_wrap( ~project, scales="free" ) +
    labs( x="Association with hypomethylation" ) +
    panel_border( size=1, color="black" ) +
    theme(axis.line=element_blank())
dev.off()

png("../output/plots/volcano_hypoblock_pancancer.png", width=12, height=12,
    unit="in", res=300)
resultsBL %>%
    dplyr::filter( project %in% hypoTumors ) %>%
    dplyr::filter( meanBaseBl > 10 ) %>%
    dplyr::arrange( inBlock ) %>%
    ggplot( aes( coefBl*-1, -log10(pvalueBl), col=inBlock ) ) +
    geom_point( size=0.4 ) +
    xlim(-17, 17) +
    scale_color_manual(values=c(`TRUE`="red", `FALSE`="#00000050")) +
    facet_wrap( ~project, scales="free" ) +
    labs( x="Association with hypomethylation" ) +
    panel_border( size=1, color="black" ) +
    theme(axis.line=element_blank())
dev.off()

tumorKeep <- read.delim("../input/tcga_tumors.txt", sep=",", header=FALSE)
colnames( tumorKeep ) <- c("project", "tumor", "filt1", "filt2")
tumorKeep <- as.character(tumorKeep$project[tumorKeep$filt2 == "YES"])
tumorKeep <- tumorKeep[!tumorKeep %in% c("SARC", "SKCM")]

png("../output/plots/boxplot_hypoblock_pancancer.png", width=6, height=8,
    unit="in", res=300)
resultsBL %>%
    dplyr::arrange( inBlock ) %>%
    dplyr::mutate( inBlock=factor(ifelse(inBlock, "in", "out"), levels=c("out", "in"))) %>%
    dplyr::filter( project %in% hypoTumors ) %>%
    dplyr::filter( meanBaseBl > 10 ) %>%
    ggplot( aes( inBlock, coefBl*-1 ) ) +
    geom_boxplot(outlier.shape=NA) +
    facet_wrap( ~project, scales="free" ) +
    coord_cartesian(y=c(-12, 12)) +
    panel_border(color="black", size=1) +
    theme(axis.line=element_blank()) +
    geom_hline(yintercept=0, col="red") +
    labs(y="Association with hypomethylation", x="Position within blocks")
dev.off()

png("../output/plots/boxplot_hypoblock_pancancer_sub.png", width=8, height=6,
    unit="in", res=300)
resultsBL %>%
    dplyr::arrange( inBlock ) %>%
    dplyr::mutate( inBlock=factor(ifelse(inBlock, "in", "out"), levels=c("out", "in"))) %>%
    dplyr::filter( project %in% hypoTumors ) %>%
    dplyr::filter( project %in% tumorKeep ) %>%
    dplyr::filter( meanBaseBl > 10 ) %>%
    ggplot( aes( inBlock, coefBl*-1 ) ) +
    geom_boxplot(outlier.shape=NA) +
    facet_wrap( ~project, scales="free", ncol=5) +
    coord_cartesian(y=c(-12, 12)) +
    panel_border(color="black", size=1) +
    theme(axis.line=element_blank()) +
    geom_hline(yintercept=0, col="red") +
    labs(y="Association with hypomethylation", x="Position within blocks")
dev.off()


normFac <- resultsDE %>%
    dplyr::filter( !inBlock ) %>%
    dplyr::group_by( project ) %>%
    dplyr::summarize( norm=median( log2(meanTumor/meanNormal), na.rm=TRUE) )

resultsDE$inBlock <- factor( ifelse( resultsDE$inBlock, "block", "non" ),
                            levels=c("non", "block") )

png("../output/plots/boxplot_log2fc_blocks_projects.png", height=10,
    width=7, unit="in", res=300)
resultsDE %>%
    dplyr::left_join( normFac ) %>%
    dplyr::filter( meanBaseDe > 10 ) %>%
    dplyr::mutate( lfc=log2( meanTumor/meanNormal )-norm ) %>%
    ggplot( aes( inBlock, lfc ) ) +
    geom_boxplot( outlier.shape=NA ) +
    facet_wrap( ~project, ncol=4 ) +
    coord_cartesian( ylim=c(-4, 4) ) +
    geom_hline( yintercept=0, col="red" ) +
    panel_border( color="black", size=1 ) +
    theme(axis.line=element_blank()) +
    labs( x="", y=bquote("Expression"~log[2]~"("*frac(Tumor, Normal)*")"))
dev.off()


### Hypometh survival ###

exprSurv <- lapply( hypoTumors, function(proj){
    print(proj)
    gns <- as.character(resultsBL %>%
                        dplyr::filter( project == proj, inBlock,##inHypometh,
                                      qvalueBl < 0.1,
                                      coefBl > 0 ) %>%
                        dplyr::pull( geneID ))
    if( length(gns) == 0  ){
        return(NULL)
    }
    path <- "/data/aryee/areyes/tcga/scripts/rnaseq/SE/"
    dsd <- readRDS( sprintf( "%s/dsd_TCGA-%s.rds", path, proj ))
    dsd <- dsd[!duplicated(rownames(dsd)),]
    sampleType <- gsub( "\\S$", "",
                       sapply( strsplit( colData( dsd )$cases, "-" ), "[[", 4 ) )
    colData(dsd)$sample_type <-
                   ifelse( sampleType %in% c("01", "02", "03", "04", "09", "40"),
                          "Tumor", "Normal" )
    colData(dsd)$sample_type <-
                   factor( colData(dsd)$sample_type, levels=c("Normal", "Tumor") )
    dsd <- dsd[,colData(dsd)$sample_type == "Tumor"]
    spBarcode <- strsplit( colData( dsd )$cases, "-" )
    colData( dsd )$patient <- paste( sapply( spBarcode, "[[", 1 ),
                                    sapply( spBarcode, "[[", 2 ),
                                    sapply( spBarcode, "[[", 3 ), sep="-" )
    dsd <- dsd[,!duplicated( colData(dsd)$patient )]
    colnames(dsd) <- colData(dsd)$patient
    survDat <-
        read.delim(sprintf("/data/aryee/areyes/tcga/scripts/survival/%s_survival.txt.gz",
                           proj), header=TRUE)
    survDat <- data.frame(
        patient=survDat$`X_PATIENT`,
        fustat=as.numeric(survDat$`OS` == 1),
        futime=survDat$`OS.time` )
    survDat <- survDat[!duplicated( survDat ),]
    dsd <- dsd[,colnames(dsd) %in% survDat$patient]
    dsd <- estimateSizeFactors( dsd )
    exprDat <- counts(dsd, normalized=TRUE)[rownames(dsd) %in% gns,,drop=FALSE]
    colnames( exprDat ) <- colData(dsd)$patient
    exprDat <- exprDat - rowMeans( exprDat ) 
    exprDat <- t(t(exprDat)/rowSds(exprDat))
    dfExpr <- data.frame(
        patient=colnames( dsd ),
        exprScore=colMedians( exprDat, na.rm=TRUE ))
    thr <- quantile( dfExpr$exprScore, 0.85 )
    dfExpr$exprBin <- ifelse( dfExpr$exprScore > thr, "high", "low" )
    dfExpr <- dfExpr[dfExpr$exprBin %in% c("high", "low"),]
    survDat <- dplyr::left_join( survDat, dfExpr )
    years=50
    idx <- survDat$futime > years*365
    survDat$futime[idx] <- years*365
    survDat$fustat[idx] <- 0
    survDat <- na.omit(survDat)
    if(nrow(survDat) == 0){
        return(NULL)
    }
    surv_object <- Surv(time = survDat$futime, event = survDat$fustat)
    survDat$exprBin <- factor( survDat$exprBin, levels=c("low", "high") )
    fit.coxph <- coxph(
        surv_object ~ exprBin,
        data = survDat )
    summ <- broom::tidy(fit.coxph)
    summ$project <- proj
    summ
} )

pp <- dplyr::bind_rows(exprSurv) %>%
    dplyr::filter( abs(estimate) < 10 ) %>%
    dplyr::mutate(
               project=forcats::fct_reorder( project, estimate ),
               between=ifelse( sign(conf.low)==sign(conf.high), "no", "yes" ) ) %>%
        ggplot( aes( project, estimate, col=between ) ) +
        geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
        geom_point() +
        geom_hline(yintercept=0, linetype="dashed") +
        labs(x="Tumor type", y="High expression hazard ratio") +
        theme( legend.pos="none") +
        scale_color_manual(values=c(no="red", yes="black")) +
    coord_flip()

png("../output/plots/coxp_expr_pancancer.png", height=4.2, width=3.5, unit="in", res=300)
print(pp)
dev.off()
        
panCox$proj <- factor( panCox$proj )
panCox <- panCox %>%
    dplyr::mutate(
               proj=forcats::fct_reorder( proj, estimate ),
               between=ifelse( sign(conf.low)==sign(conf.high), "yes", "no" ) )
panCox

####################################################
### get hazard ratio for each gene in each tumor ###
####################################################

nJobs <- 50
resourcesList <- list("queue" = "medium", "memory"="8000", "ncpus"="1", "memorylimit"="8000" )
source("~/cluster/useCluster.R")

processed <- gsub("_survPerGene.rds", "", list.files("/data/aryee/areyes/tcga/scripts/survivalPerGene"))

hasBlocks2 <- hypoTumors
hasBlocks2 <- hypoTumors[!hypoTumors %in% processed]


survPerGeneAll <- lapply( hasBlocks2, function(proj, BP){
    print(proj)
    path <- "/data/aryee/areyes/tcga/scripts/rnaseq/SE/"
    dsd <- readRDS( sprintf("%s/dsd_TCGA-%s.rds", path, proj))
    dsd <- dsd[rowRanges( dsd )$gene_type == "protein_coding",]
    sampleType <- gsub( "\\S$", "",
                       sapply( strsplit( colData( dsd )$cases, "-" ), "[[", 4 ) )
    colData(dsd)$sample_type <-
                   ifelse( sampleType %in% c("01", "02", "03", "04", "09", "40"),
                          "Tumor", "Normal" )
    dsd <- dsd[,colData(dsd)$sample_type == "Tumor"]
    blockScoreSub <- blockScoreDf %>%
        dplyr::filter( proj == proj, type == "Tumor" )    
    spBarcode <- strsplit( colData( dsd )$cases, "-" )
    colData( dsd )$patient <- paste( sapply( spBarcode, "[[", 1 ),
                                    sapply( spBarcode, "[[", 2 ),
                                    sapply( spBarcode, "[[", 3 ), sep="-" )
    dsd <- dsd[,!duplicated(colData( dsd )$patient)]
    head(  blockScoreSub[,c("patient", "blockScore", "stromalFraction")] )
    colData(dsd) <- DataFrame(
        dplyr::left_join(
                   as.data.frame(colData( dsd )[,c("patient", "sample_type")]),
                   blockScoreSub[,c("patient", "blockScore", "stromalFraction")] )
    )
    dsd <- dsd[,which(colData(dsd)$stromalFraction < .60)]
    dsd <- estimateSizeFactors(dsd)
    dfExpr <- colData(dsd)[,c("blockScore", "stromalFraction")]
    dfExpr$patient <- colData(dsd)$patient
    survDat <- read.delim(sprintf("/data/aryee/areyes/tcga/scripts/survival/%s_survival.txt.gz", proj),
                          header=TRUE)
    survDat <- data.frame(
        patient=survDat$`X_PATIENT`,
        fustat=as.numeric(survDat$`OS` == 1),
        futime=survDat$`OS.time` )
    survDat <- survDat[!duplicated( survDat ),]
    dfExpr <- dplyr::left_join( as.data.frame(dfExpr), survDat )
    dfExpr <- na.omit(dfExpr)
    if( nrow( dfExpr ) < 50 ){
        return(NULL)
    }
    colnames(dsd) <- colData(dsd)$patient
    dsd <- dsd[,dfExpr$patient]
    out <- sprintf("/data/aryee/areyes/tcga/scripts/survivalPerGene/%s_survPerGene.rds", proj)
    rs <- bplapply( rownames( dsd ), function( y, dsd, dfExpr ){
        print(y)
        library(SummarizedExperiment)
        library(DESeq2)
        library(survival)
        library(survminer)
        dfExpr$gene <- log2(counts(dsd, normalized=TRUE)[y,]+1)
        if( all(dfExpr$gene == 0) ){
            return(NULL)
        }
        dfExpr$gene <- (dfExpr$gene - mean( dfExpr$gene ))/sd( dfExpr$gene )
        dfExpr$geneBin <- ifelse(dfExpr$gene >= quantile( dfExpr$gene, 0.5 ), "high", "low" )
        dfExpr$geneBin <- factor( dfExpr$geneBin, levels=c("low", "high"))
        coxReg <- coxph(
            formula= Surv(time = dfExpr$futime, event = dfExpr$fustat) ~ stromalFraction + gene,
            data=dfExpr )
        sumCoxp <- summary( coxReg )
        nms <- colnames(sumCoxp$coef)
        rr <- c( sumCoxp$coef["gene",nms], sumCoxp$conf.int["gene",c("lower .95", "upper .95")] )
##        rr <- sumCoxp$coef["gene",c("coef", "Pr(>|z|)")]
##        coxReg <- coxph(
##            formula= Surv(time = dfExpr$futime, event = dfExpr$fustat) ~ stromalFraction + geneBin,
##            data=dfExpr )
##        summary(coxReg)$coef
        rr
    }, dsd=dsd, dfExpr=dfExpr, BPPARAM=BP )
    names( rs ) <- rownames(dsd)
    rs <- as.data.frame(do.call(rbind, rs))
    rs$geneID <- rownames(rs)
    rownames( rs ) <- NULL
    saveRDS( rs, file=out )
    rs
}, BP=bjp )

resFiles <- list.files("/data/aryee/areyes/tcga/scripts/survivalPerGene", full.names=TRUE)
names(resFiles) <- gsub("_survPerGene.rds", "", basename(resFiles))
survPerGene <- lapply( resFiles, readRDS )
names(survPerGene) <- names(resFiles)
survPerGene <- dplyr::bind_rows( survPerGene, .id="project" )

survPerGene <-
    dplyr::left_join( survPerGene,
                     resultsBL[,c("geneID", "meanBaseBl", "project", "qvalueBl",
                                  "coefBl", "qvalueBl2", "coefBl2", "inBlock")],
                     by=c("project", "geneID"))

survPerGene %>%
    dplyr::filter( project=="COAD", `Pr(>|z|)` < 0.95 ) %>%
    dplyr::mutate( hypometh=ifelse( qvalueBl < 0.1, "yes", "no") ) %>%
    dplyr::filter( geneID == "RIMKLB" )

png( "../output/plots/volcanos_hazardratio_pancancer.png", width=8, height=8, unit="in", res=300 )
survPerGene %>%
    dplyr::filter( `Pr(>|z|)` < 0.95, meanBaseBl > 20 ) %>%
    dplyr::mutate( hypometh=ifelse( qvalueBl2 < 0.1 & inBlock &
                                    coefBl2 > 0, "yes", "no") ) %>%
    dplyr::arrange( hypometh ) %>%
    dplyr::select( coef, `Pr(>|z|)`, inBlock, project, hypometh ) %>%
    na.omit() %>%
    ggplot( aes( coef, -log10(`Pr(>|z|)`), col=hypometh ) ) +
    geom_point(alpha=0.2, size=0.5) +
    scale_color_manual(values=c(yes="red", no="black")) +
    facet_wrap( ~project, scales="free") +
    labs(x="Hazard ratio coefficient", y="p-value (survival)")
dev.off()

png( "../output/plots/boxplot_hazardratio_pancancer.png", width=6, height=8, unit="in", res=300 )
survPerGene %>%
    dplyr::filter( `Pr(>|z|)` < 0.95, meanBaseBl > 10 ) %>%
    dplyr::mutate( hypometh=ifelse( qvalueBl2 < 0.1 & inBlock &
                                    coefBl2 > 0, "yes", "no") ) %>%
    dplyr::select( hypometh, `Pr(>|z|)`, coef, project ) %>%
    na.omit() %>%
    ggplot( aes( hypometh, coef ) ) +
    geom_boxplot( outlier.shape=NA ) +
    facet_wrap( ~project ) +
    coord_cartesian(ylim=c(-1, 1)) +
    panel_border(colour="black", size=1) +
    theme(axis.line=element_blank()) +
    geom_hline(yintercept=0, col="red") +
    labs(y="Hazard ratio coefficient", x="Inside blocks and downregulated")
dev.off()


png( "../output/plots/boxplot_hazardratio_pancancer_epi.png", width=6, height=8, unit="in", res=300 )
survPerGene %>%
    dplyr::filter( `Pr(>|z|)` < 0.95, meanBaseBl > 10, project %in% tumorKeep ) %>%
    dplyr::mutate( hypometh=ifelse( qvalueBl2 < 0.1 & inBlock &
                                    coefBl2 > 0, "yes", "no") ) %>%
    dplyr::select( hypometh, `Pr(>|z|)`, coef, project ) %>%
    na.omit() %>%
    ggplot( aes( hypometh, coef ) ) +
    geom_boxplot( outlier.shape=NA ) +
    facet_wrap( ~project ) +
    coord_cartesian(ylim=c(-1, 1)) +
    panel_border(colour="black", size=1) +
    theme(axis.line=element_blank()) +
    geom_hline(yintercept=0, col="red") +
    labs(y="Hazard ratio coefficient", x="Inside blocks and downregulated")
dev.off()

gns <- survPerGene %>%
    dplyr::mutate( hypometh=ifelse( qvalueBl < 0.5 & inBlock & meanBaseBl > 20 &
                                    coefBl > 0, "yes", "no") ) %>%
    dplyr::group_by( geneID ) %>%
    dplyr::summarize( hypoTumors=sum( hypometh == "yes" ), tot=dplyr::n() ) %>%
    dplyr::filter( hypoTumors >= 12, tot == 16 ) %>%
    as.data.frame %>%
    dplyr::pull( geneID )

oncogenes <- read.delim("../input/oncogenes/ongene_human.txt")
oncogenes <- oncogenes[oncogenes$GeneType == "protein-coding",]

sum( tumorKeep %in% survPerGene$project )

tumorKeep

head(resultsBL)

gns2 <- resultsBL %>%
    dplyr::filter( project %in% hypoTumors ) %>%
    dplyr::filter( project %in% tumorKeep ) %>%
    dplyr::mutate( hypometh=ifelse( qvalueBl2 < 0.5 & inBlock & meanBaseBl > 20 &
                                    coefBl > 0, "yes", "no") ) %>%
    dplyr::group_by( geneID ) %>%
    dplyr::summarize( hypoTumors=sum( hypometh == "yes" ), tot=dplyr::n() ) %>%
##    dplyr::filter( hypoTumors >= 6, tot==7 ) %>%    
    dplyr::filter( hypoTumors >= 7, tot==10 ) %>%
    as.data.frame %>%
    dplyr::pull( geneID )
gnsBack <- unique(survPerGene %>%
    dplyr::filter( meanBaseBl > 10 ) %>%
    dplyr::pull( geneID ))
back <- table( unique(gnsBack) %in% as.vector(oncogenes$OncogeneName))
fore <- table(gns2 %in% as.vector(oncogenes$OncogeneName))
mt <- rbind( fore, back )[,2:1]
mt["back",] <- mt["back",] - mt["fore",]
fisher.test( mt )

head( resultsBL )

gns2 <- resultsBL %>%
    dplyr::filter( project %in% hypoTumors ) %>%
    dplyr::filter( project %in% tumorKeep ) %>%
    dplyr::mutate( hypometh=ifelse( qvalueBl2 < 0.5 & inBlock & meanBaseBl > 20 &
                                    coefBl > 0, "yes", "no") ) %>%
    dplyr::group_by( geneID ) %>%
    dplyr::summarize( hypoTumors=sum( hypometh == "yes" ), tot=dplyr::n() ) %>%
##    dplyr::filter( hypoTumors >= 6, tot==7 ) %>%    
    dplyr::filter( hypoTumors >= 7, tot==10 ) %>%
    as.data.frame %>%
    dplyr::pull( geneID )
gnsBack <- unique(survPerGene %>%
    dplyr::filter( meanBaseBl > 10 ) %>%
    dplyr::pull( geneID ))
back <- table( unique(gnsBack) %in% as.vector(oncogenes$OncogeneName))
fore <- table(gns2 %in% as.vector(oncogenes$OncogeneName))
mt <- rbind( fore, back )[,2:1]
mt["back",] <- mt["back",] - mt["fore",]
fisher.test( mt )


resultsBL %>%
    dplyr::filter( project %in% hypoTumors ) %>%
    dplyr::filter( project %in% tumorKeep ) %>%    
    dplyr::filter( geneID %in% gns2 ) %>%
    dplyr::filter( qvalueBl2 < 0.5 ) %>%
    dplyr::group_by( geneID ) %>%
    dplyr::summarize( number_of_tumors=dplyr::n(), tumors=paste(as.character(project), collapse=",") ) %>%
    dplyr::mutate( is_oncogene=ifelse(geneID %in% as.vector(oncogenes$OncogeneName), "yes", "no" )) %>%
    as.data.frame %>%
    write.table(file="../output/pancancer_genelist.txt", quote=FALSE, sep="\t", row.names=FALSE)


length(gns2)
gns2[gns2 %in% as.vector(oncogenes$OncogeneName)]


head(gns2)

data.frame( `gene name`=gns, 

mt

candidateGenes <- readRDS("../output/rds/candidate_IB_list.rds")

intersect(intersect( gns2, oncogenes$OncogeneName ), candidateGenes)

intersect( gns2, candidateGenes)


prop.table( mt, 1 )

library(UpSetR)

oncoTable <- survPerGene %>%
    dplyr::filter( project %in% tumorKeep ) %>%
    dplyr::mutate( hypometh=ifelse( qvalueBl < 0.5 & inBlock & meanBaseBl > 20 &
                                    coefBl > 0, "yes", "no") ) %>%
    dplyr::filter( hypometh == "yes" ) %>%
    dplyr::filter( geneID %in% oncogenes$OncogeneName ) %>%
    dplyr::select( project, geneID )
oncoTable <- split( oncoTable$geneID, oncoTable$project )

png("prueba.png")
upset( fromList(oncoTable), nsets=100, nintersects=100 )
dev.off()


png( "../output/plots/boxplot_hazardratio_allHyp_pancancer.png", width=6, height=8, unit="in", res=300 )
survPerGene %>%
    dplyr::filter( `Pr(>|z|)` < 0.95, project %in% tumorKeep ) %>%
    dplyr::mutate( hypometh=ifelse( geneID %in% gns, "yes", "no" ) ) %>%
    dplyr::select( hypometh, `Pr(>|z|)`, coef, project ) %>%
    na.omit() %>%
    ggplot( aes( hypometh, coef ) ) +
    geom_boxplot( outlier.shape=NA ) +
    facet_wrap( ~project ) +
    coord_cartesian(ylim=c(-1, 1)) +
    panel_border(colour="black", size=1) +
    theme(axis.line=element_blank()) +
    geom_hline(yintercept=0, col="red") +
    labs(y="Hazard ratio coefficient", x="Inside blocks and downregulated")
dev.off()

png( "../output/plots/volcanos_hazardratio_allHyp_pancancer.png", width=8, height=8, unit="in", res=300 )
survPerGene %>%
    dplyr::filter( `Pr(>|z|)` < 0.95 ) %>%
    dplyr::mutate( hypometh=ifelse( geneID %in% gns, "yes", "no") ) %>%
    dplyr::arrange( hypometh ) %>%
    dplyr::select( coef, `Pr(>|z|)`, inBlock, project, hypometh ) %>%
    na.omit() %>%
    ggplot( aes( coef, -log10(`Pr(>|z|)`), col=hypometh ) ) +
    geom_point(alpha=0.2, size=0.5) +
    scale_color_manual(values=c(yes="red", no="black")) +
    facet_wrap( ~project, scales="free") +
    labs(x="Hazard ratio coefficient", y="p-value (survival)")
dev.off()

survPerGene %>%
    dplyr::group_by( project ) %>%
    dplyr::mutate( padjSurv=p.adjust( `Pr(>|z|)`, method="BH") ) %>%
    as.data.frame %>%
    dplyr::mutate( signSurv=padjSurv < 0.5  ) %>%
    dplyr::filter( geneID %in% gns, signSurv ) %>%
    dplyr::group_by( project ) %>%
    dplyr::summarize( sum(coef > 0), dplyr::n() )


survPerGene %>%
    dplyr::filter( geneID %in% "BCL11B" )
