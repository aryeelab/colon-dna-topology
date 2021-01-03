
library(GenomicDataCommons)
library(magrittr)
library(TCGAbiolinks)
library(curatedTCGAData)
library(MultiAssayExperiment)
library(TCGAutils)

pquery <- projects() %>% results(size=100000)
tcgaProjects <- names(pquery$disease_type)[grepl("TCGA", names(pquery$disease_type))]

query <- GDCquery(
    project = tcgaProjects,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "HTSeq - Counts")

gdc_set_cache( file.path( getwd(), "cache" ) )

gdcdata( query$results[[1]]$id )

for(i in unique(query$results[[1]]$project)){
    dir.create( file.path("rnaseq", i) )
}

for( i in seq_len( nrow( query$results[[1]] ) ) ){
    fl <- file.path( "cache", query$results[[1]]$id[i] )
    system(sprintf( "mv %s %s", fl, file.path("rnaseq", query$results[[1]]$project[i], "/" ) ))
}

query_meth <- GDCquery( project = tcgaProjects,
                  data.category = "DNA Methylation",
                  data.type = "Methylation Beta Value",
                  platform = "Illumina Human Methylation 450" )

library(BiocParallel)

bplapply( query_meth$results[[1]]$id, gdcdata,
         BPPARAM=MulticoreParam(5, tasks=length(query_meth$results[[1]]$id) ))


for(i in unique(query_meth$results[[1]]$project)){
    dir.create( file.path("meth", i), recursive=TRUE )
}

for( i in seq_len( nrow( query_meth$results[[1]] ) ) ){
    outFile <- file.path( "meth", query_meth$results[[1]]$project[i], query_meth$results[[1]]$file_name[i]  )
    inFile <- file.path( "cache", query_meth$results[[1]]$file_id[i], query_meth$results[[1]]$file_name[i] )
    system( sprintf( "mv %s %s", inFile, outFile ) )
}

library(minfi)

for( xx in tcgaProjects ){
    print(xx)
    querySub <- query_meth$results[[1]][query_meth$results[[1]]$project %in% xx,]
    methFiles <- file.path("meth", querySub$project, querySub$file_name )
    stopifnot(all( file.exists( methFiles ) ))
    nms <- as.data.frame(data.table::fread( methFiles[[1]] ) )[,1]
    methValsAll <- bplapply( methFiles, function(x){
        methVals <- data.table::fread( x )
        methVals <- as.data.frame( methVals[,c("Composite Element REF", "Beta_value")] )
        colnames( methVals ) <- c( "cg_name", "beta" )
        stopifnot( all( methVals$cg_name == nms ) )
        methVals
    }, BPPARAM=MulticoreParam(5) )
    betas <- as.data.frame( dplyr::bind_cols(lapply( methValsAll, function(x){x[,2]} )) )
    rownames(betas) <- nms
    colnames(betas) <- querySub$id
    tcgaMeth <- RatioSet( as.matrix(betas) )
    colData( tcgaMeth ) <- DataFrame(querySub)
    annotation( tcgaMeth ) <- c(
        array="IlluminaHumanMethylation450k",
        annotation="ilmn12.hg19")
    tcgaMeth <- mapToGenome( tcgaMeth )
    seqlevelsStyle( tcgaMeth ) <- "UCSC"
    saveRDS( tcgaMeth, file=file.path("meth", "SE", sprintf("meth_%s.rds", xx) ) )
}

for( i in gsub("TCGA-", "", tcgaProjects) ){
    download.file(
        sprintf("https://tcga.xenahubs.net/download/survival/%s_survival.txt.gz", i),
        destfile=sprintf("survival/%s_survival.txt.gz", i) )
}

