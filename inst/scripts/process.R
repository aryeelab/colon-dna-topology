
library(GenomicDataCommons)
library(magrittr)
library(TCGAbiolinks)
##library(DESeq2)
library(minfi)
library(data.table)

pquery = projects() %>% GenomicDataCommons::results(size=100000)

tcgaProjects <- names(pquery$disease_type)[grepl("TCGA", names(pquery$disease_type))]

sort( tcgaProjects )

query <- GDCquery(
    project = tcgaProjects,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "HTSeq - Counts")

querySub <- query$results[[1]]

geneInfo <- read.delim("rnaseq/gencode.gene.info.v22.tsv")

grGeneInfo <- GRanges( geneInfo$seqname, IRanges( geneInfo$start, geneInfo$end ), geneInfo$strand )
names( grGeneInfo ) <- geneInfo$`gene_id`

mcols( grGeneInfo ) <- DataFrame(geneInfo[,c("gene_name", "gene_type", "gene_status", "havana_gene", "full_length", "exon_length", "exon_num")] )


rnaFiles <- file.path( "rnaseq", querySub$project, querySub$id, querySub$file_name )

for( proj in tcgaProjects ){
    print( proj )
    keep <- querySub$project == proj
    rnaSubFiles <- rnaFiles[which(keep)]
    countDat <- lapply( rnaSubFiles, function(x){
        read.delim( x, header=FALSE )
    } )
    gnNames <- lapply( countDat, function(x) as.vector(x$V1)  )
    countDat <- lapply( countDat, function(x) as.vector(x$V2)  )
    all( apply( dplyr::bind_cols( gnNames ), 1, function(x){length(unique(x))} ) == 1 )
    countDat <- as.data.frame( dplyr::bind_cols( countDat ) )
    rownames( countDat ) <- gnNames[[1]]
    colnames( countDat ) <- basename(dirname(rnaSubFiles))
    stopifnot(all(querySub$id[keep] == colnames( countDat )))
    dsd <- DESeqDataSetFromMatrix( countDat, querySub[keep,], design=~1 )
    dsd <- dsd[!grepl("__", rownames( dsd )),]
    rowRanges( dsd ) <- grGeneInfo[rownames(dsd)]
    rownames(dsd) <- mcols( dsd )$`gene_name`
    saveRDS( dsd, file=file.path("rnaseq", "SE", sprintf("dsd_%s.rds", proj)) )
}


methAll <- readRDS("/data/aryee/areyes/GB/glioma_topology/output/rds/tcgaMethAnnoAllCR.rds")


