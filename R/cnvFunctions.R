#' @export
doSegmentation <- function( sampName, compBins, covData, gcContent ){
  print( sampName )
  df <- data.frame( Chr=seqnames(compBins),
                    Pos=start(compBins),
                    as.data.frame( covData )[,c( sampName, "BRD3170N" )] )
  colnames( df )[3:4] <- c("Test", "Norm")
  df$GC <- gcContent
  CN <- dataFrame2object( df )
  CN <- gcNorm( CN )
  CN <- addSmooth( CN, lambda=7 )
  CN <- peakPloidy( CN , method='closest' )
  CN <- validation( CN )
  CN <- addDNACopy( CN )
  CN <- discreteNorm( CN )
  list(
    tumorContent=as.numeric( CN@Res@suggested.tumContent ),
    segmentation=segMean.n( CN ) )
}

#' @importFrom BiocParallel bplapply
#' @importFrom biovizBase GCcontent
#' @import CNAnorm
#' @export
copyNumberForRes <- function( meth, res=40000, onlyTumors=TRUE, BPPARAM=SerialParam() ){
  compBins <- suppressWarnings(getCompBinning( res=res, paste0( "chr", c( 1:22, "X") ) ))
  names( compBins ) <- sprintf("bin%0.6d", seq_len(length(compBins)))
  ## covData <- getCoverage( methAll, compBins, type="Cov", what="perRegionTotal" )
  covData <- getSumCov( meth, compBins )
  covData <- covData[,colnames(covData) %in% colnames(meth)]
  #   covData <- covData[,!(colSums(covData, na.rm=TRUE) < 1000000)]
  gcContent <- GCcontent( BSgenome.Hsapiens.UCSC.hg19, compBins )[,1]
  if( onlyTumors ){
    tumorSamples <- colnames(meth)[!grepl( "N$", rownames( colData( meth ) ) )]
  }else{
    tumorSamples <- colnames(meth)
  }
  cnaNormRes <- bplapply( tumorSamples,
                          function(x){ try( doSegmentation(x, compBins, covData, gcContent ) ) },
                          BPPARAM=BPPARAM )
  names( cnaNormRes ) <- tumorSamples
  keep <- sapply( cnaNormRes, function(x) !inherits(x, "try-error") )
  cnaNormRes <- cnaNormRes[keep]
  mcols( compBins ) <- sapply( cnaNormRes, "[[", "segmentation" )
  colnames( mcols( compBins ) ) <- names(cnaNormRes)
  compBins
}

#' @export
getCompBinning <- function(res=150000, chrs=paste0("chr", 1:22) ){
  require(BSgenome.Hsapiens.UCSC.hg19)
  seqLengths <- seqlengths( BSgenome.Hsapiens.UCSC.hg19 )
  compartmentGR <- Reduce(c, sapply( chrs, function(chr){
    st <- seq( 1, seqLengths[chr], res )
    end <- st + res - 1
    end <- pmin(end, seqLengths[chr])
    GRanges( chr, IRanges( st, end ) )
  }) )
  mcols(compartmentGR) <- NULL
  compartmentGR
}
