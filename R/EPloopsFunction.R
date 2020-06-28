
library(GenomicRanges)

keepEPnew <- function(lto, enhancer, promoter){
  lto.df <- diffloop::summary(lto)
  Ranchors <- GRanges(lto.df$chr_1, IRanges(lto.df$start_1, lto.df$end_1))
  Lanchors <- GRanges(lto.df$chr_2, IRanges(lto.df$start_2, lto.df$end_2))

  # Determine if right anchor is near promoter region
  Rhits.p <- suppressWarnings(findOverlaps(promoter, Ranchors,
                                           maxgap = 0))
  Rvalues.p <- rep(FALSE, dim(lto.df)[1])
  Rvalues.p[unique(subjectHits(Rhits.p))] <- TRUE

  # Determine if left anchor is near promoter region
  Lhits.p <- suppressWarnings(findOverlaps(promoter, Lanchors,
                                           maxgap = 0))
  Lvalues.p <- rep(FALSE, dim(lto.df)[1])
  Lvalues.p[unique(subjectHits(Lhits.p))] <- TRUE

  # Aggregate TSS
  Rtss <- data.frame(Rhits.p)
  Rtss <- cbind(Rtss, mcols(promoter[Rtss$queryHits]))
  Ltss <- data.frame(Lhits.p)
  Ltss <- cbind(Ltss, mcols(promoter[Ltss$queryHits]))
  ttss <- unique(rbind(Rtss, Ltss)[,c(2,3)])
  tss <- aggregate(gene~subjectHits,paste,collapse=",",data=ttss)

  #######

  # Determine if right anchor is near enhancer peak
  Rhits.e <- suppressWarnings(findOverlaps(enhancer, Ranchors, maxgap = 0))
  Rvalues.e <- rep(FALSE, dim(lto.df)[1])
  Rvalues.e[unique(subjectHits(Rhits.e))] <- TRUE

  # Determine if left anchor is near enhancer peak
  Lhits.e <- suppressWarnings(findOverlaps(enhancer, Lanchors,
                                           maxgap = 0))
  Lvalues.e <- rep(FALSE, dim(lto.df)[1])
  Lvalues.e[unique(subjectHits(Lhits.e))] <- TRUE

  #######

  #Add annotation and subset
  ep.loops <- (Lvalues.e & Rvalues.p) | (Lvalues.p & Rvalues.e)
  gene.tss <- rep("none", dim(lto)[2])

  gene.tss[tss$subjectHits] <- tss$gene
  lto@rowData$loop.type <- "e-p"
  lto@rowData$gene.tss <- gene.tss
  new.loops <- subsetLoops(lto, ep.loops)
  return(new.loops)
#  lto
}
