
#' @export
fixLabels <- function( x, labels ){
  ifelse( grepl("MGH|BRD", x ), labels[x], x )
}

