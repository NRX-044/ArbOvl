
#' Title transform bed-format file into genomic ranges object
#'
#' @param file can be a string representing file path or a data frame object
#'
#' @import 'IRanges'
#' @importFrom 'GenomicRanges' 'GRanges'
#' @return a genomic ranges object
#' @export
#'
setGeneric("bed2gr",
           signature = c("file"),
           function(file){
             standardGeneric("bed2gr")
})

setMethod("bed2gr",signature("character"),function(file){
  bedtable <- read.table(file,header = FALSE)
  if(ncol(bedtable) >= 6){
    gr <- GRanges(seqnames = bedtable[,1],
                  ranges = IRanges::IRanges(bedtable[,2]+1,bedtable[,3]),
                  strand = bedtable[,6])
  }else{
    gr <- GRanges(seqnames = bedtable[,1],
                  ranges = IRanges::IRanges(bedtable[,2]+1,bedtable[,3]))
  }
  return(gr)
}
)

setMethod("bed2gr",signature("data.frame"),function(file){
  if(ncol(file) >= 6){
    gr <- GRanges(seqnames = file[,1],
                  ranges = IRanges::IRanges(file[,2]+1,file[,3]),
                  strand = file[,6])
  }else{
    gr <- GRanges(seqnames = file[,1],
                  ranges = IRanges::IRanges(file[,2]+1,file[,3]))
  }
  return(gr)
}
)


#' Title transform genomic ranges object into data frame which can be directly
#' saved into bed format file
#' @param gr single genomic range object or genomic range list object or just
#' a list
#'
#' @return  a data frame
#' @export
#'
setGeneric("gr2df",
           signature = c("gr"),
           function(gr){
  standardGeneric("gr2df")
})

setMethod("gr2df",signature("IntegerRanges"),function(gr){
  df <- data.frame(seqnames(gr),
                   start(gr)-1,
                   end(gr),
                   rep('-',length(seqnames(gr))),
                   rep('-',length(seqnames(gr))),
                   strand(gr))
  names(df) <- c()
  return(df)
})

grl2df <- function(gr){
  df <- as.data.frame(c())
  lapply(gr,function(x){
    tmp_df <- data.frame(seqnames(x),
                         start(x)-1,
                         end(x),
                         rep('-',length(seqnames(x))),
                         rep('-',length(seqnames(x))),
                         strand(x))
    df <<- rbind(df,tmp_df)
  })
  names(df) <- c()
  return(df)
}

setMethod("gr2df",signature("IntegerRangesList"),
          grl2df
)

setMethod("gr2df",signature("list"),
          grl2df
)



