
#' Title visualizing multiple intersection by venn plot
#'
#' @param oldf overlap matrix demonstrating overlapping counts of
#' each combination, can be provided by other methods
#' @param outputfig  file path for saving figure, if one don't want to
#' save it, keep it as default
#' @param method options of visualization. Currently gplots and ggvenn are
#' supported
#' @param ...
#'
#' @importFrom 'gplots' 'plot.venn'
#' @importFrom 'ggvenn' 'ggvenn'
#' @return plotting object
#' @export
#'
olm2venn <- function(oldf,outputfig = NULL,method = c('gplots','ggvenn'),...){
  method <- match.arg(method)
  if(method=='gplots'){
    rownames(oldf) <- apply(oldf,1,function(x){
      paste(x[1:(ncol(oldf)-1)],sep='',collapse='')
    })
    olmat <- as.matrix(oldf)
    olmat <- olmat[,c(ncol(olmat),1:ncol(olmat)-1)]
    class(olmat) <- 'venn'
    if(!is.null(outputfig)){
      pdf(file=outputfig,bg='transparent')
      plot.venn(olmat,...)
      dev.off()
    }else{
      plot.venn(olmat,...)
    }
  }else if(method=='ggvenn'){
    tmp_oldf <- as.data.frame(lapply(oldf,rep,oldf$overlapnum))
    tmp_oldf <- tmp_oldf[,-ncol(tmp_oldf)]
    tmp_oldf <- as.data.frame(apply(tmp_oldf,2,function(x){
      as.logical(x)
    }))
    val_col <- seq(1:nrow(tmp_oldf))
    tmp_oldf <- cbind(val_col,tmp_oldf)
    if(!is.null(outputfig)){
      plotobj <- ggvenn(tmp_oldf,...)
      pdf(file=outputfig,bg='transparent')
      print(plotobj)
      dev.off()
    }else{
      ggvenn(tmp_oldf,...)
    }
  }
}

#' Title visualize multiple intersection by upset plot
#'
#' @param oldf overlap matrix demonstrating overlapping counts of
#' each combinations, can be provided by other methods
#' @param outputfig file path for saving figure, if one don't want to
#' save it, keep it as default
#' @param method options of visualization. Currently only UpSetR is supported
#' @param ...
#'
#' @importFrom 'UpSetR' 'upset'
#' @return plotting object
#' @export
#'
olm2upset <- function(oldf,outputfig = NULL,method = 'UpSetR',...){
  if(method=='UpSetR'){
    bin_df <- oldf2cmat(oldf)
    if(!is.null(outputfig)){
      plotobj <- upset(bin_df)
      pdf(file = outputfig,bg='transparent',onefile=FALSE)
      print(plotobj)
      dev.off()
    }else{
      upset(bin_df)
    }
  }
}

#converting overlap matrix to combination matrix that suits the requirement
# of upset plot
oldf2cmat <- function(oldf){
  tmp_oldf <- oldf
  tmp_oldf <- as.data.frame(lapply(tmp_oldf,rep,oldf$overlapnum))
  tmp_oldf <- tmp_oldf[,1:ncol(oldf)-1]
  return(tmp_oldf)
}

