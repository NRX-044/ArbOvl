#' Title performing arbitrary overlap calculation among a set of genomic ranges
#' @description provide various options to define overlapping
#' @param xset a list of GenomicRanges object
#' @param min.gapwidth minimum gap to merge hits, default is 0
#' @param mode  decide counting overlapping intervals or overlapping bases as
#' intersection size, default is interval
#' @param ...
#'
#' @importFrom 'gtools' 'permutations'
#' @importFrom 'IRanges' 'findOverlapPairs'
#' @importFrom 'IRanges' 'pintersect'
#' @importFrom 'IRanges' 'reduce'
#' @importFrom 'IRanges' 'setdiff'
#' @return An overlap object containing a overlap matrix and a list containing
#' corresponding overlapping intervals
#' @examples
#'  \dontrun{
#'   gr1 <- IRanges(c(15,6,44),c(20,8,50))
#'   gr2 <- IRanges(c(8,18,20),c(16,22,24))
#'
#'   #running under default parameters
#'   grlist <- list(gr1,gr2)
#'   res <- arbtr_ol(grlist)
#'   res[[0]]
#'   res[[1]]
#'
#'   #alternating rules defining intersection
#'   res <- arbtr_ol(grlist,minoverlap = 3L)
#'   res[[0]]
#'   res[[1]]
#'
#'   #merging overlapping regions with distance less than n bases counts
#'   res <- arbtr_ol(grlist,min.gapwidth = 3L)
#'   res[[0]]
#'   res[[1]]
#'
#'   #counting overlapping bases instead of number of overlapping intervals,
#'   #which affects the calculation of intersection sizes
#'   res <- arbtr_ol(grlist,mode = 'base')
#'   res[[0]]
#'   res[[1]]
#'  }
#'
#' @export
#'
arbtr_ol <- function(xset,min.gapwidth = 0L,
                     mode = c('interval','base'),...){
  wp <- t(apply(permutations(2,
                             length(xset),
                             0:1,
                             repeats.allowed = TRUE),1,rev))
  rs <- rowSums(wp)
  wpdf <- as.data.frame(wp)
  wpdf$overlapnum <- NA
  inrlist <- lapply(1:nrow(wpdf),function(x){GRanges(seqnames = NULL,
                                                     ranges = NULL,
                                                     strand = NULL)})
  mode <- match.arg(mode)

  for(i in length(xset):0){
    ridx <- which(rs==i)
    if (i==0){
      wpdf$overlapnum[ridx] <- 0
    }else if(i==length(xset)){
      res <- get_setol(xset,mode,...)
      inrlist[[ridx]] <- reduce(res[[2]],min.gapwidth = min.gapwidth)
      if(mode=='interval'){
        wpdf$overlapnum[ridx] <- length(inrlist[[ridx]])
      }else if(mode=='base'){
        wpdf$overlapnum[ridx] <- sum(width(inrlist[[ridx]]))
      }
    }else{
      for(ii in ridx){
        subxset <- xset[as.logical(wp[ii,])]
        res <- get_setol(subxset,mode,...)
        olnum <- res[[1]]
        olinr <- res[[2]]
        # subtract the overlap part to get the exclusive overlap number

        # get index of zero in each combination
        tmp_wp <- wp[ii,]
        zidx <- which(wp[ii,]==0)
        zp <- permutations(2,
                           length(zidx),
                           0:1,
                           repeats.allowed = TRUE)
        for(iii in 2:nrow(zp)){
          # get the index of set which is not presented in each combination
          zset <- zidx[as.logical(zp[iii,])]
          tmp_wp[zset] <- tmp_wp[zset] + 1
          idx_mat <- t(apply(wp,1,function(x) x==tmp_wp))
          upperinr <- inrlist[[which(rowSums(idx_mat)==length(tmp_wp))]]
          olinr <- setdiff(olinr,upperinr)
          if(mode=='interval'){
            olnum <- length(olinr)
          }else if(mode=='base'){
            olnum <- sum(width(olinr))
          }
          tmp_wp <- wp[ii,]
        }

        reduced_olinr <- reduce(olinr,min.gapwidth=min.gapwidth)
        if(mode=='interval'){
          olnum <-length(reduced_olinr)
        }else if(mode=='base'){
          olnum <- sum(width(reduced_olinr))
        }
        wpdf$overlapnum[ii] <- olnum
        inrlist[[ii]] <- reduced_olinr
      }
    }
  }
  res <- list('oldf' = wpdf,
              'overlap.interval' = inrlist)
  return(res)
}


# wrapper for running findOVerlapPairs and pintersect
get_setol <- function(group,mode,...){
  init_ol <- group[[1]]
  if(length(group)==1){
    if(mode=='interval'){
      overlap_num <- length(group[[1]])
    }else if(mode=='base'){
      overlap_num <- sum(width(group[[1]]))
    }
    overlap_merge <- group[[1]]
    return(list(overlap_num,overlap_merge))
  }
  for(i in 2:length(group)){
    init_ol <- pintersect(findOverlapPairs(init_ol,group[[i]],...))
    if(mode=='interval'){
      overlap_num <- length(init_ol)
    }else if(mode=='base'){
      overlap_num <- sum(width(init_ol))
    }
  }
  return(list(overlap_num,init_ol))
}
