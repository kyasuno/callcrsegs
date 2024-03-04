#' @title Find the nearest copy number state
#'
#' @description This function finds the nearest copy number state in R-HDS plane.
#' Original function is findBestXYP in BubbleTree.
#'
#' @param lrr a numeric log R ratio
#' @param hds a numeric Heterozygous deviation score
#' @param cfg, configuration parameters
#' @param p a numeric vector. If not NULL, distance to the branches are calculated only for given p's.
#' @return A tbl
#' @export
#'
find_best_xyp <- function(cfg, lrr, hds, p=NULL) {

  na.out <- tibble::tibble(x=1L, y=1L, p=NA_real_, dist=NA_real_,
                   lrr.pred=0.0, hds.pred=0.0)

  if (is.na(hds)) {
    return(na.out)
  }

  ###!!!
  ### What is the rationale behind to select a scale factor WW=5?
  ### The range of hds is (0, 0.5) while the range of R is (0, max.ploidy/2).
  ### What distance of (hds - hds.p) and (R - R.p) (where R = 2^lrr) are considered to be
  ### equivalent?
  ###!!!

  WW = 5 # BubbleTree: weight to reduce the effect [of] R score

  ###!!!
  ### For each curve of copy number state (x,y), we can find a prevalence that
  ### minimizes the distance in HDS-R plane.
  ### We should not need to use grid search because there are indistinguishable states.
  ###!!!

  # out <- plyr::ddply(cfg$xypGrid, .(x,y,p), function(df){
  #   x <- df[1,"x"]
  #   y <- df[1,"y"]
  #   p <- df[1,"p"]
  #   ### tau = copy number
  #   tau <- (x+y)*p + 2*(1-p)
  #   lrr.p <- log2(tau) - 1
  #   hds.p <- ifelse(y==x, 0, p*(y-x)/2/tau)
  #   dist <- sqrt( (2^lrr - tau/2)^2/WW + (hds - hds.p)^2 )
  #   return(c(dist=dist, lrr.pred = lrr.p, hds.pred = hds.p))
  # })

  ### function to find the optimal prevalence value by minimizing distance to the branch
  ### search the best solution from min_p = 0.2 to 1.
  optimal_p <- function(x, y, WW, lrr, hds, min_p) {
    optim(
      par=min_p,
      fn=function(p, x, y, WW, lrr, hds) {
        (2^lrr - (x+y)*p/2 - (1-p))^2 / WW + (hds - 0.5*p*(y-x)/((x+y)*p + 2*(1-p)))^2
      },
      x=x, y=y, WW=WW, lrr=lrr, hds=hds,
      method="Brent",
      lower=min_p, upper=1
    )$par
  }

  if (!is.null(p)) {
    out <- purrr::map2_dfr(cfg$xyGrid$x, cfg$xyGrid$y, function(x,y) {
      if (x == 1 & y == 1) {
        p <- 1
        tau <- 2
        lrr.p <- 0
        hds.p <- 0
        dist <- sqrt( (2^lrr - 1)^2 / WW + hds^2 )
      } else {
        tau <- (x+y)*p + 2*(1-p)
        lrr.p <- log2(tau) - 1
        if (x == y) {
          hds.p <- 0
        } else {
          hds.p <- p*(y-x)/2/tau
        }
        dist <- sqrt( (2^lrr - tau/2)^2 / WW + (hds - hds.p)^2 )
      }
      tibble::tibble(
        x=x, y=y, p=round(p, digits=6),
        dist=round(dist, digits=6),
        lrr.pred = round(lrr.p, digits=4),
        hds.pred = round(hds.p, digits=4)
      )
    })
  } else {
    out <- purrr::map2_dfr(cfg$xyGrid$x, cfg$xyGrid$y, function(x,y) {
      if (x == 1 & y == 1) {
        p <- 1
        tau <- 2
        lrr.p <- 0
        hds.p <- 0
        dist <- sqrt( (2^lrr - 1)^2 / WW + hds^2 )
      } else {
        p <- optimal_p(x, y, WW, lrr, hds, cfg$min.prev)
        ### tau = copy number: R = tau/2
        tau <- (x+y)*p + 2*(1-p)
        lrr.p <- log2(tau) - 1
        if (x == y) {
          hds.p <- 0
        } else {
          hds.p <- p*(y-x)/2/tau
        }
        dist <- sqrt( (2^lrr - tau/2)^2 / WW + (hds - hds.p)^2 )
      }
      tibble::tibble(
        x=x, y=y, p=round(p, digits=3),
        dist=round(dist, digits=6),
        lrr.pred = round(lrr.p, digits=4),
        hds.pred = round(hds.p, digits=4)
      )
    })
  }
  out <- out |> dplyr::filter(dist == min(dist))

  ###!!!
  ### xypGrid contains degenerate (i.e. indistinguishable) states in terms of expected
  ### R, HDS and prevalence. That is, different x,y,p combinations can produce identical
  ### expected values of lrr.p and hds.p.
  ### (x, y) = (1,2), (1,3), (1,4), (1,5) are on the same branch in the R-HDS plot
  ### (x, y) == (2, 2), (3, 3) are on the branch in the R-HDS plot
  ### Therefore, we assume
  ### 1) hds.p of (1,k+1) is higher than the max hds.p (reached at p = 1) of (1, k) for k in (2, ..., 5).
  ###  For k >= 3, keep only the states (1,k) with hds.p > 0.5*((k - 1) - 1)/(1 + (k - 1)) = 0.5*(k-2)/k
  ### 2) lrr.pred of (3,3) is higher than the max lrr.p of (2, 2)
  ###  For k >= 3, keep only the states (k,k) with lrr.p > log2(2k-2) - 1
  ###  equivalently, 2^lrr.p > (2k-2)/2
  ### If max.ploidy = 8
  ### AABBB and AAABBBBB
  #### max.ploidy = 10
  ### AABBBB and AAABBBBBBB
  ### Generally, we can choose the state with minimum y.
  ###!!!

  ###!!!
  ### It is unclear why we need to use cfg$max.distL and cfg$max.distR to
  ### subset candidate copy number states (x, y, p).
  ###!!!
  #only keep the best one for each x y combo
  # out <- subset(out,
  #               (lrr < 0.1 & dist < cfg$max.distL) | (lrr>=0.1 & dist < cfg$max.distR))


  ### 1) Remove cnv states that are closest to the selected (by cfg$min.prev) root of the branch.
  ### 2) We often see the noise in the copy ratios without showing LOH. Remove homozygous losses or
  ### balanced amplification if prevalence is too low (< cfg$min.prev.hom = 0.3-0.4).

  out <- out |> dplyr::filter(p > cfg$min.prev)
  if (nrow(out) > 0) {
    out <- out |> dplyr::filter(!(x != 1 & x == y & p < cfg$min.prev.hom))
  }

  if (nrow(out) == 0) {
    return(na.out)
  }

  return(out)
}
