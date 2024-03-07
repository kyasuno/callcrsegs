#' @title Infer sample ploidy and estimate adjustment factor for copy ratio
#'
#' @description This function infers sample ploidy and estimate adjustment factor for copy ratio.
#' Original function is getPloidyAdj in BubbleTree.
#' @param rbd, segment data obtained by [make_input()].
#' @param cfg, config data obtained by [make_input()].
#' @returns A named vector.
#' @export
#'
get_ploidy_adjustment <- function(rbd, cfg) {

  # parameters
  centralSegInterval <- cfg$centralSegInterval # central segments: abs(lrr) < 0.1 or (abs(lrr - all.seg.mean) < 0.1)
  uc.hds.cutoff <- cfg$uc.hds.cutoff     # upper central segments: subset of central segments with hds > 0.07
  uc.cov.cutoff <- cfg$uc.cov.cutoff      # minimum coverage of upper central segments required for triploidy (0.5)
  lowHds.cutoff <- cfg$lowHds.cutoff     # define low HDS segments by hds < 0.15

  ### total segment size in units of number of targets
  total.size <- sum(rbd$seg.size, na.rm=TRUE)

  ### weighted median of lrr for all segments
  all.seg.mean <- with(rbd, limma::weighted.median(lrr, seg.size, na.rm=TRUE))

  ### central segments
  ###!!!
  ### It is not reasonable to use lrr - all.seg.mean when there are large number of gains
  ###!!!
  # if (abs(all.seg.mean) > 0.3) {
  #   # special abnormal case like 1500
  #   central.segs <- subset(rbd, abs(lrr) < centralSegInterval)
  #   all.seg.mean <- with(central.segs,
  #                        limma::weighted.median(lrr, seg.size, na.rm=TRUE))
  # } else {
  #   central.segs <- subset(rbd, abs(lrr - all.seg.mean) < centralSegInterval)
  # }

  if (cfg$use.seg.median.for.center) { ### better to turn it on especially when the data is noisy
    if (abs(all.seg.mean) > 0.3) { ### we need this condition
      central.segs <- rbd |> dplyr::filter(abs(lrr) < centralSegInterval)
      all.seg.mean <- with(central.segs, limma::weighted.median(lrr, seg.size, na.rm=TRUE))
    } else {
      central.segs <- rbd |> dplyr::filter(abs(lrr - all.seg.mean) < centralSegInterval)
    }
  } else {
    central.segs <- rbd |> dplyr::filter(abs(lrr) < centralSegInterval)
  }

  if (nrow(central.segs) == 0) {
    #!!!
    #!!! do something: should we stop here?
    #!!!
    warning("get_ploidy_adjustment: ",
            "There is no segment that can be used to define copy-neutral normal state.\n",
            "Shift factor of lrr is calculated by using all the segments and ploidy is assumed 2.")
    ploidy.adj <- c(
      upperCentral.cov = NA_real_,
      lowHDS.cov = NA_real_,
      adj=2^all.seg.mean,
      ploidy=2,
      purity=NA_real_
    )
    return(ploidy.adj)
  }

  ### upper central segments and its coverage fraction
  uc.segs <- central.segs |> dplyr::filter(!is.na(hds) & hds > uc.hds.cutoff)
  uc.cov <- sum(uc.segs$seg.size) / sum(central.segs$seg.size)

  ### weighted median of HDS across upper central segments
  uc.hds <- 0
  if (nrow(uc.segs) > 0) {
    uc.hds <- limma::weighted.median(uc.segs$hds, uc.segs$seg.size, na.rm=TRUE)
  }

  ### coverage of low HDS (< 0.15) segments
  low.hds <- with(rbd, sum(seg.size[!is.na(hds) & hds < lowHds.cutoff])) / total.size
  has.high.hds <- with(rbd, sum(seg.size > cfg$min.segSize & !is.na(hds) & hds >= lowHds.cutoff)) > 0

  ### determine purity only when ploidy = 3. Otherwise, set NA
  ### replaced hard coded value of 0.5 by uc.cov.cutoff
  ###
  if (uc.cov > uc.cov.cutoff & uc.hds < 1/6 & cfg$allow.triploidy) {

    ### adj = ( ploidy * p + 2*(1 - p) ) / 2, where p = purity = max(prev)
    ### R = (1-p) + 3*p/2 = p/2 + 1
    ### R = R' * adj
    ### hds = p / (2*(3*p + 2*(1-p))) = p / (2*(p + 2))
    ### 2*hds = p / (p+2)
    ### 2*hds*(p+2) = p
    ### (1 - 2*hds)*p = 4*hds
    ploidy <- 3
    uc.lrr <- with(uc.segs, limma::weighted.median(lrr, seg.size, na.rm=TRUE)) # na.rm=TRUE is added
    purity <- 4*uc.hds/(1-2*uc.hds)
    adj <- (1+purity/2) /2^uc.lrr

    ###!!!
    ### do we need this (it is almost likely that there is no central.segs.new)??
    ###!!!
    # also check the new central stages to make the fine adjustment
    # central.segs.new <- subset(rbd, abs(lrr + log2(adj)-all.seg.mean) < centralSegInterval)
    # change the above not to use all.seg.mean
    if (cfg$use.seg.median.for.center) {
      central.segs.new <- rbd |> dplyr::filter(abs(lrr + log2(adj) - all.seg.mean) < centralSegInterval)
    } else {
      central.segs.new <- rbd |> dplyr::filter(abs(lrr + log2(adj)) < centralSegInterval)
    }

    if (nrow(central.segs.new) > 0 ) { # unlikely to exist
      adj <-  with(central.segs.new, 2^(-limma::weighted.median(lrr, seg.size, na.rm=TRUE)))
    }

  } else {

    ploidy <- 2

    ### when ploidy = 2, we cannot estimate purity here
    purity <- NA

    ###!!! Bug in BubbleTree: should use weighted.median !!!###
    # lrr <- with(central.segs, weighted.mean(lrr, seg.size))
    lrr <- with(central.segs, limma::weighted.median(lrr, seg.size, na.rm=TRUE))
    adj <- 1 / 2^lrr

    ###!!!
    ### BubbleTree: ploidy = 4 when low.hds > 95% and
    ### the coverage of upper right segments (hds > 0.1 & 2^lrr > 1.25) < 1%
    ### => we allow tetraploidy only when
    ### the user wants to include tetraploidy case,
    ### low.hds > 0.9995,
    ### no segment with seg.size > cfg$min.segSize and hds >= lowHds.cutoff, and
    ### no upper right segment (like ABB) with seg.size > cfg$min.segSize
    ###!!!
    if ( low.hds > 0.9995 & !has.high.hds & cfg$allow.tetraploidy ) {
      ### ABB, p = 0.35 (lrr.p = 1.175, hds.p = 0.07446809)
      upright.segs <- rbd |> dplyr::filter(!is.na(hds), hds > 0.075 & 2^lrr > 1.175, seg.size >= cfg$min.segSize)
      if (nrow(upright.segs) == 0) {
        ploidy <- 4
      }
    }
  }

  ### output
  ploidy.adj <- c(
    upperCentral.cov = uc.cov,
    lowHDS.cov = low.hds,
    adj=adj,
    ploidy=ploidy, #as.integer(ploidy),
    purity=purity
  )

  return(ploidy.adj)
}
