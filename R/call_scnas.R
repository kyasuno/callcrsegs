#' @title Call SCNAs
#'
#' @description This function infers ploidy of the sample, tumor purity and subclones
#' @param rbd, segment data obtained by [make_input()].
#' @param config, a list of configuration parameters obtained by [make_input()].
#' @param sex female or male.
#' @returns A list.
#' * ploidy.adjust. a numeric vector containing upperCentral.cov, lowHDS.cov,
#' adj (a factor to shift lrr), ploidy and purity (the maximum prevalence).
#' * est.ploidy. a numeric value. ploidy estimated from adj and purity.
#' * prevalence. a tbl with cloneID and prevalence.
#' * result. a tbl.
#'     - lrr50 = original lrr while lrr is adjusted one (lrr50 + log2(adj))
#'     - x0, y0, p0, dist0, lrr.pred0, hds.pred0. The best configuration obtained by
#'     unconstrained values of p (config$min.prev <= p <= 1).
#'     - x, y, p, dist, lrr.pred, hds.pred. The final configuration derived by using
#'     observed prevalence of clones. They are calculated only when !(x0 == 1 & y0 == 1).
#' * deviation. sum(seg.size &times; dist)/100.
#' @export
#'
call_scnas <- function(rbd, config, sex) {
  # make_input has already removed chrX for male samples
  if (sex == "male") {
    rbd <- rbd |> dplyr::filter(seqnames != "chrX")
  }
  pladj <- get_ploidy_adjustment(rbd=rbd, cfg=config)

  rbd.adj <- rbd |> dplyr::mutate(lrr = lrr + log2(pladj["adj"]))
  clones <- find_clones(rbd.adj=rbd.adj, cfg=config, pladj=pladj)

  ## In BubbleTree, when ploidy is predicted to be 4 and
  ## there is no low-prevalence tumor subpopulation (p < 0.5),
  ## BubbleTree considers that ploidy of 2 should be the answer.
  ## It is unclear if it is reasonable.
  ## We use different criteria.

  ## In clones$result, check if there is any SCNAs of type B or BB with prevalence >= 0.3 and
  ## seg.size >= cfg$min.segSize. If there is none, we accept the tetraploidy.

  if (pladj["ploidy"] == 4) {
    p <- clones$prevalence$prevalence
    if (is.na(p[1])) { # there is no candidate SCNA
      pladj["ploidy"] <- 2
    } else if (any(p >= 0.5)) { # there is any p >= 0.5
      pladj["ploidy"] <- 2
    }
    # else {
    #   BorBB <- clones$result |>
    #     dplyr::filter((x == 0 & y == 1) | (x == 0 & y == 2)) |>
    #     dplyr::filter(p >= 0.3 & seg.size > config$min.segSize)
    #   if (nrow(BorBB) > 0) {
    #     pladj["ploidy"] <- 2
    #   }
    # }
    # if the ploidy is still 4,
    if (pladj["ploidy"] == 4) {
      # shift lrr and find clones
      new.adj <- 1/(1-max(p))
      pladj4 <- pladj

      pladj4["adj"] <- pladj4["adj"] * new.adj
      message("call_scna: re-running find_clones for ploidy = 4.")
      rbd.adj4 <- rbd |> dplyr::mutate(lrr=lrr+log2(pladj4["adj"]))
      clones4 <- find_clones(rbd.adj=rbd.adj4, cfg=config, pladj=pladj4)

      if (config$report.tetraploidy) {
        message("call_scna: ploidy = 4 is reported.")
        pladj <- pladj4
        rbd.adj <- rbd.adj4
        clones <- clones4
      } else {
        p <- clones4$prevalence$prevalence
        if (any(!is.na(p))) {
          pladj4["purity"] <- max(p, na.rm=TRUE)
          est.ploidy <- (2*pladj4["adj"] - 2) / pladj4["purity"] + 2
          ## if estimated ploidy is significantly lower than 4, it is unlikely to be tetraploidy
          if (est.ploidy < 3.75) {
            message("call_scna: ploidy = 2 is selected for low estimated ploidy than expected 4: ", round(est.ploidy, digits=4))
            pladj["ploidy"] <- 2
          } else {
            ## check if lrr change is consistent with gain and loss (from AABB)
            loss <- clones4$result |> dplyr::filter(x + y < 4)
            if (nrow(loss) > 0) {
              OK <- all(loss$lrr - pladj4["adj"] < 0)
            } else {
              OK <- TRUE
            }
            gain <- clones4$result |> dplyr::filter(x + y > 4)
            if (nrow(gain) > 0) {
              OK <- OK & all(gain$lrr - pladj4["adj"] > 0)
            }
            if (OK) {
              message("call_scna: ploidy = 4 is selected.")
              pladj <- pladj4
              rbd.adj <- rbd.adj4
              clones <- clones4
            } else {
              message("call_scna: ploidy = 2 is selected for inconsistent SCNA types.")
              pladj["ploidy"] <- 2
            }
          }
        } else {
          message("call_scna: ploidy = 2 is selected.")
          pladj["ploidy"] <- 2
        }
      }
    }
  }

  original_lrr <- rbd |>
    dplyr::select(seg.id, lrr) |> dplyr::rename(lrr50=lrr)
  if (any(colnames(clones$result) == "x0")) {
    result <- clones$result |>
      dplyr::left_join(original_lrr, by="seg.id") |>
      dplyr::select(
        seg.id:width,
        seg.size, num.mark, num.hets, hds.cnt,
        lrr10, lrr50, lrr, lrr90, hds10:hds90,
        x:prevalence, usedForCloneEst, x0:hds.pred0, hds.median:hds.sd
      )
  } else {
    result <- clones$result |>
      dplyr::left_join(original_lrr, by="seg.id") |>
      dplyr::select(
        seg.id:width,
        seg.size, num.mark, num.hets, hds.cnt,
        lrr10, lrr50, lrr, lrr90, hds10:hds90,
        x:prevalence, usedForCloneEst, everything()
      )
  }

  total.size <- sum(rbd$seg.size)
  deviation <- with(
    dplyr::filter(result, seg.size > config$min.segSize, !is.na(dist)),
    sum(seg.size*dist) / 100
  )
  if (is.na(pladj["purity"])) {
    p <- clones$prevalence$prevalence
    if (any(!is.na(p))) {
      pladj["purity"] <- max(p, na.rm=TRUE)
    }
  }
  ## estimated ploidy
  est.ploidy <- (2*pladj["adj"] - 2) / pladj["purity"] + 2


  list(
    ploidy.adj=pladj,
    est.ploidy=est.ploidy,
    prevalence=clones$prevalence,
    result=result,
    deviation=deviation
  )
}


