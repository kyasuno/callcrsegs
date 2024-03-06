#' @title Estimate tumor purity and subclones
#'
#' @description This function estimates tumor purity and subclones as well as the best copy
#' number states. Original function is findClones in BubbleTree.
#' @param rbd.adj, segment data in which lrr has been adjusted.
#' @param cfg, config data obtained by [make_input()].
#' @param pladj a numeric vector containing the result of [get_ploidy_adjustment].
#' @returns A list.
#' @export
#'

find_clones <- function(rbd.adj, cfg, pladj) {

  ### find the best SCNA state for each segment
  rbd.dist <- get_scna_states(rbd.adj, cfg)

  ### remove segments that are likely to be in a copy-neutral normal state
  ### (for tri/tetra-ploidy cases, they are SCNAs but not informative for prevalence)
  rbd.scna <- rbd.dist |> dplyr::filter(!(x == 1 & y == 1))

  if(nrow(rbd.scna) == 0) { # no prediction
    message("find_clones: There are no segments that are likely to be SCNA.")
    prev <- tibble::tibble(cloneID=1L, prevalence=NA_real_)
    rbd.dist <- rbd.dist |>
      dplyr::mutate(
        cloneID=NA_integer_,
        prevalence=NA_real_
      )
    return(list(prevalence=prev, result=rbd.dist))
  }

  ### select "large" segments with "large" abs(lrr) or hds
  rbd.scna.large <- rbd.scna |>
    dplyr::filter(
      seg.size >= cfg$min.segSize & (hds >= cfg$root.hdsInt | abs(lrr) >= cfg$root.rInt)
    )

  if(nrow(rbd.scna.large) == 0) {
    message("find_clones: There is no segment with seg.size >= ", cfg$min.segSize, ", ",
            "hds >= ", cfg$root.hdsInt, " and abs(lrr) >= ", round(cfg$root.rInt, digits=3), ".")
    prev <- tibble::tibble(cloneID=1L, prevalence=NA_real_)
    rbd.dist <- rbd.dist |>
      dplyr::mutate(
        cloneID=NA_integer_,
        prevalence=NA_real_
      )
    return(list(prevalence=prev, result=rbd.dist))
  }

  ### Exclude segments with large distance from any branches or copy ratio is too high
  rbd.scna.large <- rbd.scna.large |>
    dplyr::filter(dist < 0.11 & lrr < log2(cfg$max.ploidy))

  if(nrow(rbd.scna.large) == 0) { # no prediction
    message("find_clones: There are no segments for which a copy number states are reliably predicted")
    prev <- tibble::tibble(cloneID=1L, prevalence=NA_real_)
    rbd.dist <- rbd.dist |>
      dplyr::mutate(
        cloneID=NA_integer_,
        prevalence=NA_real_
      )
    return(list(prevalence=prev, result=rbd.dist))
  }

  ### Although it is unclear whether the following filtering based on the copy number,
  ### we include it. We can disable this by setting max.ploidy.clone = max.ploidy (default)
  if (cfg$max.ploidy.clone < cfg$max.ploidy) {
    high.cn.segs <- rbd.scna.large |> dplyr::filter( x + y > cfg$max.ploidy.clone)
    if (nrow(high.cn.segs) > 0) {
      message("find_clones: Remove segments with high copy number > ", cfg$max.ploidy.clone, ".")
      low.cn.segs <- rbd.scna.large |> dplyr::filter(x + y <= cfg$max.ploidy.clone)
      if(nrow(low.cn.segs) == 0) {
        message("find_clones: There is no informative segment with copy number <= ", cfg$max.ploidy.clone, ".\n",
                "  We use segments inferred to have higher copy numbers.")
      } else {
        rbd.scna.large <- low.cn.segs
      }
    }
  }

  ## split uniquely identifiable and unidentifiable states
  # ABB, ABBB, ABBBB, ...
  # AABBB and AAABBBBB
  # AABBBB and AAABBBBBBB
  # AABB, AAABBB, ...
  uniqstates <- rbd.scna.large |>
    dplyr::filter(
      !(x == 1 & y >=2),
      !(x == 2 & y == 3), !(x == 3 & y == 5),
      !(x == 2 & y == 4), !(x == 3 & y == 7),
      !(x >= 2 & x == y)
    )
  ambiguous <- rbd.scna.large |>
    dplyr::filter(
      (x == 1 & y >=2) |
        (x == 2 & y == 3) | (x == 3 & y == 5) |
        (x == 2 & y == 4) | (x == 3 & y == 7) |
        (x >= 2 & x == y)
    )

  n_uniq <- nrow(uniqstates)
  n_ambi <- nrow(ambiguous)
  if (n_ambi == 0) {
    message("find_clones: Use identifiable states to determine clones (no unidentifiable states).")
  } else if (n_ambi > 0 & n_uniq == 0) {
    message("find_clones: There is no informative segment without including unidentifiable states. ",
            "Use the unidentifiable states.")
  } else if (n_ambi > 0 & n_uniq > 0) {
    if (max(uniqstates$p) < cfg$lowest.purity) {
      message("find_clones: Highest prevalence (", round(max(uniqstates$p), digits=3),
              ") is lower than required purity: ", cfg$lowest.purity, ". ",
              "Use also unidentifiable states.")
    } else if (!cfg$useABB) {
      message("find_clones: Use only identifiable states to determine clones (unidentifiable states are discarded).")
      rbd.scna.large <- uniqstates
    } else if (cfg$useABB) {
      message("find_clones: Use both identifiable and unidentifiable states to determine clones.")
    }
  } else {
    stop("find_clones: Something went wrong")
  }

  if (nrow(rbd.scna.large) == 1) {
    prev <- tibble::tibble(cloneID=1L, prevalence=rbd.scna.large$p |> round(digits=6))
    rbd.scna.large <- rbd.scna.large |>
      dplyr::mutate(cloneID=1L, prevalence=round(p, digits=6))
  } else {
    # prevalence is calculated as BubbleTree did (weighted.mean). Not limma:weighted.median
    mem <- cutree(hclust(dist(rbd.scna.large$p)), h=cfg$cutree.h)
    rbd.scna.large <- rbd.scna.large |>
      dplyr::mutate(tmpID=mem)
    prev <- rbd.scna.large |>
      dplyr::group_by(tmpID) |>
      dplyr::summarise(
        prevalence=weighted.mean(p, seg.size) |> round(digits=6)
      ) |>
      dplyr::arrange(desc(prevalence)) |>
      dplyr::mutate(cloneID=1:dplyr::n())
    rbd.scna.large <- rbd.scna.large |>
      dplyr::left_join(prev, by="tmpID") |>
      dplyr::select(-tmpID)
    prev <- prev |> dplyr::select(cloneID, prevalence)
  }
  ## prevalence is rounded to be consistent with find_best_xyp
  prev <- prev |>
    dplyr::mutate(prevalence=round(prevalence, digits=6))

  ## find nearest clone for remaining segments
  rbd.scna.other <- rbd.scna |>
    dplyr::filter(!(seg.id %in% rbd.scna.large$seg.id))
  if (nrow(rbd.scna.other) > 0) {
    rbd.scna.other <- rbd.scna.other |>
      dplyr::mutate(
        cloneID=purrr::map_int(p, function(z) {
          dst <- abs(z - prev$prevalence)
          if (any(dst < cfg$cutree.h)) {
            which.min(dst)
          } else {
            -1L
          }
        })
      )
    rbd.scna.other <- rbd.scna.other |>
      dplyr::left_join(prev, by="cloneID")
    rbd.scna <- dplyr::bind_rows(rbd.scna.large, rbd.scna.other)
  } else {
    rbd.scna <- rbd.scna.large
  }
  rbd.dist <- rbd.dist |>
    dplyr::left_join(rbd.scna |> dplyr::select(seg.id, cloneID, prevalence), by="seg.id")

  ## Now recalculate the SCNA states given the prevalence (following BubbleTree)

  if (pladj["ploidy"] == 3) {
    # add presumed prevalence if it is much higher than the estimated purity (max. prevalence)
    pp0 <- pladj["purity"]
    if (prev$prevalence[1] < pp0 - cfg$cutree.h) {
      pp0df <- tibble::tibble(cloneID=0L, prevalence=round(pp0, digits=6))
      prev <- dplyr::bind_rows(pp0df, prev)
    }
  }

  recal <- get_scna_states(rbd.scna |> dplyr::select(-(x:hds.pred), -cloneID, -prevalence),
                           cfg, p=prev$prevalence) |>
    dplyr::rename(x0=x, y0=y, p0=p, dist0=dist, lrr.pred0=lrr.pred, hds.pred0=hds.pred)
  ## merge data
  rbd.dist <- rbd.dist |>
    dplyr::left_join(recal |> dplyr::select(seg.id, x0:hds.pred0), by="seg.id")

  ## fill data for the CN normal segments
  cnsegs <- !(rbd.dist$seg.id %in% rbd.scna$seg.id)
  rbd.dist <- rbd.dist |>
    dplyr::mutate(
      x0=dplyr::if_else(cnsegs, x, x0),
      y0=dplyr::if_else(cnsegs, y, y0),
      p0=dplyr::if_else(cnsegs, p, p0),
      dist0=dplyr::if_else(cnsegs, dist, dist0),
      lrr.pred0=dplyr::if_else(cnsegs, lrr.pred, lrr.pred0),
      hds.pred0=dplyr::if_else(cnsegs, hds.pred, hds.pred0)
    )

  rbd.dist <- rbd.dist |>
     dplyr::select(seg.id:seg.size, x:hds.pred, cloneID, prevalence, everything())


  list(prevalence=prev, result=rbd.dist)
}
