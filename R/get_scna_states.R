#' @title Estimate tumor purity and subclones
#'
#' @description This function estimates the most likely configuration of copy number states
#' @param rbd.adj, segment data obtained by [get_ploidy_adjustment], in which lrr has been adjusted.
#' @param cfg, config data obtained by [make_input()].
#' @param p a numeric vector. If not NULL, distance to the branches are calculated only for given p's.
#' @returns A list.
#' @export
#'

get_scna_states <- function(rbd.adj, cfg, p=NULL) {
  rbd.dist <- purrr::pmap_dfr(list(rbd.adj$seg.id, rbd.adj$lrr, rbd.adj$hds), function(seg.id, lrr, hds) {
    find_best_xyp(cfg, lrr, hds, p=p) |> dplyr::mutate(seg.id=seg.id)
  })
  rbd.dist <- dplyr::right_join(rbd.adj, rbd.dist, by="seg.id")
  ### Remove indistinguishable states if exist (prefer higher purity state <=> lower y)
  rbd.dist <- rbd.dist |>
    dplyr::group_by(seg.id) |>
    dplyr::filter(y == min(y)) |>
    dplyr::ungroup()
  return(rbd.dist)
}
