#' @title Get a grid of copy number states: (x, y, p)
#'
#' @param max.ploidy maximum ploidy to consider.
#' @param min.prev minimum prevalence.
#' @returns A tbl
#' @export
#'

get_xypGrid <- function(max.ploidy, min.prev) {
  xypGrid <- expand.grid(x=0:floor(max.ploidy/2), y=0:max.ploidy, p=seq(min.prev, 1, by=0.01)) |>
    tibble::as_tibble() |>
    dplyr::filter (x <= y & x+y <= max.ploidy) |>
    dplyr::filter( ! (x==1 & y==1) ) |>
    unique()
  xypGrid <- xypGrid |>
    dplyr::mutate(
      p = round(p, digits=3),
      tau = (x+y)*p + 2*(1-p),
      lrr.pred = log2(tau) - 1,
      hds.pred =  dplyr::if_else(y==x, 0, p*(y-x)/2/tau),
      GT=purrr::map2_chr(x, y, function(x,y) { paste(c(rep("A", x), rep("B", y)), collapse="") }),
      ploidy=x+y
    ) |>
    dplyr::mutate(
      GT=dplyr::if_else(GT == "", "0", GT),
      lrr.pred=round(lrr.pred, digits=4),
      hds.pred=round(hds.pred, digits=4)
    ) |>
    dplyr::select(GT, ploidy, x:hds.pred)
  return(xypGrid)
}

