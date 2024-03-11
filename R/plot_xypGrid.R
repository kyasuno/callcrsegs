#' @title Plot a grid of copy number states: (x, y, p)
#'
#' @param max.ploidy maximum ploidy to consider.
#' @param min.prev minimum prevalence.
#' @returns A ggplot object
#' @import ggplot2
#' @export
#'
plot_xypGrid <- function(max.ploidy, min.prev) {
  cols <- RColorBrewer::brewer.pal(7,"Dark2")
  xypGrid <- get_xypGrid(max.ploidy, min.prev) |>
    dplyr::mutate(R.pred=2^lrr.pred, Prev=round(100*p))
  label_data <- xypGrid |>
    dplyr::filter(p == 1) |>
    dplyr::mutate(
      R.pred=round(R.pred, digits=3),
      rx=R.pred,
      hx=hds.pred
    ) |>
    dplyr::mutate(
      hds.pred=dplyr::if_else(rx == max.ploidy / 2 | hx == 0 | hx == 0.5,
                              hds.pred + 0.01, hds.pred - 0.01)#,
      #R.pred=if_else(rx == max.ploidy / 2 | hx == 0 | hx == 0.5, R.pred, R.pred + 0.15)
    )
  lv <- label_data |> dplyr::filter(rx == max.ploidy / 2 | hx == 0 | hx == 0.5)
  lh <- label_data |> dplyr::filter(!(rx == max.ploidy / 2 | hx == 0 | hx == 0.5))
  xypGrid |>
    dplyr::filter(p %in% unique(c(min.prev, 0.2, 0.4, 0.6, 0.8, 1))) |>
    ggplot(aes(x=R.pred, y=hds.pred, group=GT)) +
    theme_minimal(base_size=10) +
    geom_text(aes(label=Prev), size=7/.pt, colour="black") +
    geom_line(linetype="dotted", colour="gray") +
    geom_text(data=lv, aes(label=GT), size=7/.pt) +
    geom_text(data=lh, aes(label=GT), size=7/.pt, hjust=-0.5) +
    scale_x_continuous(breaks=seq(0, ceiling(max.ploidy/2), by=0.25)) +
    theme(
      axis.text=element_text(size=7),
      axis.title=element_text(size=8)
    )
}
