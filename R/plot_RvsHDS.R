#' @title Plot copy ratio vs HDR
#'
#' @param res result of [call_scnas].
#' @param tumorID Tumor ID.
#' @param min.prev minimum prevalence (default: 0.2).
#' @returns A ggplot object
#' @import ggplot2
#' @export
plot_RvsHDS <- function(res, tumorID, min.prev=0.2, max.size=20) {
  cols1 <- wesanderson::wes_palette("Darjeeling1", n=23, type="continuous")
  cols1 <- c(cols1, "#808080")
  names(cols1) <- c(paste0("chr", c(1:22,"X")), "CN")

  d <- res$result
  ploidy.adj <- res$ploidy.adj

  max.cn <- with(d, x + y) |> max()
  zoom.size <- max(5, max.size/max(d$seg.size))

  d <- d |>
    dplyr::filter(!is.na(hds)) |>
    dplyr::mutate(
      bubble.size = pmin(seg.size * zoom.size, max.size)
    ) |>
    dplyr::mutate(
      SCNAchromosome=as.character(seqnames)
    )

  if (ploidy.adj["ploidy"] == 2) {
    d <- d |>
      dplyr::mutate(
        SCNAchromosome=dplyr::if_else(x == 1 & y == 1, "CN", SCNAchromosome)
      )
  } else if (ploidy.adj["ploidy"] == 3) {
    d <- d |>
      dplyr::mutate(
        SCNAchromosome=dplyr::if_else(x == 1 & y == 2, "CN", SCNAchromosome)
      )
  } else if (ploidy.adj["ploidy"] == 4) {
    d <- d |>
      dplyr::mutate(
        SCNAchromosome=dplyr::if_else(x == 2 & y == 2, "CN", SCNAchromosome)
      )
  }

  chromsObs <- d$SCNAchromosome |> unique()
  chroms <- levels(d$seqnames)[levels(d$seqnames) %in% chromsObs]
  cols1 <- cols1[c(chroms, "CN")]

  d <- d |>
    dplyr::mutate(
      SCNAchromosome=factor(SCNAchromosome, levels=c(chroms, "CN"))
    ) |>
    dplyr::arrange(desc(SCNAchromosome))

  # nSCNA is different between assumed ploidy
  if (res$ploidy.adj["ploidy"] == 2) {
    nSCNAs <- res$result |> dplyr::filter(!(x == 1 & y == 1)) |> nrow()
  } else if (res$ploidy.adj["ploidy"] == 3) {
    nSCNAs <- res$result |> dplyr::filter(!(x == 1 & y == 2)) |> nrow()
  } else if (res$ploidy.adj["ploidy"] == 4) {
    nSCNAs <- res$result |> dplyr::filter(!(x == 2 & y == 2)) |> nrow()
  }

  nOut <- sum(2^d$lrr > max(3, max.cn/2))
  if (nOut == 0) {
    captext <- ""
  } else if (nOut == 1) {
    captext <- paste0(nOut, " event with R > ", max(3, max.cn/2), " is not shown")
  } else {
    captext <- paste0(nOut, " events with R > ", max(3, max.cn/2), " are not shown")
  }

  p <- plot_xypGrid(max.ploidy=max(6, max.cn), min.prev=min.prev)
  p <- p +
    geom_point(data=d |> dplyr::filter(2^lrr <= max(3, max.cn/2)),
               aes(x=2^lrr, y=hds, size=bubble.size, fill=SCNAchromosome),
               shape=21, alpha=0.75,
               inherit.aes=FALSE) +
    scale_fill_manual(values=cols1) +
    guides(
      size="none",
      fill = guide_legend(nrow=12, override.aes=list(shape = 21,alpha=0.5, size=4))
    ) +
    labs(
      fill="SCNA\nChromosome",
      title= paste0(tumorID, "\n",
                    sprintf(
                      "Purity: %s\nPloidy: %3.1f; Deviation: %6.4f; nSCNAs: %3.0f",
                      paste(round(res$prevalence$prevalence, 2), collapse=", "),
                      round(res$est.ploidy, 2),
                      res$deviation,
                      nSCNAs)),
      caption=captext,
      x="R", y="HDS"
    ) +
    theme(title=element_text(size=8))
  p
}
