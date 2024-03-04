#' @title Plot LOH data
#'
#' @param res result of [call_scnas].
#' @param centromere.hg38 Use `data("centromere.hg38")`.
#' @param hg38.seqinfo Use `data("centromere.hg38")`.
#' @param hets A tbl containing heterozygous site data (seqnames, start, end, freq).
#' @param min.prev minimum prevalence (default: 0.2).
#' @param baseline.adj Upward (y = n(B)) and downward (x = n(A)) shift of A and B allele (default: 0.1).
#' @param space.skip Space between chromosomes (default: 0.01)
#' @returns A ggplot object
#' @import ggplot2
#' @export
#'

plot_HetLOH <- function(res, centromere.hg38, hg38.seqinfo, hets=NULL,
                        min.prev=0.2, baseline.adj=0.1, space.skip=0.005) {

  cols <- RColorBrewer::brewer.pal(8,"Dark2")
  col.p <- RColorBrewer::brewer.pal(8,"Paired")

  # check if the sample is female by looking at chrX segments
  chroms <- res$result$seqnames |> as.character() |> unique()
  if (any(chroms == "chrX")) {
    is.female <- TRUE
  } else {
    is.female <- FALSE
  }

  if (!is.female) {
    hg38.seqinfo <- hg38.seqinfo[paste0("chr", c(1:22))]
  } else {
    hg38.seqinfo <- hg38.seqinfo[paste0("chr", c(1:22,"X"))]
  }

  ploidy.adj <- res$ploidy.adj
  prevalence <- res$prevalence
  deviation <- res$deviation
  d <- res$result |>
    dplyr::mutate(
      seqnames=factor(as.character(seqnames), levels=GenomeInfoDb::seqlevels(hg38.seqinfo))
    )

  n.clones <- nrow(prevalence)
  main.prev <- prevalence$prevalence[1]

  # convert 1,1,1 to 1,1, prev
  d <- d |>
    dplyr::mutate(
      y.adj=y + baseline.adj,
      x.adj=x - baseline.adj
    ) |>
    dplyr::mutate(
      # p = if_else(x == 1 & y == 1 & p == 1, main.prev, p)
      prevalence = dplyr::if_else(x == 1 & y == 1, main.prev, prevalence)
    )

  ## convert chromosomal coordinates to genomic coordinates
  d.genome <- get_genomic_coord(d, hg38.seqinfo, space.skip=space.skip)
  d <- d.genome$data

  # centromere data
  centromere.hg38 <- centromere.hg38 |>
    dplyr::filter(chrom %in% GenomeInfoDb::seqlevels(hg38.seqinfo)) |>
    dplyr::rename(seqnames=chrom, start=chromStart, end=chromEnd) |>
    dplyr::mutate(start=start + 1)
  c.genome <- get_genomic_coord(centromere.hg38, hg38.seqinfo, space.skip=space.skip)
  centromere.hg38 <- c.genome$data

  p <- ggplot(inherit.aes=FALSE)
  if (!is.null(hets)) {
    hets <- hets |> dplyr::filter(seqnames %in% GenomeInfoDb::seqlevels(hg38.seqinfo))
    hets.genome <- get_genomic_coord(hets, hg38.seqinfo, space.skip=space.skip)
    hets <- hets.genome$data
    p <- p +
      geom_point(data=hets, aes(x=.start, y=freq),
                 colour=col.p[3], size=1/ggplot2::.pt, alpha=0.2,
                 inherit.aes=FALSE)
  }

  ## observed BAF
  p <- p +
    geom_rect(
      data=d,
      aes(xmin=.start, xmax=.end, ymin=0.5 - hds - 0.01, ymax=0.5 - hds + 0.01),
      colour=cols[2], fill=cols[2], na.rm=TRUE, inherit.aes=FALSE) +
    geom_rect(
      data=d,
      aes(xmin=.start, xmax=.end, ymin = 0.5 + hds - 0.01, ymax=0.5 + hds + 0.01),
      colour=cols[2], fill=cols[2], na.rm=TRUE, inherit.aes=FALSE)

  ## predicted BAF
  p <- p +
    geom_rect(
      data=d,
      aes(xmin=.start, xmax=.end, ymin=0.5 - hds.pred - 0.01, ymax=0.5 - hds.pred + 0.01),
      colour="#000000", fill="#000000", na.rm=TRUE, inherit.aes=FALSE) +
    geom_rect(
      data=d,
      aes(xmin=.start, xmax=.end, ymin = 0.5 + hds.pred -0.01, ymax=0.5 + hds.pred + 0.01),
      colour="#000000", fill="#000000", na.rm=TRUE, inherit.aes=FALSE) +
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,by=0.25), expand=c(0,0)) +
    labs(y="BAF", x="")

  # centromere
  p <- p +
    scale_x_continuous(breaks = (centromere.hg38$.start + centromere.hg38$.end) / 2,
                       labels=as.character(centromere.hg38$seqnames),
                       expand=c(0,0)) +
    theme(
      panel.background = element_blank(),
      panel.grid.major.y = element_line(colour="gray90", linewidth=0.3),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10)
    )

  if(n.clones > 1) { # highlight subclones
    sp <- prevalence$prevalence[-1] |> sort(decreasing = TRUE)
    dd <- d |>
      dplyr::filter(prevalence <= sp[1]) |>
      dplyr::mutate(p_factor=factor(prevalence, levels=sp))
    cc <- RColorBrewer::brewer.pal(8,"Set2")[1:length(sp)]

    p <- p + geom_rect(data=dd,
                       aes(xmin=.start, xmax=.end, ymin=-Inf, ymax=Inf, fill=p_factor),
                       alpha=0.2, inherit.aes = FALSE, show.legend = FALSE) +
      scale_fill_manual(values=cc)
  }
  return(p)
}
