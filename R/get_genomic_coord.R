#' @title Plot allele-specific copy number changes across genome
#'
#' @param df A tbl.
#' @param hg38.seqinfo Use `data("centromere.hg38")`
#' @param space.skip Space between chromosomes (default: 0.01)
#' @returns A ggplot object
#' @import ggplot2
#' @export

get_genomic_coord <- function(df, hg38.seqinfo, space.skip=0.1) {
  # this is a data.frame version of biovizBase::transformToGenome
  contig_names <- GenomeInfoDb::seqlevels(hg38.seqinfo)
  contig_lengths <- GenomeInfoDb::seqlengths(hg38.seqinfo)
  total_lengths <- sum(contig_lengths)
  space.skip <- space.skip * total_lengths

  skps <- space.skip * ((1:length(contig_names))-1)
  names(skps) <- contig_names

  contig_ends <- cumsum(as.numeric(contig_lengths))
  contig_starts <- c(0, head(contig_ends, -1))
  names(contig_starts) <- contig_names
  df <- df |>
    dplyr::mutate(
      .start=contig_starts[match(seqnames, contig_names)] + skps[match(seqnames, contig_names)] + start,
      .end=contig_starts[match(seqnames, contig_names)] + skps[match(seqnames, contig_names)] + end
    )

  max.chr <- rev(contig_names)[1]
  x.max <- contig_lengths[max.chr] + skps[max.chr] + contig_starts[max.chr] + space.skip

  breaks <- contig_starts + contig_lengths/2  + skps

  return(
    list(data=df, x.max=x.max, breaks=breaks, labels=contig_names, space.skip=space.skip)
  )
}
