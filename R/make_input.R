#' @title Make input data
#'
#' @description This function generates dataset to be used in this package.
#' @param input_dir Input directory where <tumorID>.modelFinal.seg, <tumorID>.hets.tsv and
#' <tumorID>.denoisedCR.tsv can be found.
#' @param tumorID Sample ID of tumor.
#' @param sex female, male or a filename that contains a word, female or male.
#' If missing, sex is assumed to be male and only autosomes are analyzed.
#' @returns A list containing:
#' * tumorID.
#' * sex.
#' * rbd, input for the SCNA calling.
#' * cnv.nohets, segments without heterozygous sites.
#' * crs, copy ratio data used for plotting.
#' * segs, segmentation result from GATK where seg.id has been added and
#' posterior MAFs have been converted to HDS.
#' * mafs, MINOR_ALLELE_FRACTION_POSTERIOR_50 for plotting.
#' * hets, ALT allele frequency at heterozygous sites in normal sample.
#' @export
#'

make_input <- function(input_dir, tumorID, sex) {
  if (missing(sex)) {
    sex <- "male"
  }
  if (!sex %in% c("female", "male")) {
    sex <- read_lines(sex)
  }
  if (!dir.exists(input_dir)) {
    stop("input_dir: ", input_dir, " was not found.")
  }
  ## check files
  crs.file <- file.path(input_dir, paste0(tumorID, ".denoisedCR.tsv"))
  if (!file.exists(crs.file)) {
    stop(crs.file, " was not found.")
  }
  mfsegs.file <- file.path(input_dir, paste0(tumorID, ".modelFinal.seg"))
  if (!file.exists(mfsegs.file)) {
    stop(mfsegs.file, " was not found.")
  }
  hets.file <- file.path(input_dir, paste0(tumorID, ".hets.tsv"))
  if (!file.exists(hets.file)) {
    stop(hets.file, " was not found.")
  }
  crs <- read_tsv(crs.file, comment="@", progress=FALSE, show_col_types = FALSE)
  names(crs) <- c("seqnames", "start", "end", "lrr")
  mfsegs <- read_tsv(mfsegs.file, comment="@", progress=FALSE, show_col_types = FALSE) |>
    mutate(seg.id=1:n())
  hets <- read_tsv(hets.file, comment="@", progress=FALSE, show_col_types = FALSE) |>
    mutate(freq=ALT_COUNT / (REF_COUNT + ALT_COUNT)) |>
    filter(!is.na(freq))

  cnv.gr <- GRanges(
    seqnames=mfsegs$CONTIG,
    IRanges(start=mfsegs$START, end=mfsegs$END),
    seg.id=mfsegs$seg.id,
    num.mark=mfsegs$NUM_POINTS_COPY_RATIO,
    lrr10=mfsegs$LOG2_COPY_RATIO_POSTERIOR_10,
    lrr=mfsegs$LOG2_COPY_RATIO_POSTERIOR_50,
    lrr90=mfsegs$LOG2_COPY_RATIO_POSTERIOR_90,
    num.hets=mfsegs$NUM_POINTS_ALLELE_FRACTION,
    hds10=if_else(!is.na(mfsegs$MINOR_ALLELE_FRACTION_POSTERIOR_90),
                  abs(mfsegs$MINOR_ALLELE_FRACTION_POSTERIOR_90 - 0.5), NA_real_),
    hds=if_else(!is.na(mfsegs$MINOR_ALLELE_FRACTION_POSTERIOR_50),
                abs(mfsegs$MINOR_ALLELE_FRACTION_POSTERIOR_50 - 0.5), NA_real_),
    hds90=if_else(!is.na(mfsegs$MINOR_ALLELE_FRACTION_POSTERIOR_10),
                  abs(mfsegs$MINOR_ALLELE_FRACTION_POSTERIOR_10 - 0.5), NA_real_)
  )
  snp.gr <- GRanges(
    seqnames=hets$CONTIG,
    IRanges(start=hets$POSITION, end=hets$POSITION),
    freq=hets$freq
  )

  # combine cnv.gr and snp.gr
  hits <- GenomicRanges::findOverlaps(snp.gr, cnv.gr) |> as_tibble()
  cnv.df <- cnv.gr |> as_tibble() |> dplyr::rename(cnv.start=start, cnv.end=end)
  snp.df <- snp.gr |> as_tibble()
  rbd <- cbind(
    snp.df[hits$queryHits, ],
    cnv.df[hits$subjectHits, ] |>
      dplyr::select(seg.id, num.mark:lrr90, cnv.start, cnv.end, num.hets:hds90)
  ) |>
    as_tibble()

  # keep records of segments without heterozygous sites
  cnv.missing <- cnv.df |> filter(!(seg.id %in% rbd$seg.id))

  # calculate HDS
  rbd <- rbd |>
    group_by(seg.id, seqnames, cnv.start, cnv.end,
             num.mark, lrr10, lrr, lrr90,
             num.hets, hds10, hds, hds90) |>
    filter(!is.na(num.mark) & !is.na(seg.id) & !is.na(freq)) |>
    summarise(
      hds.median=median(abs(freq - 0.5)),
      hds.mad=if_else(length(freq) > 1, mad(abs(freq - 0.5)), 0),
      hds.mean=mean(abs(freq - 0.5)),
      hds.sd=sd(abs(freq - 0.5)),
      hds.cnt=length(freq),
      .groups="drop"
    ) |>
    dplyr::rename(
      start=cnv.start, end=cnv.end
    )

  hets <- hets |>
    dplyr::rename(seqnames=CONTIG) |>
    mutate(start=POSITION, end=POSITION) |>
    dplyr::select(seqnames, start, end, freq)

  mafs <- tibble(
    seqnames=mfsegs$CONTIG,
    start=mfsegs$START,
    end=mfsegs$END,
    seg.id=mfsegs$seg.id,
    n.hets=mfsegs$NUM_POINTS_ALLELE_FRACTION,
    maf=mfsegs$MINOR_ALLELE_FRACTION_POSTERIOR_50
  ) |> filter(!is.na(maf))

  list(
    tumorID=tumorID,
    sex=sex,
    rbd=rbd,
    cnv.nohets=cnv.missing,
    crs=crs,
    segs=cnv.gr |> as_tibble(),
    mafs=mafs,
    hets=hets
  )
}

