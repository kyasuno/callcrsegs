#' @title Make input data
#'
#' @description This function generates dataset to be used in this package.
#' @param input_dir Input directory where <tumorID>.modelFinal.seg, <tumorID>.hets.tsv and
#' <tumorID>.denoisedCR.tsv can be found.
#' @param tumorID Sample ID of tumor.
#' @param sex female, male or a filename that contains a word, female or male.
#' @param use.physical.length use physical segment lengths to calculate proportional lengths.
#' By default (FALSE), number of targets is used.
#' If missing, sex is assumed to be male and only autosomes are analyzed.
#' @returns A list containing:
#' * tumorID.
#' * sex.
#' * rbd, input for the SCNA calling.
#' * config, A list of configuration parameters.
#' * cnv.nohets, segments without heterozygous sites.
#' * crs, copy ratio data used for plotting.
#' * segs, segmentation result from GATK where seg.id has been added and
#' posterior MAFs have been converted to HDS.
#' * mafs, MINOR_ALLELE_FRACTION_POSTERIOR_50 for plotting.
#' * hets, ALT allele frequency at heterozygous sites in normal sample.
#' @export
#'

make_input <- function(input_dir, tumorID, sex, use.physical.length=FALSE) {
  if (missing(sex)) {
    sex <- "male"
  }
  if (!sex %in% c("female", "male")) {
    sex <- readr::read_lines(sex)
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
  crs <- readr::read_tsv(crs.file, comment="@", progress=FALSE, show_col_types = FALSE)
  names(crs) <- c("seqnames", "start", "end", "lrr")
  mfsegs <- readr::read_tsv(mfsegs.file, comment="@", progress=FALSE, show_col_types = FALSE) |>
    dplyr::mutate(seg.id=1:dplyr::n())
  hets <- readr::read_tsv(hets.file, comment="@", progress=FALSE, show_col_types = FALSE) |>
    dplyr::mutate(freq=ALT_COUNT / (REF_COUNT + ALT_COUNT)) |>
    dplyr::filter(!is.na(freq))

  cnv.gr <- GenomicRanges::GRanges(
    seqnames=mfsegs$CONTIG,
    IRanges::IRanges(start=mfsegs$START, end=mfsegs$END),
    seg.id=mfsegs$seg.id,
    num.mark=mfsegs$NUM_POINTS_COPY_RATIO,
    lrr10=mfsegs$LOG2_COPY_RATIO_POSTERIOR_10,
    lrr=mfsegs$LOG2_COPY_RATIO_POSTERIOR_50,
    lrr90=mfsegs$LOG2_COPY_RATIO_POSTERIOR_90,
    num.hets=mfsegs$NUM_POINTS_ALLELE_FRACTION,
    hds10=dplyr::if_else(!is.na(mfsegs$MINOR_ALLELE_FRACTION_POSTERIOR_90),
                  abs(mfsegs$MINOR_ALLELE_FRACTION_POSTERIOR_90 - 0.5), NA_real_),
    hds=dplyr::if_else(!is.na(mfsegs$MINOR_ALLELE_FRACTION_POSTERIOR_50),
                abs(mfsegs$MINOR_ALLELE_FRACTION_POSTERIOR_50 - 0.5), NA_real_),
    hds90=dplyr::if_else(!is.na(mfsegs$MINOR_ALLELE_FRACTION_POSTERIOR_10),
                  abs(mfsegs$MINOR_ALLELE_FRACTION_POSTERIOR_10 - 0.5), NA_real_)
  )
  snp.gr <- GenomicRanges::GRanges(
    seqnames=hets$CONTIG,
    IRanges::IRanges(start=hets$POSITION, end=hets$POSITION),
    freq=hets$freq
  )

  # combine cnv.gr and snp.gr
  hits <- GenomicRanges::findOverlaps(snp.gr, cnv.gr) |> tibble::as_tibble()
  cnv.df <- cnv.gr |> tibble::as_tibble() |> dplyr::rename(cnv.start=start, cnv.end=end)
  snp.df <- snp.gr |> tibble::as_tibble()
  rbd <- cbind(
    snp.df[hits$queryHits, ] |> dplyr::select(-width),
    cnv.df[hits$subjectHits, ] |>
      dplyr::select(seg.id, num.mark:lrr90, cnv.start, cnv.end, width, num.hets:hds90)
  ) |>
    tibble::as_tibble()

  # keep records of segments without heterozygous sites
  cnv.missing <- cnv.df |> dplyr::filter(!(seg.id %in% rbd$seg.id))

  # calculate HDS
  rbd <- rbd |>
    dplyr::group_by(seg.id, seqnames, cnv.start, cnv.end, width,
             num.mark, lrr10, lrr, lrr90,
             num.hets, hds10, hds, hds90) |>
    dplyr::filter(!is.na(num.mark) & !is.na(seg.id) & !is.na(freq)) |>
    dplyr::summarise(
      hds.median=median(abs(freq - 0.5)),
      hds.mad=dplyr::if_else(length(freq) > 1, mad(abs(freq - 0.5)), 0),
      hds.mean=mean(abs(freq - 0.5)),
      hds.sd=sd(abs(freq - 0.5)),
      hds.cnt=length(freq),
      .groups="drop"
    ) |>
    dplyr::rename(
      start=cnv.start, end=cnv.end
    )

  ## percentage of the length (%) in terms of the number of targets
  if (use.physical.length) {
    rbd <- rbd |>
      dplyr::mutate(width=as.numeric(end - start + 1))
    total.mark <- sum(rbd$width)
    rbd <- rbd |>
      dplyr::mutate(seg.size = width / total.mark * 100) |>
      dplyr::select(-width)
  } else {
    total.mark <- sum(rbd$num.mark, na.rm=TRUE)
    rbd <- rbd |>
      dplyr::mutate(seg.size = num.mark / total.mark * 100)
  }

  hets <- hets |>
    dplyr::rename(seqnames=CONTIG) |>
    dplyr::mutate(start=POSITION, end=POSITION) |>
    dplyr::select(seqnames, start, end, freq)

  mafs <- tibble::tibble(
    seqnames=mfsegs$CONTIG,
    start=mfsegs$START,
    end=mfsegs$END,
    seg.id=mfsegs$seg.id,
    n.hets=mfsegs$NUM_POINTS_ALLELE_FRACTION,
    maf=mfsegs$MINOR_ALLELE_FRACTION_POSTERIOR_50
  ) |> dplyr::filter(!is.na(maf))

  ## configuration paramters
  max.ploidy <- 10
  config <- list()
  config$xyGrid <- expand.grid(x=0:floor(max.ploidy/2), y=0:max.ploidy) |>
    tibble::as_tibble() |>
    dplyr::filter (x <= y & x+y <= max.ploidy) |>
    #filter( ! (x==1 & y==1) ) |>
    unique()
  # root.hdsInt and root.rInt corresponds to x=1, y=2, p=0.35
  config <- c(config,
              list(max.ploidy=max.ploidy,
                   max.ploidy.clone=max.ploidy,
                   # parameters for get_ploidy_adjustment
                   centralSegInterval = log2(1.1), # central segments: abs(lrr) < 0.1 or (abs(lrr - all.seg.mean) < 0.1)
                   uc.hds.cutoff = 0.07,     # upper central segments: subset of central segments with hds > uc.hds.cutoff
                   uc.cov.cutoff = 0.5,      # minimum coverage of upper central segments required for triploidy
                   lowHds.cutoff = 0.15,     # define low HDS segments by hds < lowHds.cutoff
                   # whether we consider triploidy and tetraploidy (get_ploidy_adjustment)
                   allow.triploidy=TRUE,
                   allow.tetraploidy=TRUE,
                   # maximum distance allowed by find_clones (changed from find_best_xyp)
                   # max.distL = 0.06,
                   # max.distR = 0.03,
                   # there are some cases where observed values are
                   # inbetween two branches.
                   # AAAABBBBB vs. AABB/AAABBB...
                   # B[k] vs. ABB/ABBB ....
                   #
                   best.tol = 0.0025,
                   # selection criteria for segments used to identify clones
                   root.hdsInt=0.075, #0.05,
                   root.rInt=log2(1.175), # 0.3,
                   min.segSize = 0.01, # 0.5, # minimum segment size used for identifying clones
                        # (this small number is necessary to detect chromothripsis and focal events)
                   min.prev=0.2, # 0.15, # this is used in find_best_xyp, find_clones, ...
                   min.prev.hom=0.4, # this is used in find_clones to exclude "noise" segments
                   # find_clones: whether we use unidentifiable ABB[.|B|BB|BBB] segments
                   useABB=TRUE,
                   #high.purity=TRUE,
                   lowest.purity=0.4,
                   # find_clones: use best xyp
                   useBestXYP=TRUE,
                   # clustering cutoff (complete linkage)
                   cutree.h=0.2,
                   #total.mark=NA,
                   #cnv.gr=NULL,
                   #max.hds.sd=0.3,
                   #min.het.cnt = 20,
                   verbose=TRUE
                   )
  )

  list(
    tumorID=tumorID,
    sex=sex,
    rbd=rbd,
    config=config,
    cnv.nohets=cnv.missing,
    crs=crs,
    segs=cnv.gr |> tibble::as_tibble(),
    mafs=mafs,
    hets=hets
  )
}

