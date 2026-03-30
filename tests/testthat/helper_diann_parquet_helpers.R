nd_write_parquet <- function(df, path) {
  arrow::write_parquet(df, sink = path)
  invisible(path)
}

nd_balanced_diann_design <- function() {
  data.frame(
    sample = c("S1", "S2", "S3", "S4", "S5", "S6"),
    group = c("A", "A", "A", "B", "B", "B"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

nd_diann_precursor_input_options <- function() {
  diannInputOptions(
    level = "precursor",
    quantityCol = "Precursor.Quantity",
    qCols = "Q.Value",
    qCutoffs = 0.01,
    rt = FALSE
  )
}

nd_make_diann_precursor_report <- function(
  design = nd_balanced_diann_design(),
  proteins = paste0("P", seq_len(20)),
  precursors_per_protein = 2L
) {
  precursor_ids <- unlist(lapply(
    proteins,
    function(protein) paste0(protein, "_pep", seq_len(precursors_per_protein))
  ))

  report <- expand.grid(
    Run = as.character(design$sample),
    `Precursor.Id` = precursor_ids,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  report$Run <- as.character(report$Run)
  report$`Precursor.Id` <- as.character(report$`Precursor.Id`)

  sample_index <- match(report$Run, design$sample)
  precursor_index <- match(report$`Precursor.Id`, precursor_ids)
  protein_group <- sub("_pep[0-9]+$", "", report$`Precursor.Id`)
  group <- design$group[sample_index]

  report$`Protein.Group` <- protein_group
  report$`Protein.Names` <- paste0("Prot_", protein_group)
  report$`Precursor.Quantity` <- 1000 +
    precursor_index * 25 +
    sample_index * 10 +
    ifelse(group == "B", 200, 0)
  report$`Q.Value` <- 0.001

  report[, c(
    "Run",
    "Precursor.Id",
    "Protein.Group",
    "Protein.Names",
    "Precursor.Quantity",
    "Q.Value"
  )]
}

nd_make_ambiguous_diann_report <- function(
  design = nd_balanced_diann_design(),
  proteins = paste0("P", seq_len(20)),
  precursors_per_protein = 2L
) {
  report <- nd_make_diann_precursor_report(
    design = design,
    proteins = proteins,
    precursors_per_protein = precursors_per_protein
  )
  report$`Precursor.Normalised` <- report$`Precursor.Quantity`

  protein_qty <- stats::aggregate(
    report$`Precursor.Quantity`,
    by = list(Run = report$Run, Protein.Group = report$`Protein.Group`),
    FUN = sum
  )
  colnames(protein_qty)[3] <- "PG.MaxLFQ"

  out <- merge(
    report,
    protein_qty,
    by = c("Run", "Protein.Group"),
    all.x = TRUE,
    sort = FALSE
  )

  out[, c(
    "Run",
    "Precursor.Id",
    "Protein.Group",
    "Protein.Names",
    "Precursor.Quantity",
    "Precursor.Normalised",
    "PG.MaxLFQ",
    "Q.Value"
  )]
}

nd_write_diann_parquet_fixture <- function(
  tmp_dir,
  design = nd_balanced_diann_design(),
  file_name = "diann_report.parquet"
) {
  data_path <- file.path(tmp_dir, file_name)
  nd_write_parquet(nd_make_diann_precursor_report(design = design), data_path)
  data_path
}

nd_load_diann_contrast_object <- function(data_path, design) {
  design_path <- withr::local_tempfile(pattern = "design_", fileext = ".tsv")
  nd_write_table(design, design_path)

  NormalyzerDE:::setupRawContrastObject(
    dataPath = data_path,
    designPath = design_path,
    sampleColName = "sample",
    inputFormat = "diann",
    inputOptions = nd_diann_precursor_input_options()
  )
}

nd_run_diann_contrast <- function(
  data_path,
  design,
  type = c("limma", "limpa")
) {
  type <- match.arg(type)
  se <- nd_load_diann_contrast_object(data_path, design)
  nst <- NormalyzerStatistics(se, logTrans = TRUE)

  args <- list(
    nst = nst,
    comparisons = "A-B",
    condCol = "group",
    type = type,
    leastRepCount = 1
  )
  if (identical(type, "limpa")) {
    args$limpaProteinIdCol <- "Protein.Group"
    args$limpaKeep <- "elist"
    args$limpaQuantArgs <- list(chunk = 10L, verbose = FALSE)
  }

  calculateContrasts(
    args$nst,
    comparisons = args$comparisons,
    condCol = args$condCol,
    type = args$type,
    leastRepCount = args$leastRepCount,
    limpaProteinIdCol = args$limpaProteinIdCol,
    limpaKeep = args$limpaKeep,
    limpaQuantArgs = args$limpaQuantArgs
  )
}

nd_extract_contrast_table <- function(nst, comparison, feature_col) {
  out <- data.frame(
    feature = as.character(annotMat(nst)[, feature_col]),
    p = unname(pairwiseCompsP(nst)[[comparison]]),
    fdr = unname(pairwiseCompsFdr(nst)[[comparison]]),
    ave = unname(pairwiseCompsAve(nst)[[comparison]]),
    fold = unname(pairwiseCompsFold(nst)[[comparison]]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  out[order(out$feature), , drop = FALSE]
}

nd_extract_limpa_elist_matrix <- function(nst, backend_key = ".global") {
  y <- backendData(nst)[["limpa"]]$elists[[backend_key]]
  genes <- as.data.frame(y$genes, check.names = FALSE)
  feature_ids <- as.character(genes$`Protein.Group`)

  mat <- y$E
  rownames(mat) <- feature_ids
  mat[order(rownames(mat)), order(colnames(mat)), drop = FALSE]
}
