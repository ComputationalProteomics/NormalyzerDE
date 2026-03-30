nd_write_table <- function(df, path, sep = "\t") {
  utils::write.table(
    df,
    file = path,
    sep = sep,
    row.names = FALSE,
    quote = FALSE
  )
  invisible(path)
}

nd_two_sample_design <- function() {
  data.frame(
    sample = c("S1", "S2"),
    group = c("A", "B"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

nd_make_design <- function(groups, sample_names = NULL) {
  groups <- as.character(groups)

  if (is.null(sample_names)) {
    sample_names <- paste0("s", seq_along(groups))
  }
  sample_names <- as.character(sample_names)

  stopifnot(length(sample_names) == length(groups))

  design <- data.frame(
    sample = sample_names,
    group = groups,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(design) <- design$sample
  design
}

nd_make_summarized_experiment <- function(
  assay,
  groups,
  row_data = NULL,
  sample_names = NULL
) {
  assay <- as.matrix(assay)

  if (is.null(sample_names)) {
    sample_names <- colnames(assay)
  }
  if (is.null(sample_names)) {
    sample_names <- paste0("s", seq_len(ncol(assay)))
  }
  sample_names <- as.character(sample_names)

  stopifnot(length(sample_names) == ncol(assay))

  colnames(assay) <- sample_names

  if (is.null(row_data)) {
    row_data <- data.frame(
      feature = paste0("f", seq_len(nrow(assay))),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  SummarizedExperiment::SummarizedExperiment(
    assay = assay,
    colData = nd_make_design(groups, sample_names = sample_names),
    rowData = row_data
  )
}

nd_write_data_and_design <- function(
  tmp_dir,
  data,
  design = nd_two_sample_design(),
  data_name = "data.tsv",
  design_name = "design.tsv",
  sep = "\t"
) {
  data_path <- file.path(tmp_dir, data_name)
  design_path <- file.path(tmp_dir, design_name)

  nd_write_table(data, data_path, sep = sep)
  nd_write_table(design, design_path)

  list(dataPath = data_path, designPath = design_path)
}
