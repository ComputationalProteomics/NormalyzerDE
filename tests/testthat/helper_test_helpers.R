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

nd_limpa_design <- function() {
  nd_make_design(c("A", "A", "A", "B", "B", "B"))
}

nd_limpa_de_matrix <- function() {
  mat <- matrix(
    c(
      10,
      11,
      10,
      13,
      12,
      11,
      NA,
      NA,
      NA,
      9,
      9,
      10,
      5,
      5,
      5,
      NA,
      NA,
      NA,
      7,
      8,
      7,
      7,
      7,
      7
    ),
    nrow = 4,
    byrow = TRUE
  )
  colnames(mat) <- paste0("s", seq_len(ncol(mat)))
  mat
}

nd_limpa_input_data <- function(mat = nd_limpa_de_matrix()) {
  data.frame(
    feature = paste0("f", seq_len(nrow(mat))),
    as.data.frame(mat, check.names = FALSE),
    check.names = FALSE
  )
}

nd_write_limpa_input_fixture <- function(
  tmp_dir,
  mat = nd_limpa_de_matrix(),
  design = nd_limpa_design(),
  data_name = "limpa_data.tsv",
  design_name = "limpa_design.tsv"
) {
  nd_write_data_and_design(
    tmp_dir = tmp_dir,
    data = nd_limpa_input_data(mat),
    design = design,
    data_name = data_name,
    design_name = design_name
  )
}

nd_make_limpa_prequant_fixture <- function(
  n_features = 20L,
  groups = c("A", "A", "B", "B"),
  row_data = data.frame(
    `Protein.Group` = paste0("P", seq_len(n_features)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ),
  seed = 1,
  missing_count = 8L,
  mean = 10,
  sd = 1
) {
  set.seed(seed)

  mat <- matrix(
    stats::rnorm(n_features * length(groups), mean = mean, sd = sd),
    nrow = n_features
  )
  colnames(mat) <- paste0("s", seq_len(ncol(mat)))

  if (missing_count > 0) {
    mat[sample.int(length(mat), missing_count)] <- NA_real_
  }

  stopifnot(nrow(row_data) == nrow(mat))

  list(
    mat = mat,
    se = nd_make_summarized_experiment(
      assay = mat,
      groups = groups,
      row_data = row_data
    )
  )
}

nd_run_prequant_limpa <- function(
  se,
  out_dir,
  job_name,
  limpa_options,
  quiet = TRUE,
  suppress_warnings = TRUE
) {
  run <- normalyzer(
    jobName = job_name,
    experimentObj = se,
    outputDir = out_dir,
    preQuant = "limpa",
    limpaOptions = limpa_options,
    noLogTransform = TRUE,
    normalizeRetentionTime = FALSE,
    skipAnalysis = TRUE,
    quiet = quiet,
    sampleAbundThres = 1,
    requireReplicates = FALSE
  )

  if (isTRUE(suppress_warnings)) {
    suppressWarnings(run)
  } else {
    run
  }
}

nd_run_limpa_de_from_paths <- function(
  paths,
  out_dir,
  job_name,
  limpa_options,
  quiet = TRUE
) {
  suppressWarnings(normalyzerDE(
    jobName = job_name,
    comparisons = "A-B",
    designPath = paths$designPath,
    dataPath = paths$dataPath,
    outputDir = out_dir,
    type = "limpa",
    logTrans = FALSE,
    leastRepCount = 1,
    limpaOptions = limpa_options,
    quiet = quiet
  ))
}

nd_capture_conditions <- function(expr) {
  messages <- character()
  warnings <- character()
  error <- NULL

  value <- tryCatch(
    withCallingHandlers(
      expr,
      message = function(m) {
        messages <<- c(messages, conditionMessage(m))
        invokeRestart("muffleMessage")
      },
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      error <<- e
      NULL
    }
  )

  list(
    value = value,
    error = error,
    messages = messages,
    warnings = warnings
  )
}

nd_expect_messages <- function(messages, patterns) {
  for (pattern in patterns) {
    testthat::expect_true(
      any(grepl(pattern, messages)),
      info = paste("Missing message matching:", pattern)
    )
  }
}

nd_expect_error_cases <- function(cases, class = "normalyzerde_error") {
  for (case_name in names(cases)) {
    testthat::expect_error(
      cases[[case_name]](),
      class = class,
      info = case_name
    )
  }
}

nd_with_pdf <- function(path, code) {
  grDevices::pdf(path)
  device_id <- grDevices::dev.cur()

  on.exit(
    {
      open_devices <- grDevices::dev.list()
      if (!is.null(open_devices) && device_id %in% open_devices) {
        grDevices::dev.off(device_id)
      }
    },
    add = TRUE
  )

  force(code)
}
