context("limpa option helpers")

test_that("limpaOptions and normalizeLimpaOptionsInput validate structure", {
  opts <- limpaOptions(
    proteinIdCol = "protein_id",
    quantArgs = list(chunk = 25L),
    keep = "elist"
  )

  expect_s3_class(opts, "normalyzerde_limpa_options")
  expect_equal(opts$proteinIdCol, "protein_id")
  expect_equal(opts$quantArgs$chunk, 25L)
  expect_equal(opts$keep, "elist")

  expect_error(limpaOptions(unknown = TRUE), "unused argument")
  expect_error(
    NormalyzerDE:::normalizeLimpaOptionsInput("bad"),
    class = "normalyzerde_error"
  )
  expect_error(
    NormalyzerDE:::normalizeLimpaOptionsInput(list(TRUE)),
    class = "normalyzerde_error"
  )
  expect_error(
    NormalyzerDE:::normalizeLimpaOptionsInput(list(unknown = TRUE)),
    class = "normalyzerde_error"
  )

  empty_opts <- NormalyzerDE:::normalizeLimpaOptionsInput(list())
  expect_s3_class(empty_opts, "normalyzerde_limpa_options")
})

test_that("limpa option builders validate scalar inputs", {
  expect_error(
    NormalyzerDE:::buildLimpaOptionsInternal(proteinIdCol = ""),
    class = "normalyzerde_error"
  )
  expect_error(
    NormalyzerDE:::buildLimpaOptionsInternal(byRow = NA),
    class = "normalyzerde_error"
  )
  expect_error(
    NormalyzerDE:::buildLimpaOptionsInternal(quantifiedRds = ""),
    class = "normalyzerde_error"
  )
})

test_that("limpa arg sanitizers require named lists and drop forbidden args", {
  expect_equal(NormalyzerDE:::sanitizeLimpaDpcArgs(NULL), list())
  expect_error(
    NormalyzerDE:::sanitizeLimpaDpcArgs(1),
    class = "normalyzerde_error"
  )
  expect_error(
    NormalyzerDE:::sanitizeLimpaDpcArgs(list(1)),
    class = "normalyzerde_error"
  )
  expect_equal(
    NormalyzerDE:::sanitizeLimpaDpcArgs(list(y = 1, keep = "x")),
    list(keep = "x")
  )

  expect_equal(NormalyzerDE:::sanitizeLimpaQuantArgs(NULL), list())
  expect_error(
    NormalyzerDE:::sanitizeLimpaQuantArgs(1),
    class = "normalyzerde_error"
  )
  expect_error(
    NormalyzerDE:::sanitizeLimpaQuantArgs(list(1)),
    class = "normalyzerde_error"
  )
  expect_equal(
    NormalyzerDE:::sanitizeLimpaQuantArgs(
      list(y = 1, protein.id = "p", dpc = c(0, 1), chunk = 10L)
    ),
    list(chunk = 10L)
  )

  expect_equal(NormalyzerDE:::sanitizeLimpaDEArgs(NULL), list())
  expect_error(
    NormalyzerDE:::sanitizeLimpaDEArgs(1),
    class = "normalyzerde_error"
  )
  expect_error(
    NormalyzerDE:::sanitizeLimpaDEArgs(list(1)),
    class = "normalyzerde_error"
  )
  expect_equal(
    NormalyzerDE:::sanitizeLimpaDEArgs(
      list(y = 1, design = "x", plot = TRUE, robust = TRUE)
    ),
    list(robust = TRUE)
  )
})

test_that("limpa quant and DE defaults are applied and validated", {
  quant_defaults <- NormalyzerDE:::applyLimpaQuantDefaultsAndValidate(list())
  expect_equal(quant_defaults[["dpc.slope"]], 0.8)
  expect_equal(quant_defaults[["chunk"]], 1000L)
  expect_false(quant_defaults[["verbose"]])

  expect_error(
    NormalyzerDE:::applyLimpaQuantDefaultsAndValidate(list("dpc.slope" = 0)),
    class = "normalyzerde_error"
  )
  expect_error(
    NormalyzerDE:::applyLimpaQuantDefaultsAndValidate(list(chunk = 0L)),
    class = "normalyzerde_error"
  )
  expect_error(
    NormalyzerDE:::applyLimpaQuantDefaultsAndValidate(list(verbose = NA)),
    class = "normalyzerde_error"
  )

  de_defaults <- NormalyzerDE:::applyLimpaDEDefaultsAndValidate(list())
  expect_false(de_defaults[["sample.weights"]])
  expect_error(
    NormalyzerDE:::applyLimpaDEDefaultsAndValidate(list("sample.weights" = NA)),
    class = "normalyzerde_error"
  )
})

test_that("limpa quantified helpers normalize metadata and validate protein IDs", {
  y_quant <- list(
    E = matrix(c(1, 2, 3, 4), nrow = 2, dimnames = list(NULL, c("s1", "s2"))),
    other = list(
      n.observations = matrix(1, nrow = 2, ncol = 2),
      standard.error = matrix(0.1, nrow = 2, ncol = 2)
    ),
    genes = data.frame(extra = c("x", "y"), check.names = FALSE)
  )

  normalized <- NormalyzerDE:::normalizeLimpaQuantifiedEList(
    y_quant,
    proteinIdCol = "protein",
    proteinIds = c("p1", "p2")
  )

  expect_equal(rownames(normalized$E), c("1", "2"))
  expect_equal(rownames(normalized$other$n.observations), c("1", "2"))
  expect_equal(rownames(normalized$other$standard.error), c("1", "2"))
  expect_equal(colnames(normalized$genes)[1], "protein")
  expect_equal(normalized$genes$protein, c("p1", "p2"))
  expect_equal(rownames(normalized$genes), c("1", "2"))

  expect_error(
    NormalyzerDE:::quantifyLimpaByProteinInternal(
      dataMat = matrix(seq_len(6), nrow = 3),
      genesDf = data.frame(protein = c("p1", "p2"), stringsAsFactors = FALSE),
      proteinIdCol = "protein",
      dpc = c(0, 1),
      quantArgs = list()
    ),
    class = "normalyzerde_error"
  )

  expect_error(
    NormalyzerDE:::quantifyLimpaByProteinInternal(
      dataMat = matrix(seq_len(4), nrow = 2),
      genesDf = data.frame(
        protein = c("p1", ""),
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),
      proteinIdCol = "protein",
      dpc = c(0, 1),
      quantArgs = list()
    ),
    class = "normalyzerde_error"
  )
})

test_that("limpa sample weight helpers collect stored and derived weights", {
  se <- nd_make_summarized_experiment(
    assay = matrix(c(1, 2, 3, 4), nrow = 2),
    groups = c("A", "B")
  )
  nst <- NormalyzerStatistics(se, logTrans = FALSE)

  backendData(nst) <- list(
    limpa = list(
      sampleWeights = list("A-B" = c(s1 = 0.5, s2 = 0.8))
    )
  )

  stored <- NormalyzerDE:::collectLimpaSampleWeights(nst)
  expect_equal(stored$comparison, c("A-B", "A-B"))
  expect_equal(stored$sample, c("s1", "s2"))
  expect_equal(stored$sampleWeight, c(0.5, 0.8))

  filtered <- NormalyzerDE:::collectLimpaSampleWeights(nst, comparison = "missing")
  expect_s3_class(filtered, "data.frame")
  expect_equal(nrow(filtered), 0)

  fit <- list(
    targets = data.frame(
      sample.weight = c(0.25, 0.75),
      row.names = c("s1", "s2"),
      check.names = FALSE
    ),
    EList = list(E = matrix(c(1, 2), nrow = 1, dimnames = list(NULL, c("s1", "s2"))))
  )

  backendData(nst) <- list(
    limpa = list(
      fits = list(".global" = fit)
    )
  )

  derived <- NormalyzerDE:::collectLimpaSampleWeights(nst)
  expect_equal(derived$comparison, c(".global", ".global"))
  expect_equal(derived$sample, c("s1", "s2"))
  expect_equal(derived$sampleWeight, c(0.25, 0.75))

  expect_null(NormalyzerDE:::extractLimpaSampleWeightsFromFit(NULL))
  expect_error(getLimpaSampleWeights(list()), class = "normalyzerde_error")

  backendData(nst) <- list()
  expect_error(getLimpaSampleWeights(nst), class = "normalyzerde_error")
})

test_that("limpa row/protein helpers normalize configuration and infer metadata", {
  expect_equal(
    NormalyzerDE:::normalizeLimpaByRowConfig(TRUE, "auto"),
    list(limpaByRow = TRUE, limpaProteinIdCol = NULL)
  )
  expect_equal(
    NormalyzerDE:::normalizeLimpaByRowConfig(FALSE, "Protein.Group"),
    list(limpaByRow = FALSE, limpaProteinIdCol = "Protein.Group")
  )
  expect_error(
    NormalyzerDE:::normalizeLimpaByRowConfig(NA, "auto"),
    class = "normalyzerde_error"
  )
  expect_error(
    NormalyzerDE:::normalizeLimpaByRowConfig(TRUE, "Protein.Group"),
    class = "normalyzerde_error"
  )

  annot <- as.matrix(data.frame(
    `Protein.Group` = c("P1", "P2"),
    Other = c("x", "y"),
    check.names = FALSE
  ))
  expect_equal(
    NormalyzerDE:::inferLimpaProteinIdCol(annot, "auto"),
    "Protein.Group"
  )
  expect_null(NormalyzerDE:::inferLimpaProteinIdCol(annot, NULL))
  expect_error(
    NormalyzerDE:::inferLimpaProteinIdCol(annot, "Missing"),
    class = "normalyzerde_error"
  )
  expect_null(
    NormalyzerDE:::inferLimpaProteinIdCol(
      as.matrix(data.frame(Feature = c("f1", "f2"), check.names = FALSE)),
      "auto"
    )
  )

  genes_df <- data.frame(
    protein = c("p1", "p1", "p2"),
    stable = c("A", "A", "B"),
    varying = c("x", "y", "z"),
    empty = c("", "", ""),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  expect_equal(
    NormalyzerDE:::inferStableProteinAnnotationCols(
      genes_df,
      proteinId = genes_df$protein,
      proteinIdCol = "protein"
    ),
    "stable"
  )
})

test_that("limpa DPC estimation helpers validate argument names and none mode", {
  dummy_fn <- function(y, allowed = NULL) NULL

  expect_equal(
    NormalyzerDE:::validateLimpaArgsByFormals(
      list(allowed = 1),
      fn = dummy_fn,
      methodLabel = "dpc"
    ),
    list(allowed = 1)
  )
  expect_error(
    NormalyzerDE:::validateLimpaArgsByFormals(
      list(bad = 1),
      fn = dummy_fn,
      methodLabel = "dpc"
    ),
    class = "normalyzerde_error"
  )

  expect_null(
    NormalyzerDE:::estimateLimpaDpcFromData(
      dataMat = matrix(1, nrow = 1),
      limpaDpcMethod = "none"
    )
  )
})

test_that("limpa helper leftovers handle empty names and optional metadata", {
  se <- nd_make_summarized_experiment(
    assay = matrix(c(1, 2, 3, 4), nrow = 2),
    groups = c("A", "B")
  )
  nst <- NormalyzerStatistics(se, logTrans = FALSE)

  expect_error(
    NormalyzerDE:::validateLimpaDpc("bad"),
    class = "normalyzerde_error"
  )

  backendData(nst) <- list(limpa = list(sampleWeights = list()))
  expect_null(NormalyzerDE:::collectLimpaSampleWeights(nst))

  backendData(nst) <- list(limpa = list(sampleWeights = list(c(0.2, 0.8))))
  unnamed_weights <- NormalyzerDE:::collectLimpaSampleWeights(nst)
  expect_equal(unnamed_weights$comparison, c(".global", ".global"))
  expect_equal(unnamed_weights$sample, c("1", "2"))
  expect_equal(unnamed_weights$sampleWeight, c(0.2, 0.8))

  fit_named <- list(
    targets = list(sample.weight = c(s1 = 0.3, s2 = 0.7))
  )
  expect_equal(
    NormalyzerDE:::extractLimpaSampleWeightsFromFit(fit_named),
    c(s1 = 0.3, s2 = 0.7)
  )

  fit_rows <- list(
    targets = data.frame(
      sample.weight = c(0.4, 0.6),
      row.names = c("s1", "s2"),
      check.names = FALSE
    )
  )
  expect_equal(
    NormalyzerDE:::extractLimpaSampleWeightsFromFit(fit_rows),
    c(s1 = 0.4, s2 = 0.6)
  )

  expect_null(
    NormalyzerDE:::inferLimpaProteinIdCol(
      matrix("x", nrow = 2, ncol = 1),
      "auto"
    )
  )
  expect_equal(
    NormalyzerDE:::inferStableProteinAnnotationCols(
      data.frame(check.names = FALSE),
      proteinId = character(),
      proteinIdCol = "protein"
    ),
    character()
  )

  expect_equal(
    NormalyzerDE:::normalizeLimpaByRowConfig(TRUE, NULL),
    list(limpaByRow = TRUE, limpaProteinIdCol = NULL)
  )

  normalized <- NormalyzerDE:::normalizeLimpaQuantifiedEList(
    list(E = matrix(c(1, 2), nrow = 1), other = list()),
    proteinIdCol = NULL,
    proteinIds = NULL
  )
  expect_null(normalized$genes)
})
