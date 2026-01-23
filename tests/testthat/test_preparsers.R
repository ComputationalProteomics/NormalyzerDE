context("preparsers")

test_that("proteiosToNormalyzer reads a Proteios matrix", {
  fp <- tempfile(fileext = ".tsv")
  on.exit(unlink(fp), add = TRUE)

  mat <- matrix(c(1, 2, 3, 4), nrow = 2)
  utils::write.table(
    mat,
    file = fp,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  out <- NormalyzerDE:::proteiosToNormalyzer(fp)
  expect_true(is.matrix(out))
  expect_equal(dim(out), dim(mat))
  expect_equal(as.numeric(out), as.numeric(mat))
})

test_that("maxQuantToNormalyzer parses MaxQuant peptide-level tables", {
  fp <- tempfile(fileext = ".txt")
  on.exit(unlink(fp), add = TRUE)

  in_df <- data.frame(
    Sequence = c("AAA", "BBB"),
    Mass = c(100.1, 200.2),
    Proteins = c("P1", "P2"),
    Leading.razor.protein = c("P1", "P2"),
    PEP = c(0.01, 0.02),
    Charges = c(2, 3),
    Intensity.S1 = c(1000, 2000),
    Intensity.S2 = c(1100, 2100),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  utils::write.table(
    in_df,
    file = fp,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  out <- NormalyzerDE:::maxQuantToNormalyzer(fp, protLevel = FALSE, sep = "\t")

  expect_true(is.matrix(out))
  expect_equal(nrow(out), nrow(in_df) + 1)
  expect_equal(ncol(out), 6 + 2)
  expect_equal(
    unname(out[1, ]),
    c(
      "Sequence",
      "Mass",
      "Proteins",
      "Leading.razor.protein",
      "PEP",
      "Charges",
      "S1",
      "S2"
    )
  )
  expect_equal(unname(out[2, 1]), "AAA")
})

test_that("maxQuantToNormalyzer parses MaxQuant protein-level tables", {
  fp <- tempfile(fileext = ".txt")
  on.exit(unlink(fp), add = TRUE)

  in_df <- data.frame(
    Protein.IDs = c("P1", "P2"),
    Majority.protein.IDs = c("P1", "P2"),
    Fasta.headers = c("hdr1", "hdr2"),
    Intensity.S1 = c(1000, 2000),
    Intensity.S2 = c(1100, 2100),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  utils::write.table(
    in_df,
    file = fp,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  out <- NormalyzerDE:::maxQuantToNormalyzer(fp, protLevel = TRUE, sep = "\t")

  expect_true(is.matrix(out))
  expect_equal(nrow(out), nrow(in_df) + 1)
  expect_equal(ncol(out), 3 + 2)
  expect_equal(
    unname(out[1, ]),
    c(
      "Protein.IDs",
      "Majority.protein.IDs",
      "Fasta.headers",
      "S1",
      "S2"
    )
  )
})

test_that("maxQuantToNormalyzer errors when required columns are missing", {
  fp <- tempfile(fileext = ".txt")
  on.exit(unlink(fp), add = TRUE)

  in_df <- data.frame(
    Mass = c(100.1),
    Proteins = c("P1"),
    Leading.razor.protein = c("P1"),
    PEP = c(0.01),
    Charges = c(2),
    Intensity.S1 = c(1000),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  utils::write.table(
    in_df,
    file = fp,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  expect_error(
    NormalyzerDE:::maxQuantToNormalyzer(fp, protLevel = FALSE, sep = "\t"),
    "Didn't find all of the following expected columns"
  )
})
