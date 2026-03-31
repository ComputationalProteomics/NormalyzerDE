context("analyzeResults helper functions")

test_that("calculateReplicateCV reshapes single-group results consistently", {
  methodList <- list(
    raw = matrix(c(1, 2, 3, 2, 4, 6), nrow = 2, byrow = TRUE),
    scaled = matrix(c(2, 4, 6, 1, 3, 5), nrow = 2, byrow = TRUE)
  )
  groups <- rep("A", 3)

  out <- calculateReplicateCV(methodList, groups)

  expected <- vapply(
    methodList,
    function(methodData) {
      mean(apply(
        methodData,
        1,
        function(feature) {
          stats::sd(feature) / mean(feature)
        }
      )) * 100
    },
    0
  )

  expect_equal(dim(out), c(1L, length(methodList)))
  expect_equal(unname(out[1, ]), unname(expected))
  expect_equal(colnames(out), names(methodList))
})

test_that("calculateSummarizedCorrelationVector validates correlation type", {
  methodList <- list(raw = matrix(c(1, 2, 3, 4), nrow = 2))

  expect_error(
    calculateSummarizedCorrelationVector(
      methodlist = methodList,
      allReplicateGroups = c("A", "A"),
      sampleGroupsWithReplicates = "A",
      corrType = "kendall"
    ),
    class = "normalyzerde_error"
  )
})

test_that("calculateANOVAPValues supports categorical group labels", {
  methodList <- list(
    raw = matrix(
      c(
        1.0, 1.1, 2.0, 2.1, 5.0, 5.1,
        4.0, 4.1, 4.0, 4.1, 4.0, 4.1
      ),
      nrow = 2,
      byrow = TRUE
    )
  )
  groups <- c("low", "low", "mid", "mid", "high", "high")

  expected <- apply(
    methodList$raw,
    1,
    function(sampleIndex) {
      summary(stats::aov(unlist(sampleIndex) ~ factor(groups)))[[1]][[5]][1]
    }
  )

  out <- calculateANOVAPValues(methodList, groups, categoricalANOVA = TRUE)

  expect_equal(unname(as.vector(out)), unname(expected))
})

test_that("findLowlyVariableFeaturesCVs handles warning and mismatch paths", {
  methodList <- list(raw = matrix(c(1, 2, 3, 4), nrow = 2))

  expect_warning(
    expect_null(findLowlyVariableFeaturesCVs(c(Inf, Inf), methodList)),
    class = "normalyzerde_warning"
  )

  expect_error(
    findLowlyVariableFeaturesCVs(
      seq(0.01, 0.20, by = 0.01),
      list(raw = matrix(seq_len(9), nrow = 3))
    ),
    class = "normalyzerde_error"
  )
})
