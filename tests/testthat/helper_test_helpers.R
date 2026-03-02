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
