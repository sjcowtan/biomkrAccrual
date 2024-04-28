#' Display arm closure summaries
#' @param file_prefix Consistent beginning of filename
#' @param output_path Directory where output .csvs are written
#' @export
ggclosures <- function(
  file_prefix = "arm_closures_24-04-28_",
  output_path = "output_data/"
) {
  # What output files do we have?
  filenames <- list.files(
    output_path, 
    pattern = paste0("^", file_prefix, ".*.csv"), 
    full.names = TRUE
  )

  # Read files
  closures_ls <- lapply(filenames, read.csv)

  # Summarise files
  summ <- lapply(closures_ls, summary)
  name_types <- c("1", "2", "mixed", "unbalanced")
  names(summ) <- c(
    paste0("gamma_rate_closures_", name_types),
    paste0("fixed_rate_closures_", name_types)
  )

  print(summ)

  # Standard Deviations
  sd_mx <- t(sapply(
    closures_ls, 
    function(a_df) {
      sapply(
        seq_len(ncol(a_df)),
        function(i) sd(a_df[i, ], na.rm = TRUE)
      )
    }
  ))

  sd_mx <- sd_mx[, -1]

  colnames(sd_mx) <- paste0("T", 1:3, "_sd")
  rownames(sd_mx) <- c(
   paste0("gamma_rate_closures_", name_types),
    paste0("fixed_rate_closures_", name_types)
  )

  print(sd_mx)
}