#' Write GWAS results to disk.
#'
#' @param res A results list returned by \code{\link{farmcpu}}.
#'
#' @author Aaron Kusmec
#'
#' @export

write_results <- function(res) {
  for (i in names(res)) {
    for (j in names(res[[i]])) {
      write.csv(res[[i]][[j]], file = paste0(i, ".", j, ".csv"),
                quote = FALSE, row.names = FALSE)
    }
  }
}
