#' Load a bibliometric data frame
#'
#' Obtains the "M" data frame using bibliometrix to read a vector of files
#' or all text files in the given folder if the argument path is given.
#' Files must be those downloaded from Web of Science, in the plain text format.
#'
#' @param filenames A vector with the filenames.
#' @param path Full path to the folder containing the bibliometric data.
#' @return An "M" data frame, as loaded and converted through bibliometrix.
#' @export
read_biblio_data <- function (filenames, path = NULL) {
  if (!is.null(path)) {
    filenames = paste0(folder, "/", list.files(folder, ".txt$"))
  }
  M = bibliometrix::convert2df(filenames, dbsource = "isi", format = "plaintext")
  M
}

make_wos_query <- function (M, next_name) {

}
