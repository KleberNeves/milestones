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


#' Builds queries from cited references
#'
#' For a set of Web of Science records, this function builds and saves
#' text files containing queries that can be pasted on Web of Science to
#' obtain the records for the references cited by the first set of records
#' (i.e. the next generation of cited references). The query is built twice,
#' using DOIs and then author names and journals. These are redundant to
#' obtain a better recall in the search.
#'
#' @param M A bibliometrix M data frame.
#' @param next_name The name you want to give to the generation (e.g. G1).
#' @return A list with the cited reference years, counted as indicated.
#' @export
make_wos_query <- function (M, next_name) {
  k = get_CR_dois(M)
  k = k[!is.na(k)]
  k = unique(k)
  parts = split(k, ceiling(seq_along(k)/3000))

  k2 = get_CR_names(M)
  k2 = k2[!is.na(k2)]
  k2 = unique(k2)
  parts2 = split(k2, ceiling(seq_along(k2)/500))

  save.query = function (q, ps, dois = T) {
    part = ps[[q]]

    if (dois) {
      query = paste0('DO=("',
                     paste(part, collapse = '" OR "'),
                     '")')
      f = file(paste0("./query_", next_name, "_DO_", q,".txt"))
    } else {
      query = paste(part, collapse = ' OR ')
      f = file(paste0("./query_", next_name, "_AU_", q,".txt"))
    }

    writeLines(query, f)
    close(f)
  }

  if (length(parts) > 0)
    lapply(1:length(parts), save.query, ps = parts, dois = T)

  if (length(parts2) > 0)
    lapply(1:length(parts2), save.query, ps = parts2, dois = F)
}

#' Builds parts of queries
#'
#' Helper function for making Web of Science queries using
#' the DOIs of the cited references.
#'
#' @param M A bibliometrix M data frame.
#' @return A character vector with the queries for the cited references.
#' @importFrom magrittr %>%
get_CR_dois = function (M) {
  k = stringr::str_sub(
    unlist(
      stringr::str_extract_all((M %>% dplyr::arrange(-TC))$CR, "DOI 10(.+?);")
    ), start = 5, end = -2)

  k = stringr::str_remove(k, "#")
  k = stringr::str_remove(k, "\\($")
  k = stringr::str_remove(k, "\\)$")
  k = stringr::str_remove(k, ", .+")

  k
}

#' Builds parts of queries
#'
#' Helper function for making Web of Science queries using
#' the author, journal and year of publication of the cited references.
#'
#' @param M A bibliometrix M data frame.
#' @return A character vector with the queries for the cited references.
#' @importFrom magrittr %>%
get_CR_names = function (M) {
  k = unlist(stringr::str_split((M %>% dplyr::arrange(-TC))$CR, "; "))
  k = k[which(!stringr::str_detect(k, "DOI 10"))]
  k = stringr::str_split(k, ", ")
  if (length(k) > 0) {
    k = do.call(rbind, k)[,1:3]
    k = apply(matrix(k, ncol = 3), 1, function(x) {
      # browser()
      if (x[1] == "[ANONYMOUS]") { return ("") }
      if (!is.na(as.numeric(x[1]))) { return ("") }
      if (!stringr::str_detect(x[2], "^[0-9]{4}$")) {
        x[2] = x[3]
        x[2] = NA
      }
      x = stringr::str_remove_all(x, "[.?+(),;:/'*’~^´`&-‐_=‘—“”]")
      x = stringr::str_remove_all(x, '"')
      if (is.na(as.numeric(x[2]))) {
        paste0('(AU=("',
               stringr::str_replace_all(x[1], "\\.", ""),
               '") AND SO=("',
               x[3],
               '"))')
      } else {
        paste0('(AU=("',
               stringr::str_replace_all(x[1], "\\.", ""),
               '") AND PY=("',
               x[2],
               '") AND SO=("',
               x[3],
               '"))')
      }
    })
    k = k[k != ""]
  }
  k
}
