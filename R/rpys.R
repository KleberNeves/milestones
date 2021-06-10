#' Identifies milestones through RPYS papers
#'
#' Given a bibliometric data frame, it identifies the cited reference years
#' and builds the RPYS graph. It then finds the peak years and extracts the
#' most cited references in each of those peaks. These are the documents
#' returned with the graph (the graph has the peaks highlighted).
#'
#' @param M A bibliometric data frame (from bibliometrix).
#' @return A list containing a data frame of the documents and the RPYS plot.
#' @importFrom magrittr %>%
#' @export
identify_rpys_milestones <- function (M) {

  # Clean the references, keep only the well-formatted ones
  print ("Processing references ...")
  cited_refs = tibble::tibble(
    REF = M$CR %>% stringr::str_split(";") %>% unlist(),
    CLEAN_REF = REF %>% purrr::map_chr(~stringr::str_extract(.x,'[A-Z .-]+, [0-9]{4}, .+?(,|$)')) %>%
      stringr::str_remove_all("[.]") %>%
      stringr::str_remove(",$"),
  ) %>%

  # Get all the years of the references cited
    dplyr::mutate(
      AUTHOR = stringr::str_extract(CLEAN_REF,'[A-Z .-]+,') %>% stringr::str_remove(',$'),
      YEAR = CLEAN_REF %>%
        stringr::str_extract("(, [0-9]{4},)|(^[0-9]{4})") %>%
        stringr::str_extract("[0-9]{4}") %>%
        as.numeric()
    )

  # If a year has more citations than the average of the previous 5 years, consider it a peak
  is_peak = function (variable) {
    local_speed = purrr::imap_dbl(variable, function(current, i) {
      if (i > 1) before = variable[i - 1] else before = 0
      current - before
    })

    downwards = purrr::imap_lgl(local_speed, function(current, i) {
      if (i > 1) before = local_speed[i - 1] else before = 0
      if (i < length(local_speed)) after = local_speed[i + 1] else after = 0
      current > 0 & after < 0
    })

    downwards
  }

  lvls = (min(cited_refs$YEAR, na.rm = T) - 1):(max(cited_refs$YEAR, na.rm = T))

  cited_years = cited_refs %>%
    dplyr::filter(!is.na(YEAR)) %>%
    dplyr::mutate(YEAR = factor(YEAR, levels = lvls)) %>%
    dplyr::count(YEAR, .drop = F) %>%
    dplyr::mutate(
      peak = is_peak(n)
    )

  peak_years = cited_years %>% dplyr::filter(peak) %>% dplyr::pull(YEAR) %>%
    as.character() %>% as.numeric()

  # For each peak, pick the most-cited paper of each year
  peak_papers = map_dfr(peak_years, function(peak_year) {
    year_refs = cited_refs %>%
      dplyr::filter(YEAR == peak_year)

    if (nrow(year_refs) == 0) return (tibble::tibble())

    top_authors = year_refs %>%
      dplyr::count(AUTHOR) %>%
      dplyr::slice_max(order_by = n, n = 1)

    top_refs = year_refs %>%
      dplyr::right_join(top_authors, by = "AUTHOR") %>%
      dplyr::distinct(CLEAN_REF, .keep_all = T) %>%
      dplyr::select(REF, YEAR)

    top_refs
  })

  # Return a data frame with these papers and a RPYS graph
  p = ggplot2::ggplot(cited_years %>% dplyr::filter(n > 1)) +
    ggplot2::aes(x = YEAR, y = n, fill = peak) +
    ggplot2::geom_col() +
    ggplot2::theme(legend.position = "none")

  list(documents = peak_papers, rpys = p)
}

#' Save list of RPYS milestones papers
#'
#' Saves the list of papers identified through the RPYS algorithm
#'
#' @param RPYS_MS The object returned by the *identify_rpys_milestones* function.
#' @param filename Filename to save the list of milestone references.
#' @importFrom magrittr %>%
#' @importFrom utils write.table
#' @export
save_papers_rpys <- function (RPYS_MS, filename) {
  x = RPYS_MS$documents
  colnames(x) = c("Reference", "Year")
  write.table(x, filename, sep = "\t", row.names = F)
}
