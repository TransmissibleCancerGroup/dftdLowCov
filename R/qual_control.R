###
# Common functions for quality control scripts


#' Easier extraction of data from data table columns using regex
#' @export
scan_col <- function(pattern, s, dtypes = NULL, names = NULL) {
    search <- regexec(pattern, s)
    find <- regmatches(s, search)
    extraction_cols <- 2:length(find[[1]])
    extracted <- lapply(extraction_cols, function(i) sapply(regmatches(s, search), `[`, i))
    if (!is.null(dtypes)) {
        for (n in 1:length(dtypes)) {
            if (dtypes[n] == 'i') extracted[[n]] <- as.integer(extracted[[n]])
            else if (dtypes[n] == 'c') extracted[[n]] <- as.character(extracted[[n]])
            else if (dtypes[n] == 'C') extracted[[n]] <- toupper(as.character(extracted[[n]]))
            else if (dtypes[n] == 'n') extracted[[n]] <- as.numeric(extracted[[n]])
            else if (dtypes[n] == 'f') extracted[[n]] <- as.factor(extracted[[n]])
        }
    }
    names(extracted) <- names
    extracted
}

#' parse sample cnv group and count
#' @importFrom "data.table" as.data.table
#' @export
parse_sample_cnv_group_and_count <- function(s) {
    p <- strsplit(s, ",")[[1]]
    pattern <- "([0-9Tab.]*) \\(.*\\) [>=]*([+-]?[0-9]*[.]?[0-9]+)"
    dt <- as.data.table(scan_col(pattern, p, c('c', 'n'), names = c("sample", "copy_number")))
    dt[, .(count = .N), by = copy_number][order(copy_number)]
}

#' parse_sample_cnv_group
#' @export
parse_sample_cnv_group <- function(s) {
    p <- strsplit(s, ",")[[1]]
    pattern <- "([0-9Tab.]*) \\(.*\\) [>=]*([+-]?[0-9]*[.]?[0-9]+)"
    dt <- as.data.table(scan_col(pattern, p, c('c', 'n'), names = c("sample", "copy_number")))
}

#' parse_sample_cnv_group_with_clade
#' Also capture the clade group
#' @export
parse_sample_cnv_group_with_clade <- function(s) {
    p <- strsplit(s, ",")[[1]]
    pattern <- "([0-9Tab.]*) (\\(.*\\)) [>=]*([+-]?[0-9]*[.]?[0-9]+)"
    dt <- as.data.table(scan_col(pattern, p, c('c', 'c', 'n'), names = c("sample", "clade", "copy_number")))
}

#' parse_aberrant_copy_number_states
#' @export
parse_aberrant_copy_number_states <- function(s) {
    p <- strsplit(s, ",")[[1]]
    pattern <- "[>=]*([+-]?[0-9]*[.]?[0-9]+) \\((\\d+)\\)"
    as.data.table(scan_col(pattern, p, c('n', 'n'), names = c("copy_number", "count")))
}
