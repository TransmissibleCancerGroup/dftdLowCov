# Functions for counting breakpoint reuse
#' @export
make_segmatch <- function(coordinates, ignore_direction = TRUE) {

    # Check the data has the expected columns
    expected_columns <- c("Order", "chr", "start", "end", "Type", "CNV_ID_ext", "Category")
    if (!all(expected_columns %in% colnames(coordinates))) {
        missing_columns <- paste(setdiff(expected_columns, colnames(coordinates)))
        stop(paste("Input data is missing the following necessary columns:", missing_columns))
    }

    # "by_columns" = The columns to use as data.table's "by" argument - We will compress the data down to each unique combination of these columns
    by_columns <- if (ignore_direction) {
        c("chr", "start")
    } else {
        c("chr", "start", "Type")
    }

    # This selection finds all CNVs that share a start pos (BUT not end pos - a separate segments search will account for these),
    # collapses CNVIDs and tumour types into lists, and filters out any matches that occur only within a tumour type
    shared_starts <-
        unique(coordinates[, .(commonCnvIds = paste0(sort(unique(CNV_ID_ext)), collapse=","),
                               tumourTypes = paste0(sort(unique(Category)), collapse=","),
                               nTumourTypes = uniqueN(Category),
                               Order = min(Order),
                               dft1Count = sum(Category=="DFT1"),
                               dft2Count = sum(Category=="DFT2"),
                               nonDftdCount = sum(startsWith(Category, "NonDFTD")),
                               nEnds = uniqueN(end)), by = by_columns],
               by = by_columns)[nEnds > 1]


    # This selection is as above, but for CNVs that match their end position only.
    by_columns[by_columns == "start"] <- "end"
    shared_ends <-
        unique(coordinates[, .(commonCnvIds = paste0(sort(unique(CNV_ID_ext)), collapse=","),
                               tumourTypes = paste0(sort(unique(Category)), collapse=","),
                               nTumourTypes = uniqueN(Category),
                               Order = min(Order),
                               dft1Count = sum(Category=="DFT1"),
                               dft2Count = sum(Category=="DFT2"),
                               nonDftdCount = sum(startsWith(Category, "NonDFTD")),
                               nStarts = uniqueN(start)), by = by_columns],
               by = by_columns)[nStarts > 1]

    # This selection finds CNVs that exactly match in both start and end positions, among different tumour types
    by_columns <- c(by_columns[1], "start", by_columns[2:length(by_columns)])
    shared_segments <-
        unique(coordinates[, .(commonCnvIds = paste0(sort(unique(CNV_ID_ext)), collapse=","),
                               tumourTypes = paste0(sort(unique(Category)), collapse=","),
                               nTumourTypes = uniqueN(Category),
                               dft1Count = sum(Category=="DFT1"),
                               dft2Count = sum(Category=="DFT2"),
                               nonDftdCount = sum(startsWith(Category, "NonDFTD")),
                               Order = min(Order),
                               nCNVs = uniqueN(CNV_ID_ext)),
                           by = by_columns][nCNVs > 1],
               by = by_columns)
    shared_segments[, nCNVs := NULL]

    # This selection matches the end pos of one interval to the start pos of another
    by_columns <- by_columns[by_columns != "end"]
    on_columns <- if (ignore_direction) {
        c("chr", start="end")
    } else {
        c("chr", "Type", start="end")
    }
    start_matches_end <- coordinates[coordinates, , on = on_columns, nomatch=0L]
    start_matches_end <- start_matches_end[, .(commonCnvIds = paste0(sort(unique(c(CNV_ID_ext, i.CNV_ID_ext))), collapse=","),
                                               tumourTypes = paste0(sort(unique(c(Category, i.Category))), collapse=","),
                                               nTumourTypes = uniqueN(c(Category, i.Category)),
                                               dft1Count = unique(rbindlist(list(data.table(CNV_ID_ext, Category), data.table(i.CNV_ID_ext, i.Category)), use.names=FALSE), by = "CNV_ID_ext")[, sum(Category=="DFT1")],
                                               dft2Count = unique(rbindlist(list(data.table(CNV_ID_ext, Category), data.table(i.CNV_ID_ext, i.Category)), use.names=FALSE), by = "CNV_ID_ext")[, sum(Category=="DFT2")],
                                               nonDftdCount = unique(rbindlist(list(data.table(CNV_ID_ext, Category), data.table(i.CNV_ID_ext, i.Category)), use.names=FALSE), by = "CNV_ID_ext")[, sum(startsWith(Category, "NonDFTD"))],
                                               Order = min(c(Order, i.Order))), by = by_columns]

    # Merge lists
    # First, some tidying up the columns
    shared_starts[, end := NA]
    shared_starts[, nEnds := NULL]

    shared_ends[, start := NA]
    shared_ends[, nStarts := NULL]

    start_matches_end[, end := start]

    shared_starts[, shared_breakpoint := "A_START_ONLY"]
    shared_ends[, shared_breakpoint := "B_END_ONLY"]
    shared_segments[, shared_breakpoint := "C_SEGMENT"]
    start_matches_end[, shared_breakpoint := "D_START_MATCHES_END"]

    segmatch <- rbindlist(list(shared_starts, shared_ends, shared_segments, start_matches_end), use.names = TRUE)
    if (ignore_direction) {
        setkeyv(segmatch, c("chr", "start", "end"))
    } else {
        setkeyv(segmatch, c("chr", "start", "end", "Type"))
    }
    setcolorder(segmatch)
    setorder(segmatch, shared_breakpoint, Order)
    return (segmatch)
}
