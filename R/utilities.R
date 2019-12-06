#' Calculate number of bases covered by a table of CNVs
#' Start and end coordinates must be in the format that the end coordinate
#' is the position of the first base NOT covered by the CNV (SAM format)
#' @importFrom "data.table" data.table setkey setkeyv foverlaps rbindlist copy key
#' @export
get_depth <- function(data_table, chrom_lengths = NULL) {
    get_depth_for_single_chr <- function(dt) {
        .disjoin <- function(.dt) {
            starts <- unique(.dt$start)
            ends <- unique(.dt$end)
            adjstart <- head(sort(unique(c(starts, ends))), -1)
            adjend <- tail(sort(unique(c(ends, starts))), -1)
            adj <- data.table(start=adjstart, end=adjend, width=adjend - adjstart)
            setkey(adj, start, end)
            setkeyv(.dt, key(adj))
            ol <- foverlaps(.dt, adj, nomatch=0L, minoverlap=1L)
            res <- unique(ol, by = c("start", "end"))
            res[, intervalID := .I]
            res[, .(chr, start, end, intervalID)]
        }

        intervals <- .disjoin(dt)

        # Adjust end points of all intervals so that exactly matching end+start coords don't count as a match
        orig <- copy(dt)
        orig[, end := end - 1]
        intervals[, end := end - 1]

        # Set keys for foverlaps
        setkey(intervals, chr, start, end)
        setkey(orig, chr, start, end)

        # Search the CNVs in "orig" using the disjoined intervals in "intervals" (include any missing intervals as NA)
        ol <- foverlaps(intervals, orig, nomatch=NA)

        # Annotate the depth
        ol[, depth := as.integer(0)]
        result <- ol[ol[!is.na(start), .N, by = intervalID], depth := as.integer(N), on = "intervalID"]
        result <- unique(result[, .(chr, start=i.start, end=i.end, depth)])
        result[, width := as.integer(end - start + 1)]

        return(result[1:.N])
    }

    results <- list()
    for (CHR in c(1:6, "X")) {
        result <- get_depth_for_single_chr(data_table[chr==CHR])

        # Adjust 0-coverage figure for any positions at the end of the chromosome
        # that are not entered in the table
        if (!is.null(chrom_lengths)) {
            chrom_size <- chrom_lengths[chr == CHR, as.integer(LENGTH)]
            if (length(chrom_size) > 0 ) {
                max_end <- result[, max(end)]
                if (chrom_size > max_end) {
                    new_row <- copy(result[1])
                    new_row[, c("start", "end", "depth") := .(max_end+1, chrom_size, 0)]
                    new_row[, width := end - start + 1]
                    result <- rbindlist(list(result, new_row))
                }
            }
        }
        results[[length(results) + 1]] <- result
    }

    rbindlist(results)
}

#' Returns TRUE if vectors a and b overlap
#' @export
is.overlap <- function(a, b) {
    stopifnot(length(a) == 2 & length(b) == 2)
    return (a[2] >= b[1] & a[1] < b[2])
}
