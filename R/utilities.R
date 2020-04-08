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

#' Generate a colour palette for the DFT1 clades
#' @export
clade_colours <- function() {
    rgb2hex <- function(rgb) {
        rgb2hex_inner <- function(rgb_) {
            digits <- c(0:9, LETTERS[1:6])
            leading <- floor(rgb_ / 16)
            trailing <- rgb_ %% 16
            paste0(digits[leading+1], digits[trailing+1])
        }
        paste0(c("#", sapply(rgb, rgb2hex_inner)), collapse="")
    }

    colours <- c(                      # Previous values (before March 6 2020)
        A1=rgb2hex(c(209, 38, 48)),    # 214, 43, 56
        A2=rgb2hex(c(240, 126, 167)),  # 241, 162, 177
        B=rgb2hex(c(103, 137, 206)),   # 92, 140, 196
        C=rgb2hex(c(73, 190, 114)),    # 16, 185, 127
        D=rgb2hex(c(248, 155, 15)),    # 254, 163, 56
        E=rgb2hex(c(38, 44, 113)))     # 37, 42, 113

    return (colours)
}

#' Write sequence information in a table to FASTA format
#' @export
write_fasta <- function(filename, table) {
    msg <- "Table should have labels in a 'label' column and sequences in a 'sequence' column"

    if(!("sequence" %in% colnames(table))) {
        stop(msg)
    }

    if (!("label" %in% colnames(table))) {
        stop(msg)
    }

    conn <- file(filename,
                 open = "w")

    for (i in 1:nrow(table)) {
        cat(sprintf(">%s\n%s\n\n", table[i, label], table[i, sequence]), file = conn)
    }

    close(conn)
}

#' Split a string into a vector of component characters
#' @export
split_chars <- function(s) {
    if (length(s) > 1) {
        stop("This function is for single strings only")
    }

    if (nchar(s) == 1) return (s)

    return(strsplit(s, "")[[1]])
}

#' Check if a vector of characters is a variant site
#' @export
is_variant_site <- function(chars, alphabet = c("DNA", "BINARY")) {
    dna_chars <- list(A=1,  C=2,  G=4,  T=8,
                      B=14, D=13, H=11, V=7,
                      R=5,  Y=10, S=6,  W=9,
                      K=12, M=3, N=15)

    binary_chars <- list(`0`=1, `1`=2, N=3)

    # Checks

    if (!is.character(chars)) chars <- as.character(chars)

    if (length(chars) == 1) {
        if (nchar(chars[1]) > 1) {
            chars <- split_chars(chars[1])
        }
    }

    alphabet <- match.arg(alphabet, c("DNA", "BINARY"))

    if (alphabet == "DNA") {
        translator <- dna_chars
    } else if (alphabet == "BINARY") {
        translator <- binary_chars
    }

    u <- sapply(toupper(chars), function(char) translator[[char]])

    # A variant site will reduce to zero via bitwise&, an invariant site will not
    # E.g. if A = 0001, C = 0010 and N = 1111:
    # AC (variant) -> 0001 & 0010 = 0000
    # AA (invariant) -> 0001 & 0001 = 0001
    # CN (invariant) -> 0010 & 1111 = 0010
    Reduce(bitwAnd, u) == 0
}
