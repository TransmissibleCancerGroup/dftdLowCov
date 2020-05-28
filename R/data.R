#' Loads the CNV table (supplementary table S3) in data.table format.
#' @importFrom "data.table" setDT
#' @importFrom "openxlsx" read.xlsx
#' @export
load_cnv_table <- function() {
    cnv_table_file <- system.file("extdata", "2020-05-07_TableS3_CNVs.xlsx", package = "dftdLowCov")
    cnv_table <- read.xlsx(cnv_table_file)
    setDT(cnv_table)
    cnv_table
}

#' Loads the tumour sample table (supplementary table S1A) in data.table format.
#' @importFrom "data.table" setDT
#' @importFrom "openxlsx" read.xlsx
#' @export
load_sample_table <- function() {
    table_file <- system.file("extdata", "2020-05-07_TableS1_Tumours-Hosts.xlsx", package = "dftdLowCov")
    sample_table <- read.xlsx(table_file, sheet = 1)
    setDT(sample_table)
    sample_table
}

#' load_chromosome lengths into a list
#' @importFrom "data.table" rbindlist
#' @export
load_chromosome_lengths <- function() {
    filename <- system.file("extdata", "devil_chrom_lengths.tsv", package = "dftdLowCov")
    chr.lengths <- fread(filename)
    chrlengths_truncated_x <- rbindlist(list(chr.lengths[CHROM != "Chrx"], chr.lengths[CHROM == "Chrx", .(CHROM, LENGTH = 54800000)]))
    full_chrlengths <- copy(chr.lengths)
    chrlengths_truncated_x[, seqnames := toupper(sub("Chr", "", CHROM))]
    full_chrlengths[, seqnames := toupper(sub("Chr", "", CHROM))]
    list(full = full_chrlengths,
         truncated = chrlengths_truncated_x)
}
