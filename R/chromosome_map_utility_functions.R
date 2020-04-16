# Helper functions used in drawing chromosome maps, and in making simulations

#' @export
load_subject_from_finalised_cnvs <- function(filename) {
    cnvs <- as.data.table(readxl::read_xlsx(filename))
    subject <- unique(cnvs[, .(seqnames=chr, start, end, width, State=Type, Unique.ID=New.CNV_ID_ext, `Cell-line`, `DFT2`, `Non-DFTD`)], by = "Unique.ID")
    return (subject)
}

