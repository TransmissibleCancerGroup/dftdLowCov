# Functions to decide if a CNV belongs to a binary set, e.g.:
#   1) is Cell-line specific (isCellLine())
#   2) is Carcinoma specific (isCarcinoma())
#   3) is part of a multi-copy gain or loss (isMultiGain() / isMultiLoss())
#   4) is associated with marker 5 (isMarker5())
# All functions take the cnv table as a data.table, and return a vector of
# CNV IDs for the CNVs matching the condition

#' @export
isCellLine <- function(cnvtable) {
    cellline_cnvs <- cnvtable[!is.na(`Cell-line`) & `Cell-line` == "Cell-line", New.CNV_ID_ext]
    return (cellline_cnvs)
}

#' @export
isCarcinoma <- function(cnvtable) {
    carcinoma_unique_cnvs <- (function(v) v[!is.na(v)])(
        sapply(1:nrow(cnvtable), function(i) {
            samples <- cnvtable[i, dftdLowCov::parse_sample_cnv_group(Sample_CNV_Group)[, sample]]
            if(length(setdiff(samples, c("340T", "997T1"))) == 0) {
                return(cnvtable[i, New.CNV_ID_ext])
            }
            return (NA_character_)
        })
    )
    return (carcinoma_unique_cnvs)
}

#' @export
isMultiGain <- function(cnvtable) {
    # Remove multiple copy-number gains and losses
    # 1) label rows where the New.CNV_ID_ext does have 'G'
    # 2) label rows where the sum of this quantity per New.CNV_ID is > 0
    # 3) Mark for removal rows that have G, as long as they belong to a group that also has not G :) so simple and not annoying at all :)
    cnvtable[, has_G_in_cnvext := grepl("G", New.CNV_ID_ext)]
    cnvtable[, G_and_not_G := any(has_G_in_cnvext) & any(!has_G_in_cnvext), by = New.CNV_ID]
    multigain_cnvs <- cnvtable[(has_G_in_cnvext) & (G_and_not_G), New.CNV_ID_ext]

    # 4) For the groups that only have G, only keep the first guy
    cnvtable[(has_G_in_cnvext) & !(G_and_not_G), pos_in_group := seq_len(.N), by = New.CNV_ID]
    multigain_cnvs <- c(multigain_cnvs, cnvtable[(has_G_in_cnvext) & !(G_and_not_G) & pos_in_group > 1, New.CNV_ID_ext])

    # Kill the labels
    cnvtable[, c("has_G_in_cnvext", "G_and_not_G", "pos_in_group") := NULL]

    return (multigain_cnvs)
}

#' @export
isMultiLoss <- function(cnvtable) {
    cnvtable[, has_L_in_cnvext := grepl("L", New.CNV_ID_ext)]
    cnvtable[, L_and_not_L := any(has_L_in_cnvext) & any(!has_L_in_cnvext), by = New.CNV_ID]
    multiloss_cnvs <- cnvtable[(has_L_in_cnvext) & (L_and_not_L), New.CNV_ID_ext]
    cnvtable[(has_L_in_cnvext) & !(L_and_not_L), pos_in_group := seq_len(.N), by = New.CNV_ID]
    multiloss_cnvs <- c(multiloss_cnvs, cnvtable[(has_L_in_cnvext) & !(L_and_not_L) & pos_in_group > 1, New.CNV_ID_ext])

    # Remove labels
    cnvtable[, c("has_L_in_cnvext", "L_and_not_L", "pos_in_group") := NULL]

    return (multiloss_cnvs)
}

#' @export
isMarker5 <- function(cnvtable) {
    # Remove anything to do with "Marker 5"
    m5_associated_cnvids <- c(6, 11, 16, 17, 18, 285, 1803)
    m5_cnvs <- cnvtable[New.CNV_ID %in% m5_associated_cnvids, New.CNV_ID_ext]

    return (m5_cnvs)
}

#' @export
isDFT1 <- function(cnvtable) {
    cnvtable[is.na(DFT2) & is.na(`Non-DFTD`), New.CNV_ID_ext]
}

#' @export
isDFT2 <- function(cnvtable) {
    cnvtable[!is.na(DFT2), New.CNV_ID_ext]
}

#' @export
isNonDFTD <- function(cnvtable) {
    cnvtable[!is.na(`Non-DFTD`), New.CNV_ID_ext]
}
