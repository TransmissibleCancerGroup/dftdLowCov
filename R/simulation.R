make_bins <- function(chrom_lengths, binsize = 100000) {
    bins <-
        chrom_lengths[, .(seqnames = toupper(sub("Chr", "", CHROM)),
                          start = seq(0, max(LENGTH) - binsize, binsize),
                          end = seq(0, max(LENGTH) - binsize, binsize) + binsize - 1),
                      by = CHROM][, .(seqnames, start, end)]


    setkey(bins, seqnames, start, end)
    # setkeyv(subject, key(bins))
    bins
}

#' Calculates the depth of coverage (of CNV segments) for the genome,
#' discretised into bins of size `binsize`
depth_of_coverage <- function(data, chrom_lengths, binsize = 100000) {
    bins <- make_bins(chrom_lengths, binsize)
    setkey(data, seqnames, start, end)
    ol <- foverlaps(bins, data)
    ol[, present := ifelse(is.na(start), 0L, 1L)]
    ol[, depth := sum(present), by = .(seqnames, i.start, i.end)]
    ol
}

doc_summary <- function(doc_result) {
    doc_result[, .(nBins = uniqueN(i.start)), by = .(seqnames, depth)][order(seqnames, depth)]
}

randomize_positions_within_chrom <- function(data, chromosome_ranges, wraparound = FALSE) {
    setkey(chromosome_ranges, seqnames)
    data <- chromosome_ranges[data, .(seqnames, chromsize = width, start = i.start, end = i.end, width = i.width, State, Unique.ID)]
    data[, rstart := as.integer(0)]

    if (wraparound) {
        data[, rstart := sample.int(chromsize, size = .N, replace = TRUE), by = seqnames]
    } else {
        data[chromsize - width + 1 > 0, rstart := sapply(chromsize - width + 1, sample.int, size = 1)]
    }

    data[, rend := as.integer(rstart + width - 1)]
    data[, wrapped := FALSE]
    data <- data[, .(seqnames, chromsize, start = rstart, end = rend, width = as.integer(rend - rstart + 1), Unique.ID, State, wrapped)]
    data[, c("start", "end", "width") := lapply(.SD, as.numeric), .SDcols = c("start", "end", "width")]

    if (wraparound) {
        data <- wrap_around(data, data[["chromsize"]])
    } else {
        data[end > chromsize, end := chromsize]
    }

    data[, chromsize := NULL]
    data
}

#' Takes a data frame of segments and randomises their start
#' positions (up to max start), while keeping their width constant.
#' If the end position exceeds maxstart, the segment is 'wrapped around'
#' by cutting into two segments, the second of which is mapped to position 0.
randomize_positions_over_whole_genome <- function(segments, limit) {
    segments <- copy(segments)
    setDT(segments)

    # Check input:
    # Segments must have a 'start' column and at least one of 'end' and 'width'
    if (!("start" %in% colnames(segments))) {
        stop("No 'start' column")
    }

    if (!("end" %in% colnames(segments))) {
        # need "width"
        if (!("width" %in% colnames(segments))) {
            stop("Neither 'end' nor 'width' in columns")
        } else {
            segments[, end := start + width - 1]
        }
    } else {
        if (!("width" %in% colnames(segments))) {
            if (!("end" %in% segments)) {
                stop("Neither 'end' nor 'width' in columns")
            } else {
                segments[, width := end - start + 1]
            }
        }
    } # input OK


    new_starts <- sample.int(limit, segments[, .N], replace = TRUE)
    new_ends <- new_starts + segments[["width"]] - 1

    new_segments <- data.table(start = new_starts,
                               end = new_ends,
                               width = segments[["width"]],
                               Unique.ID = segments[["Unique.ID"]],
                               State = segments[["State"]])

    wrap_around(new_segments, limit)
}

#' Wrap segments that exceed the limit value back into the start of the sequence
wrap_around <- function(dt, limit) {
    dt[, .limit := as.numeric(limit)]
    dt[, c("start", "end", "width") := lapply(.SD, as.numeric), .SDcols = c("start", "end", "width")]
    overspill <- copy(dt[end > .limit])
    dt[end > .limit, c("end", "width") := .(.limit, .limit - start + 1L)]
    overspill[, c("start", "end", "width") := .(0L, end %% .limit - 1L, end %% .limit)]
    dt <- rbindlist(list(dt, overspill))
    dt[, wrapped := FALSE]
    dt[Unique.ID %in% overspill$Unique.ID, wrapped := TRUE]
    dt[(wrapped), Unique.ID := paste0(Unique.ID, ".", letters[1:.N]), by = Unique.ID]
    dt[, .limit := NULL]
    dt
}

assign_seqnames <- function(dt, chromosome_ranges) {
    setkey(chromosome_ranges, start, end)
    ol <- foverlaps(dt, chromosome_ranges)
    ol[, N:=.N, by = Unique.ID]
    ol[N > 1, wrapped := TRUE]
    ol[, adj.start := pmax(start, i.start) - start]
    ol[, adj.end := pmin(end, i.end) - start]
    ol[, adj.width := adj.end - adj.start + 1]
    ol[, .(CHROM, seqnames, start = adj.start, end = adj.end, width = adj.width, Unique.ID, State, wrapped)]
}

randomize_positions <- function(dataframe, chrom_lengths, keep_to_original_chrom = FALSE, wraparound = TRUE) {
    chrom_ranges <- chrom_ranges_from_lengths(chrom_lengths)
    if (keep_to_original_chrom) {
        return (randomize_positions_within_chrom(dataframe, chrom_ranges, wraparound))
    } else {
        return (assign_seqnames(randomize_positions_over_whole_genome(dataframe, limit = chrom_ranges[, max(end)]),
                                chrom_ranges))
    }
}

#' Calculate chromosome ranges from chromosome lengths
chrom_ranges_from_lengths <- function(chr_lengths) {
    chr_lengths <- copy(chr_lengths)
    chr_lengths[, LENGTH.CS := cumsum(as.numeric(LENGTH))]
    chrom_ranges <- chr_lengths[, .(CHROM,
                                    seqnames = toupper(sub("Chr", "", CHROM)),
                                    start=c(0, LENGTH.CS[1:(.N-1)]),
                                    end=LENGTH.CS-1)]
    chrom_ranges[, width := end - start + 1]
    return(chrom_ranges)
}

#' Data cleaning ready for making simulations
#' @export
preprocess_for_simulation <- function(cnvtable, chrom_lengths_file) {
    #########################
    # LOAD CHROMOSOME LENGTHS
    #########################
    lengths <- load_chromosome_lengths()
    chrlengths_truncated_x <- lengths$truncated
    full_chrlengths <- lengths$full

    ##################################################
    # Process input CNVs to fit the chromosome lengths
    ##################################################
    # Data should have inclusive end points (1) i.e. the quoted end position is covered by the CNV,
    # so an interval of 10000 20000 covers bases 10000, 10001, 10002, ..., 19999, 20000.
    # As loaded, the data has exclusive end points (2), so for the above example, position 20000 is NOT
    # covered by the CNV.
    # Here I make the adjustment from type (2) to type (1). Why? This frees me from having to remember
    # which system is being used at which particular time, and makes it less likely that I'll make
    # quite so many off by one errors.

    # Make end points inclusive:
    cnvtable[, end := end - 1]

    # Remove any segments that start after the end of the chromosome
    cnvtable[chrlengths_truncated_x, chrom_end := LENGTH, on="seqnames"]
    cnvtable <- cnvtable[start <= chrom_end]

    # Adjust any end points that over shoot the end of the chromosome
    cnvtable[end > chrom_end, end := chrom_end]
    # I assume any end points that fall within a few bases of the end of the chromosome
    # are due to preexisting errors in the data, so I set them to the end of the chromosome
    cnvtable[chrom_end - end < 5, end := chrom_end]

    # Readjust widths to match new end points
    cnvtable[, width := end - start + 1]

    return (cnvtable)
}

#' Counts the number of bases covered by 0, 1, 2 etc CNV intervals
#' @param data - data.table of CNV intervals
#' @param chrlengths - data.table of chromosome lengths
#' @param endpos_is_inclusive - TRUE means an interval [start=10, end=20] includes the last position,
#' so is 11 bases long (i.e. 10, 11, 12, ..., 19, 20]. FALSE means an interval [start=10, end=20]
#' excludes the last position, and is 10 bases long, because it stops at 19.
#' @importFrom "IRanges" IRanges disjoin
#' @export
coverage_depth_summary <- function(data, chrlengths, endpos_is_inclusive) {
    genome_counts <- data.table(seqnames = character(0), depth = integer(0), nbases = integer(0))

    # Per chromosome
    for (seqname in unique(data$seqnames)) {
        total_chr_bases <- chrlengths[CHROM == paste0("Chr", tolower(seqname)), LENGTH]

        chrdata <- data[seqnames==seqname & start < total_chr_bases]

        end_adjustment <- ifelse(endpos_is_inclusive, 0, 1)
        disjoint <- as.data.table(disjoin(IRanges(start = chrdata[, start],
                                                           end = chrdata[, end-end_adjustment])))
        setkey(disjoint, start, end)
        ol <- foverlaps(chrdata, disjoint, by.x = c("start", "end"), nomatch=0L)
        stats <- ol[, .(width = end - start + 1, depth = uniqueN(Unique.ID)), by = .(start, end)]
        counts <- stats[, .(nbases = sum(width)), by = depth][order(depth)]

        # zero coverage
        total_chr_bases <- max(total_chr_bases, chrdata[, max(end)]) - ifelse(endpos_is_inclusive, 0, 1)
        total_covered_bases <- counts[, sum(nbases)]
        counts <- rbindlist(list(counts, data.table(depth=0, nbases=total_chr_bases - total_covered_bases)))

        counts[, seqnames := seqname]
        genome_counts <- rbindlist(list(genome_counts, counts), use.names = TRUE)
    }

    genome_counts[order(seqnames, depth)]
}


#' Simulation code
#' @param input_data data.table; table of CNVs
#' @param chrlengths data.table; table of the Devil chromosome lengths
#' @param simulation_name string; used in name of output files
#' @param nsims int; number of simulation replicates
#' @param base_outdir path; output will be written under this path
#' @param plot_every int; draw a chromosome map every `plot_every` reps
#' @param write_output bool; write output to disk?
#' @param keep_to_original_chrom bool; should shuffled CNV segments be restricted
#'     to their original chromosome, or be allowed to move genome wide?
#' @param hard_bounds bool; should chromosome ends act as hard boundaries for shuffled
#'     CNVs? Or should CNVs be allowed to straddle over chromosome ends?
#' @param parametric bool; if TRUE, then CNV widths are drawn from a Gamma distribution,
#'     fitted to the observed widths in the input data. If FALSE, the original widths
#'     are used.
#' @importFrom "logging" logwarn
#' @importFrom "ggplot2" ggplot aes theme_classic scale_fill_manual ylab xlab ggtitle geom_point
#' @export
do_simulation <-
    function(input_data, chrlengths, simulation_name, nsims, base_outdir,
             plot_every = 100,
             write_output = FALSE,
             keep_to_original_chrom = FALSE,
             hard_bounds = FALSE,
             parametric = FALSE) {

        # NB this function is huge, but it works, and as such I am reluctant to refactor
        # it - KG April 2020

        if (!keep_to_original_chrom & hard_bounds) {
            stop("Incompatible arguments: keep_to_original_chrom=FALSE, hard_bounds=TRUE")
        }

        GAIN.COLOUR <- plot_chromosome_map(get_plot_options = TRUE)[["GAINS.SEGMENT.COLOUR"]]
        LOSS.COLOUR <- plot_chromosome_map(get_plot_options = TRUE)[["LOSSES.SEGMENT.COLOUR"]]

        # Configure plotting options
        opts = list(CHROM.FILL.COLOUR = "grey70",
                    GAINS.SEGMENT.OUTLINE.COLOUR = "#EA6A58",
                    LOSSES.SEGMENT.OUTLINE.COLOUR = "#94A6F1",
                    CHROM.STRIPE.COLOUR = NA,
                    GAINS.HIST.COLOUR = NA,
                    LOSSES.HIST.COLOUR = NA,
                    CHROM.WIDTH = (4/3),
                    CHROM.BUFFER = (2/3),
                    MIN.X = -(80/3),
                    MAX.X = (80/3))

        SIMULATION.PLOT.EVERY <- plot_every

        do_sims <- write_output # don't accidentally overwrite simulations
        if (!do_sims) {
            logwarn("do_sims is set to FALSE, so no simulations will be performed")
        }

        OUTDIR <- file.path(base_outdir, "simulation_results", simulation_name)
        if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)
        NAME_STEM <- "samples"

        # Process a local copy of the input data
        data <- copy(input_data[!is.na(Unique.ID) & seqnames != "X-loss"])
        data <- data[!(seqnames == "X" & end > 54800000)]
        # data[, end := end - 1] # make so endpoint is inclusive
        data[, width := end - start + 1]

        # Save a copy of the data to the simulation directory
        fwrite(input_data, file = file.path(OUTDIR, "input_data.tsv"), sep = '\t')

        if (parametric) {
            generator <- fitdistrplus::fitdist(data$width, "gamma", method = "mme")
        }

        pdf(file.path(OUTDIR, paste0(NAME_STEM, ".pdf")), width = 20, height = 10)
        layout(matrix(c(1,2,3,4,
                        1,2,3,4,
                        1,2,3,4,
                        1,2,3,4,
                        5,6,7,8,
                        5,6,7,9), ncol = 4, byrow=T))
        par(mar = c(2.5,1,0.5,1), oma = c(1, 0, 4, 0))

        samples <- list()
        summaries <- list()
        pval_summaries <- list()

        # Relevant to p-val calculation
        gain_pvals <- rep(0, 100)
        input_coverage_summary_gain <- coverage_depth_summary(data[State == "gain"], chrlengths, endpos_is_inclusive=TRUE)[, .(nbases = sum(nbases)), by = depth][order(depth)]
        gain_pvals_max_entry <- input_coverage_summary_gain[, max(depth)]

        loss_pvals <- rep(0, 100)
        input_coverage_summary_loss <- coverage_depth_summary(data[State == "loss"], chrlengths, endpos_is_inclusive=TRUE)[, .(nbases = sum(nbases)), by = depth][order(depth)]
        loss_pvals_max_entry <- input_coverage_summary_loss[, max(depth)]

        input_coverage_summary_gain[, State := "gain"]
        input_coverage_summary_gain[, rep := 0]
        input_coverage_summary_loss[, State := "loss"]
        input_coverage_summary_loss[, rep := 0]
        input_coverage_summary <- rbindlist(list(input_coverage_summary_gain, input_coverage_summary_loss))
        setorder(input_coverage_summary, -depth)
        input_coverage_summary[, nbases_at_least := cumsum(nbases), by = State]
        setorder(input_coverage_summary, State, depth)
        fwrite(input_coverage_summary, file = file.path(OUTDIR, paste0("input", ".pval_summaries.tsv")), sep = '\t')

        # Plot input chromosome map
        for (chr in c(1:6, "X")) {
            dt.gains <- dftdLowCov::rectangle_packing4(data[State == "gain" & seqnames == chr])
            dt.losses <- dftdLowCov::rectangle_packing4(data[State == "loss" & seqnames == chr])
            plot_chromosome_map(dt.gains, dt.losses, chrlengths, plot_options = opts,
                                chr_override = chr) # was full_chrlengths, but now we want to avoid plotting the end of ChrX, where we never called CNVs, and is excluded from the randomisation
        }
        mtext("Input data", outer = TRUE)

        gain_doc <- depth_of_coverage(data[State == "gain"], chrlengths)
        loss_doc <- depth_of_coverage(data[State == "loss"], chrlengths)
        gain_summ <- doc_summary(gain_doc)
        loss_summ <- doc_summary(loss_doc)

        gain_summ[, rep := 0]
        gain_summ[, State := "gain"]
        loss_summ[, rep := 0]
        loss_summ[, State := "loss"]

        summary <- rbindlist(list(gain_summ, loss_summ))

        plot(nbases_at_least ~ depth, data = input_coverage_summary[State == "gain"], col = GAIN.COLOUR, main = "Gains genome wide summary", bty = "n", xlab = "Depth", ylab = "N", pch = 20, cex = 2)
        plot(nbases_at_least ~ depth, data = input_coverage_summary[State == "loss"], col = LOSS.COLOUR, main = "Losses genome wide summary", bty = "n", xlab = "Depth", ylab = "N", pch = 20, cex = 2)

        if (!do_sims) {
            dev.off()
            return()
        }

        for (i in 1:nsims) {
            if (parametric) {
                data[, width := as.integer(ceiling(rgamma(.N, shape = generator$estimate[["shape"]],
                                                          rate = generator$estimate[["rate"]])))]
            }
            rdata <- randomize_positions(copy(data),
                                         copy(chrlengths),
                                         keep_to_original_chrom = keep_to_original_chrom,
                                         wraparound = !hard_bounds)

            rdata[, rep := i]
            samples[[i]] <- rdata

            gain_doc <- depth_of_coverage(rdata[State == "gain"], chrlengths)
            loss_doc <- depth_of_coverage(rdata[State == "loss"], chrlengths)
            gain_summ <- doc_summary(gain_doc)
            loss_summ <- doc_summary(loss_doc)

            gain_summ[, rep := i]
            gain_summ[, State := "gain"]
            loss_summ[, rep := i]
            loss_summ[, State := "loss"]

            summary <- rbindlist(list(gain_summ, loss_summ))
            summaries[[i]] <- summary

            # Empirical p-value
            rep_coverage_summary_gain <- coverage_depth_summary(rdata[State=="gain"], chrlengths, endpos_is_inclusive = TRUE)[, .(nbases = sum(nbases)), by = depth][order(depth)]
            if (rep_coverage_summary_gain[, max(depth)] > gain_pvals_max_entry) {
                gain_pvals_max_entry <- rep_coverage_summary_gain[, max(depth)]
            }

            # Grow the array if necessary
            if (gain_pvals_max_entry > length(gain_pvals)) {
                gain_pvals <- c(gain_pvals, rep(0, 100))
            }

            for (j in 1:gain_pvals_max_entry) {
                input_value <- input_coverage_summary_gain[depth >= j, sum(nbases)]
                sim_value <- rep_coverage_summary_gain[depth >= j, sum(nbases)]
                if (sim_value >= input_value) {
                    gain_pvals[j] <- gain_pvals[j] + 1
                }
            }

            rep_coverage_summary_loss <- coverage_depth_summary(rdata[State=="loss"], chrlengths, endpos_is_inclusive = TRUE)[, .(nbases = sum(nbases)), by = depth][order(depth)]
            if (rep_coverage_summary_loss[, max(depth)] > loss_pvals_max_entry) {
                loss_pvals_max_entry <- rep_coverage_summary_loss[, max(depth)]
            }

            # Grow the array if necessary
            if (loss_pvals_max_entry > length(loss_pvals)) {
                loss_pvals <- c(loss_pvals, rep(0, 100))
            }

            for (k in 1:loss_pvals_max_entry) {
                input_value <- input_coverage_summary_loss[depth >= k, sum(nbases)]
                sim_value <- rep_coverage_summary_loss[depth >= k, sum(nbases)]
                if (sim_value >= input_value) {
                    loss_pvals[k] <- loss_pvals[k] + 1
                }
            }
            rep_coverage_summary_gain[, State := "gain"]
            rep_coverage_summary_gain[, rep := i]
            rep_coverage_summary_loss[, State := "loss"]
            rep_coverage_summary_loss[, rep := i]
            rep_coverage_summary <- rbindlist(list(rep_coverage_summary_gain, rep_coverage_summary_loss))
            setorder(rep_coverage_summary, -depth)
            rep_coverage_summary[, nbases_at_least := cumsum(nbases), by = State]
            setorder(rep_coverage_summary, State, depth)
            pval_summaries[[i]] <- rep_coverage_summary

            if (i %% SIMULATION.PLOT.EVERY == 0) {
                # plot
                for (chr in c(1:6, "X")) {
                    dt.gains <- dftdLowCov::rectangle_packing4(rdata[State == "gain" & seqnames == chr])
                    dt.losses <- dftdLowCov::rectangle_packing4(rdata[State == "loss" & seqnames == chr])
                    plot_chromosome_map(dt.gains, dt.losses, chrlengths, plot_options = opts,
                                        chr_override = chr) # was full_chrlengths, but now we want to avoid plotting the end of ChrX, where we never called CNVs, and is excluded from the randomisation
                }
                mtext(paste0("Simulation rep ", i), outer = TRUE)

                plot(nbases_at_least ~ depth, data = rep_coverage_summary[State == "gain"], col = GAIN.COLOUR, main = "Gains genome wide summary", bty = "n", xlab = "Depth", ylab = "N", pch = 20, cex = 2)
                plot(nbases_at_least ~ depth, data = rep_coverage_summary[State == "loss"], col = LOSS.COLOUR, main = "Losses genome wide summary", bty = "n", xlab = "Depth", ylab = "N", pch = 20, cex = 2)
                # boxplot(nBins ~ depth, data = summary[State == "gain"], col = GAIN.COLOUR, main = "Gains genome wide summary", bty = "n", xlab = "Depth", ylab = "N")
                # boxplot(nBins ~ depth, data = summary[State == "loss"], col = LOSS.COLOUR, main = "Losses genome wide summary", bty = "n", xlab = "Depth", ylab = "N")
            }
            print(sprintf("Simulation %d/%d", i, nsims))
        }
        dev.off()
        summaries.dt <- rbindlist(summaries); rm(summaries)
        samples.dt <- rbindlist(samples); rm(samples)
        pval_summaries.dt <- rbindlist(pval_summaries); rm(pval_summaries)

        fwrite(summaries.dt, file.path(OUTDIR, paste0(NAME_STEM, ".summaries.tsv")), sep = '\t')
        fwrite(samples.dt, file.path(OUTDIR, paste0(NAME_STEM, ".tsv")), sep = '\t')
        fwrite(pval_summaries.dt, file.path(OUTDIR, paste0(NAME_STEM, ".pval_summaries.tsv")), sep = '\t')

        plots <- list()
        for (chr in as.character(c(1:6, "X"))) {
            plots[[chr]] <-
                ggplot(data = summaries.dt[seqnames == chr, .(nBins = sum(nBins)), by = .(depth, rep, State)],
                       aes(x = as.factor(depth), y = nBins, fill = State, colour = State)) +
                #geom_violin(position = position_dodge(.75), trim = TRUE, scale = "width") +
                theme_classic() +
                scale_fill_manual(values = c(GAIN.COLOUR, LOSS.COLOUR)) +
                ylab("N") + xlab("Depth of Coverage") + ggtitle(paste0("Chr", chr)) +
                geom_point(position = position_dodge(0.5))
            #geom_boxplot(position = position_dodge(0.5), outlier.shape = NA)
        }
        plots[["all"]] <-
            ggplot(data = summaries.dt[, .(nBins = sum(nBins)), by = .(depth, rep, State)],
                   aes(x = as.factor(depth), y = nBins, fill = State, colour = State)) +
            #geom_violin(position = position_dodge(.75), trim = TRUE, scale = "width") +
            theme_classic() +
            scale_fill_manual(values = c(GAIN.COLOUR, LOSS.COLOUR)) +
            ylab("N") + xlab("Depth of Coverage") + ggtitle("Genome-wide") +
            geom_point(position = position_dodge(0.5))
        #geom_boxplot(position = position_dodge(0.5), outlier.shape = NA) + theme(legend.justification = "left")

        legend <- get_legend(plots[["all"]])

        MARGIN <- c(1,1,1,1)
        a <- align_plots(
            plots[["1"]] + theme(legend.position = "none", plot.margin = unit(MARGIN, "lines")),
            plots[["2"]] + theme(legend.position = "none", plot.margin = unit(MARGIN, "lines")),
            plots[["3"]] + theme(legend.position = "none", plot.margin = unit(MARGIN, "lines")),
            plots[["4"]] + theme(legend.position = "none", plot.margin = unit(MARGIN, "lines")),
            plots[["5"]] + theme(legend.position = "none", plot.margin = unit(MARGIN, "lines")),
            plots[["6"]] + theme(legend.position = "none", plot.margin = unit(MARGIN, "lines")),
            plots[["X"]] + theme(legend.position = "none", plot.margin = unit(MARGIN, "lines")),
            plots[["all"]] + theme(legend.position = "none", plot.margin = unit(MARGIN, "lines")),
            align = 'v', axis = 'l')

        top_row = plot_grid(a[[1]], a[[2]], a[[3]], a[[4]], NULL, rel_widths = c(1,1,1,1,0.2),
                            nrow = 1)
        bottom_row = plot_grid(a[[5]], a[[6]], a[[7]], a[[8]], legend, rel_widths = c(1,1,1,1,0.2),
                               nrow = 1)
        final_plot <- plot_grid(top_row, bottom_row, ncol = 1, scale = 1.0)
        save_plot(filename = file.path(OUTDIR, paste0(NAME_STEM, "_coverage_summary.pdf")),
                  plot = final_plot,
                  base_height = 10, base_width = 20)

        # Empirical pvals
        # gain_pvals <- gain_pvals[gain_pvals_max_entry] / nsims
        # loss_pvals <- loss_pvals[loss_pvals_max_entry] / nsims

        return (list(gainp = gain_pvals, lossp = loss_pvals))
    }
