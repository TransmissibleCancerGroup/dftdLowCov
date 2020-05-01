require(data.table)

#' Calculates the depth at
coverage_histogram <- function(dt, binsize = 100000) {
    bins <- data.table(start = seq(0, dt[, max(end)] - binsize, binsize))
    bins[, end := start + binsize - 1]
    setkey(bins, start, end)
    setkey(dt, NULL)
    foverlaps(dt, bins, nomatch = 0L)[, .N, by = .(start, end)]
}

#' A faster algorithm to pack aligned reads into layers. Doesn't need IRanges.
rectangle_packing4 <- function(dt, sortby = "width") {
    if (nrow(dt) == 0) {
        dt[, layer := 0]
        return (dt)
    }

    # Function to assign the lowest available layer, avoiding
    # any of the layers in `disallowed`
    assignLayer <- function(disallowed) {
        if (length(disallowed) == 0 || (length(disallowed) == 1 && disallowed == 0)) return (1)
        M <- max(disallowed)
        candidate <- setdiff(1:M, disallowed)
        ifelse(length(candidate) == 0, M+1, candidate)
    }

    dt[, origOrder := .I]
    setkey(dt, start, end)
    hits <- foverlaps(dt, dt)[Unique.ID != i.Unique.ID][, .(queryHits = Unique.ID, subjectHits = i.Unique.ID)]

    if (sortby == "width") {
        setorder(dt, -width)
    } else {
        setorder(dt, start)
    }
    dt[, layer := as.integer(0)]

    # No conflicts possible for first segment, so place in bottom layer
    dt[1, layer := 1]

    # processed <- 1
    N <- dt[, .N]
    if (N > 1) {
        for (currID in dt[2:.N, Unique.ID]) {
            # Search for any segments that overlap the current segment. If any layers are assigned
            # for the hits, then they are disallowed for our current segment.
            disallowed_layers <- dt[Unique.ID %in% hits[queryHits == currID, subjectHits], unique(layer)]
            assigned_layer <- assignLayer(disallowed_layers)
            dt[Unique.ID == currID, layer := assigned_layer]

            # Just for counting and showing progress
            # processed <- processed + 1
            # print(paste("(", processed, "/", N, ")"))
        }
    }

    # Put the dt back how we found it
    setorder(dt, origOrder)
    dt[, origOrder := NULL]
    dt
}

#' Plot an individual chromosome's map. Assumes dt.gains and dt.losses
#' have already been produced, including setting the plotting layers
#' for each segment with rectangle packing.
#' Additional segment colouring options can be set using `category_colours`.
#' @param dt_gains data.table; table of CNV gain segments on the chromosome. Must have
#'     the plotting layer annotation present (as added by rectangle_packing4).
#' @param dt_losses data.table; table of CNV loss segments on the chromosome. Must have
#'     layer annotation.
#' @param chr_lengths data.table; table of chromosome lengths.
#' @param category_colours list; mapping of Category values to colours. These will be used
#'     to override the GAINS.SEGMENT.COLOUR and LOSS.SEGMENT.COLOUR. Only applies if the
#'     input data has a Category column.
#' @param plot_options list; user can provide plotting options via this list.
#' @param get_plot_options list; returns the default plot options.
#' @param chr_override character; set the chromosome that the segments should be plotted
#'     onto. Used in the case that the chromosome can't be inferred from the data.
#' @export
plot_chromosome_map <- function(dt_gains,
                                dt_losses,
                                chr_lengths,
                                category_colours = NULL,
                                plot_options = NULL,
                                get_plot_options = FALSE,
                                chr_override = NULL) {

    default.plot.options <- list(
        GAINS.SEGMENT.COLOUR  = "#EA6A58",
        LOSSES.SEGMENT.COLOUR = "#94A6F1",
        GAINS.HIST.COLOUR     = "pink",
        LOSSES.HIST.COLOUR    = "skyblue",
        CHROM.FILL.COLOUR     = "#E1E1E1",
        CHROM.BORDER.COLOUR   = "black",
        CHROM.STRIPE.COLOUR   = "black",
        CHROM.WIDTH           = 2,
        CHROM.BUFFER          = 1,
        SEGMENT.LWD           = 0.5,
        MIN.X                 = -40,
        MAX.X                 = 40
    )

    if (is.null(plot_options)) {
        plot_options <- default.plot.options
    }

    for (key in names(default.plot.options)) {
        if (!(key %in% names(plot_options))) {
            plot_options[[key]] <- default.plot.options[[key]]
        }
    }

    if (get_plot_options) {
        return(default.plot.options)
    }

    CHR <- unique(dt_gains$seqnames)
    if(length(CHR) == 0) {
        CHR <- chr_override
    } else {
        if(length(CHR) > 1) stop("More than one chromosome in input data")
        if(nrow(dt_losses) > 0 && unique(dt_losses$seqnames) != CHR) stop("Chromosomes do not match between tables")
    }

    # Depth histograms will be plotted underneath the segments
    if (nrow(dt_gains) > 0) {
        histo_gains <- coverage_histogram(dt_gains)
    } else {
        histo_gains <- data.table(start=1, end=1, N=1)[0]
    }

    if (nrow(dt_losses) > 0) {
        histo_losses <- coverage_histogram(dt_losses)
    } else {
        histo_losses <- data.table(start=1, end=1, N=1)[0]
    }

    # Make sure losses are measured as negative, so they go on left side of plot
    if (nrow(dt_losses) > 0) {
        dt_losses[, layer := -abs(layer)]
        histo_losses[, N := -abs(N)]
    }

    # Get chromosome length
    L <- chr_lengths[CHROM==paste0("Chr", tolower(CHR)), LENGTH]

    # Constrain segments to the end of the chromosome
    dt_gains[end > L, end := L]
    dt_losses[end > L, end := L]

    # Open a new plot
    max_loss_depth <- if (nrow(histo_losses) > 0) histo_losses[, min(N)] - 1 else -1
    max_gain_depth <- if (nrow(histo_gains) > 0) histo_gains[, max(N)] + 1 else 1
    plot(c(max_loss_depth, max_gain_depth),
         c(0, L),
         type = "n", ylab = "Position", xlab = "Coverage",
         xaxt = "n", yaxt = "n",
         bty = "n",
         xlim = c(plot_options$MIN.X, plot_options$MAX.X),
         ylim = c(ifelse(CHR %in% c("5", "6", "X"),
                         -3.0e8,
                         -6.2e8), 2e7))

    # Add a y-axis
    if (CHR %in% c("1", "5")) {
        axis(side = 2, at=seq(0, -(L + 1e8 - L %% 1e8), -1e8),
             labels=abs(seq(0, -(L + 1e8 - L %% 1e8), -1e8)), line = -2)
    }


    # Plot the chromosome
    chrom_halfwidth = plot_options$CHROM.WIDTH / 2 # Half the width of the chromosome rectangle
    chrom_buffer = plot_options$CHROM.BUFFER # additional buffer space between the chromosome and any plotted elements
    chrom_spacing = chrom_halfwidth + chrom_buffer # Total offset used to position plot elements away from central chromosome
    if (chrom_halfwidth == 0) {
        segments(0, L, lwd = 2, lend = 1,
                 col = plot_options$CHROM.BORDER.COLOUR)
    } else {
        rect(-chrom_halfwidth, 0, chrom_halfwidth, -L,
             border = NA,
             col = plot_options$CHROM.FILL.COLOUR)
        if (nrow (dt_gains) > 0) {
            segments(-chrom_halfwidth, unique(dt_gains[start < L & end <= L, c(-start, -end)]), chrom_halfwidth,
                     col = plot_options$CHROM.STRIPE.COLOUR, lend = 1)
        }
        if (nrow (dt_losses) > 0) {
            segments(-chrom_halfwidth, unique(dt_losses[start < L & end <= L, c(-start, -end)]), chrom_halfwidth,
                     col = plot_options$CHROM.STRIPE.COLOUR, lend = 1)
        }
        rect(-chrom_halfwidth, 0, chrom_halfwidth, -L,
             border = plot_options$CHROM.BORDER.COLOUR,
             col = NA)
        text(0, 2e7, CHR, cex = 1.2, font = 2)
    }

    # Gains histogram
    if (nrow (histo_gains) > 0) {
        rect(chrom_spacing,
             -histo_gains$start,
             chrom_spacing + histo_gains$N,
             -histo_gains$end,
             col = plot_options$GAINS.HIST.COLOUR, border = NA)
    }

    # Losses histogram
    if (nrow(histo_losses) > 0) {
        rect(-chrom_spacing,
             -histo_losses$start,
             histo_losses$N - chrom_spacing,
             -histo_losses$end,
             col = plot_options$LOSSES.HIST.COLOUR, border = NA)
    }

    # Plot gains segments
    # 1) Decide on the fill colour for the segments. If category colours have been specified,
    #    and there is a Category annotation in the data, then prioritise this colour.
    #    Otherwise, or if there is any error or missing value in the category colours,use
    #    whatever is in the plot options list for GAINS.SEGMENT.COLOUR.
    if (nrow(dt_gains) > 0 & "Category" %in% colnames(dt_gains) & !is.null(category_colours)) {
        gains_fill_colours <- category_colours[dt_gains$Category]
        gains_fill_colours[sapply(gains_fill_colours, is.null)] <- plot_options$GAINS.SEGMENT.COLOUR
        gains_fill_colours <- unlist(gains_fill_colours)
        stopifnot(length(gains_fill_colours) == nrow(dt_gains))
    } else {
        gains_fill_colours <- plot_options$GAINS.SEGMENT.COLOUR
    }
    # 2) Plot the segments using rect
    rect(dt_gains$layer - 0.8 + chrom_spacing,
         -dt_gains$start,
         dt_gains$layer - 0.2 + chrom_spacing,
         -dt_gains$end, col = gains_fill_colours,
         lwd = plot_options$SEGMENT.LWD)

    if ("mip_position" %in% colnames(dt_gains)) {
        points(dt_gains[, .(layer -0.5 + chrom_spacing, -mip_position)], pch = "*", col = "goldenrod")
    }

    # Losses segments
    # 1) Select the fill colour(s)
    if (nrow(dt_losses) > 0 & "Category" %in% colnames(dt_losses) & !is.null(category_colours)) {
        losses_fill_colours <- category_colours[dt_losses$Category]
        losses_fill_colours[sapply(losses_fill_colours, is.null)] <- plot_options$LOSSES.SEGMENT.COLOUR
        losses_fill_colours <- unlist(losses_fill_colours)
        stopifnot(length(losses_fill_colours) == nrow(dt_losses))
    } else {
        losses_fill_colours <- plot_options$LOSSES.SEGMENT.COLOUR
    }
    # 2) Plot the segments with rect
    rect(dt_losses$layer + 0.8 - chrom_spacing,
         -dt_losses$start,
         dt_losses$layer + 0.2 - chrom_spacing,
         -dt_losses$end, col = losses_fill_colours,
         lwd = plot_options$SEGMENT.LWD)

    if ("mip_position" %in% colnames(dt_losses)) {
        points(dt_losses[, .(layer + 0.5 - chrom_spacing, -mip_position)], pch = "*", col = "goldenrod")
    }

    axlim <- if (nrow(dt_losses) > 0) dt_losses[, min(layer)] else -1
    axlim <- axlim + 5 - (axlim %% 5)
    axis(side = 1,
         at = seq(axlim, 0, 5) - chrom_spacing,
         labels = abs(seq(axlim, 0, 5)),
         cex.axis = .8, lwd = .5,
         pos = -L-2e7,
         col = plot_options$LOSSES.SEGMENT.COLOUR,
         col.ticks = plot_options$LOSSES.SEGMENT.COLOUR,
         col.axis = plot_options$LOSSES.SEGMENT.COLOUR)

    axlim <- if (nrow(dt_gains) > 0) dt_gains[, max(layer)] else 1
    axlim <- axlim - (axlim %% 5)
    axis(side = 1,
         at = seq(0, axlim, 5) + chrom_spacing,
         labels = seq(0, axlim, 5),
         pos = -L - 2e7,
         cex.axis = .8, lwd = 1,
         col = plot_options$GAINS.SEGMENT.COLOUR,
         col.ticks = plot_options$GAINS.SEGMENT.COLOUR,
         col.axis = plot_options$GAINS.SEGMENT.COLOUR)
}
