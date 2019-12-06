# Some functions to pack CNV intervals into layers for plotting chromosome maps.
# There are 4 versions of the same function - rectangle_packing4 is the only one that
# should be used (and is the only one exported into the package).
# The others are slow (1 and 2) or have unnecessary dependencies (3).

#' Algorithm to pack segments into layers for plotting. Greedily assigns widest segments first.
rectangle_packing <- function(dt) {
    dt[, layer := 0]
    setkey(dt, seqnames, start, end)
    remaining <- dt$Unique.ID
    total <- length(remaining)
    assigned <- 0

    # Pick out the longest remaining segment and remove from remaining list
    current <- dt[Unique.ID %in% remaining][width == max(width), Unique.ID][1]
    remaining <- remaining[remaining != current]

    # First time - all layers are empty - assign to bottom layer
    if (dt[, max(layer)] == 0) {
        dt[Unique.ID == current, layer := 1]
        assigned <- assigned + 1
    }

    while (length(remaining) > 0) {
        # Pick out the longest remaining segment and remove from remaining list
        current <- dt[Unique.ID %in% remaining][width == max(width), Unique.ID][1]
        remaining <- remaining[remaining != current]

        for (curr_layer in 1:dt[, max(layer)]) {
            dt_layer <- dt[layer == curr_layer]

            ol <- foverlaps(dt[Unique.ID == current], dt[layer == curr_layer], nomatch = 0L)
            # Doesn't overlap with current layer, so assign here
            if (ol[, .N] == 0) {
                dt[Unique.ID == current, layer := curr_layer]
                assigned <- assigned + 1
                message <- paste("Assigned", current, " to layer", curr_layer, "(", assigned, "/", total, ")")
                break
            }
        }

        # Couldn't find a compatible layer, so assign to a new layer
        if (dt[Unique.ID == current, layer] == 0) {
            new_layer <- dt[, max(layer) + 1]
            dt[Unique.ID == current, layer := new_layer]
            assigned <- assigned + 1
            message <- paste("Assigned", current, " to new layer", new_layer, "(", assigned, "/", total, ")")
        }
        dt[layer>0]
        print(message)
    }
    dt
}

#' Algorithm to pack segments into layers for plotting. Greedily assigns widest segments first.
rectangle_packing2 <- function(dt) {
    dt[, layer := 0]
    setkey(dt, seqnames, start, end)
    remaining <- dt$Unique.ID
    total <- length(remaining)
    assigned <- 0

    # Incompatibilities
    ol <- foverlaps(dt,dt)[Unique.ID != i.Unique.ID]
    incomp <- list()
    for (id in ol$Unique.ID) {
        incomp[[id]] <- ol[Unique.ID == id, i.Unique.ID]
    }

    # Pick out the longest remaining segment and remove from remaining list
    current <- dt[Unique.ID %in% remaining][width == max(width), Unique.ID][1]
    remaining <- remaining[remaining != current]

    # First time - all layers are empty - assign to bottom layer
    if (dt[, max(layer)] == 0) {
        dt[Unique.ID == current, layer := 1]
        assigned <- assigned + 1
    }

    while (length(remaining) > 0) {
        # Pick out the longest remaining segment and remove from remaining list
        current <- dt[Unique.ID %in% remaining][width == max(width), Unique.ID][1]
        remaining <- remaining[remaining != current]

        for (curr_layer in 1:dt[, max(layer)]) {
            dt_layer <- dt[layer == curr_layer]

            already_placed <- dt[layer == curr_layer, Unique.ID]
            incompatible_with_current <- incomp[[current]]

            # If no intersection, then the current segment is allowed in the current layer
            allowed_in_layer <- length(intersect(incompatible_with_current, already_placed)) == 0

            # No incomptibility with current layer, so assign here
            if (allowed_in_layer) {
                dt[Unique.ID == current, layer := curr_layer]
                assigned <- assigned + 1
                message <- paste("Assigned", current, " to layer", curr_layer, "(", assigned, "/", total, ")")
                break
            }
        }

        # Couldn't find a compatible layer, so assign to a new layer
        if (dt[Unique.ID == current, layer] == 0) {
            new_layer <- dt[, max(layer) + 1]
            dt[Unique.ID == current, layer := new_layer]
            assigned <- assigned + 1
            message <- paste("Assigned", current, " to new layer", new_layer, "(", assigned, "/", total, ")")
        }
        dt[layer>0]
        print(message)
    }
    dt
}

#' A faster algorithm to pack aligned reads into layers
rectangle_packing3 <- function(dt, sortby = "width") {
    if (nrow(dt) == 0) return (dt)
    require(IRanges)

    # Function to assign the lowest available layer, avoiding
    # any of the layers in `disallowed`
    assignLayer <- function(disallowed) {
        if (length(disallowed) == 0 || (length(disallowed) == 1 && disallowed == 0)) return (1)
        M <- max(disallowed)
        candidate <- setdiff(1:M, disallowed)
        ifelse(length(candidate) == 0, M+1, candidate)
    }

    # Take advantage of fast findOverlaps function in IRanges
    # NB, data is sorted so that widest segments are assigned first - a greedy
    # strategy to try to make the assignment look nicer
    if (sortby == "width") {
        ir <- dt[order(-width), IRanges(start = start, end = end, names = Unique.ID)]
    } else {
        ir <- dt[order(start), IRanges(start = start, end = end, names = Unique.ID)]
    }
    mcols(ir)$layer <- as.integer(0)
    hits <- findOverlaps(ir, drop.self = TRUE)

    # No conflicts possible for first segment, so place in bottom layer
    mcols(ir[1, ])$layer <- 1

    N <- dt[, .N]
    if (N > 1) {
        for (n in 2:N) {
            # Search for any segments that overlap the current segment. If any layers are assigned
            # for the hits, then they are disallowed for our current segment.
            disallowed_layers <- unique(mcols(ir[subjectHits(hits)[queryHits(hits)==n], ])$layer)
            assigned_layer <- assignLayer(disallowed_layers)
            mcols(ir[n, ])$layer <- assigned_layer
            #print(paste("(", n, "/", N, ")"))
        }
    }

    # Transfer the layer info back from the IRanges to the data.table
    dt[order(match(dt$Unique.ID, as.character(names(ir)))), layer := as.integer(mcols(ir)$layer)]
    dt
}


#' A faster algorithm to pack aligned reads into layers. Doesn't need IRanges.
#' @importFrom "data.table" setorder
#' @export
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

    # Put the dt back how we found it
    setorder(dt, origOrder)
    dt[, origOrder := NULL]
    dt
}
