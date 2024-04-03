#' Get three prime position
#'
#' A function that from a GRanges
#' gives the 3' position
#' @param input_gr A GRanges or GRangeList
#' @return A vector of integers
three_prime_pos <- function(input_gr) {
    return(GenomicRanges::start(
        GenomicRanges::resize(input_gr, width = 1, fix = "end")
    ))
}

#' Get five prime position
#'
#' A function that from a GRanges
#' gives the 5' position
#' @param input_gr A GRanges or GRangeList
#' @return A vector of integers
five_prime_pos <- function(input_gr) {
    return(GenomicRanges::start(
        GenomicRanges::resize(input_gr, width = 1, fix = "start")
    ))
}

#' Extract last exons
#'
#' A function that from a GRanges from gtf
#' select only entries for the last exons
#' If multiple exons overlap the last base of the groupping_variable,
#' they will all be reported.
#' @param input_gr A GRanges from a gtf
#' @param groupping_variable A string with the name of the metadata
#'                           which should be used to group
#' @param invert A boolean that indicates if you want all except the last exons
#' @return A GRanges which contains a subset of `input_gr`
extract_last_exon <- function(
        input_gr, groupping_variable = "transcript_id", invert = FALSE
) {
    # I only work with exons
    exons_gr <- base::subset(input_gr, type == "exon")
    # I split by groupping_variable
    exons_gr_group <- GenomicRanges::mcols(exons_gr)[, groupping_variable]
    exons_gr_split <-
        GenomicRanges::split(
            exons_gr,
            exons_gr_group
        )
    # I get the three prime extremity
    group_range <- unlist(range(exons_gr_split))
    group_extremity <- three_prime_pos(group_range)
    names(group_extremity) <- names(group_range)
    indices_last_exons <- which(
        three_prime_pos(exons_gr) ==
            group_extremity[exons_gr_group]
    )
    # Return the corresponding GRanges object
    if (invert) {
        return(exons_gr[-indices_last_exons])
    } else {
        return(exons_gr[indices_last_exons])
    }
}

#' Overlap exons and extend three prime end
#'
#' A function that from 2 GRanges returns a subset of the first
#' GRanges which have been extended using the second GRanges
#' @param input_gr_to_extend A GRanges with exons to extend
#'                           (strand * are excluded)
#' @param input_gr_to_overlap A GRanges with intervals to overlap
#' @return A GRanges which is a subset of `input_gr_to_extend` where
#'         3' end have been modified to match the 3' end of
#'         `input_gr_to_overlap` if they overlap
#'         (initial width have been stored into old_width)
#' @importFrom GenomicRanges strand
extend_using_overlap <- function(input_gr_to_extend, input_gr_to_overlap) {
    # Remove strands which are not in + - and non exonic features
    input_gr_to_extend <- subset(
        input_gr_to_extend,
        type == "exon" &
            as.character(strand(input_gr_to_extend)) != "*"
    )
    # Remove non-exonic features from input_gr_to_overlap
    input_gr_to_overlap <- subset(input_gr_to_overlap, type == "exon")
    rlang::inform(
        paste("Compute overlap between", length(input_gr_to_extend),
              "exons and", length(input_gr_to_overlap), "exons.\n")
    )
    # Compute overlap:
    full_overlap <- suppressWarnings(GenomicRanges::findOverlaps(
        input_gr_to_extend,
        input_gr_to_overlap
    ))
    progression_msg("Split by overlap.\n")
    # Group all overlaps into a single range
    input_gr_to_overlap_intersect <- GenomicRanges::split(
        input_gr_to_overlap[S4Vectors::subjectHits(full_overlap)],
        S4Vectors::queryHits(full_overlap)
    )
    input_gr_to_overlap_intersect_range <-
        unlist(range(input_gr_to_overlap_intersect))
    # Select the corresponding to_extend
    input_gr_to_extend_intersect <-
        input_gr_to_extend[unique(S4Vectors::queryHits(full_overlap))]
    # Calculate the potential_extension
    potential_extention <- ifelse(
        strand(input_gr_to_extend_intersect) == "+",
        GenomicRanges::end(input_gr_to_overlap_intersect_range) -
            GenomicRanges::end(input_gr_to_extend_intersect),
        GenomicRanges::start(input_gr_to_extend_intersect) -
            GenomicRanges::start(input_gr_to_overlap_intersect_range)
    )
    # Select only true extension
    input_gr_to_extend_intersect_selected <-
        input_gr_to_extend_intersect[potential_extention > 0]
    potential_extention_selected <-
        potential_extention[potential_extention > 0]
    # Store old width
    input_gr_to_extend_intersect_selected$old_width <-
        GenomicRanges::width(input_gr_to_extend_intersect_selected)
    # Extend
    input_gr_to_extend_intersect_selected_extended <-
        GenomicRanges::resize(
            input_gr_to_extend_intersect_selected,
            width =
                input_gr_to_extend_intersect_selected$old_width +
                potential_extention_selected
        )
    return(input_gr_to_extend_intersect_selected_extended)
}

#' Get overlaps from different genes
#'
#' A function that from 2 GRanges generates a dataframe
#' With queryHits, subjectHits when the gene_id is different
#' @param gr1 A GRanges with 'gene_id'
#' @param gr2 A GRanges with 'gene_id'
#' @return a data.frame with overlaps between gr1 and gr2 when gene_id from gr1
#'         is different from gene_id from gr2. The data.frame has 4 columns:
#'         `queryHits`, `subjectHits`, `query_gene_id` and `subject_gene_id`
overlap_different_genes <- function(gr1, gr2) {
    # Compute overlap
    ov <- suppressWarnings(GenomicRanges::findOverlaps(gr1, gr2))
    ov_df <- as.data.frame(ov)
    ov_df$query_gene_id <- gr1$gene_id[ov_df$queryHits]
    ov_df$subject_gene_id <- gr2$gene_id[ov_df$subjectHits]
    return(subset(ov_df, query_gene_id != subject_gene_id))
}



#' Adjust for collision
#'
#' A function that from a GRanges with 'old_width'
#' Change the starts and ends to prevent collisions larger than
#' with old coordinates
#' @param input_gr A GRanges with 1 meta: 'old_width'
#' @return A list with:
#'         - 'pot_issues': A dataframe with exons which overlaps between
#'                        `input_gr` and itself while gene_ids are different
#'         - 'new_gr': A GRanges identical to `input_gr` except that start/end
#'                     have been adjusted to prevent collisions.
adjust_for_collision <- function(input_gr) {
    # Use an 'id':
    input_gr$id <- paste0(input_gr$exon_id, "_", input_gr$transcript_id)
    # Select entries which have been extended:
    extended_gr <- subset(input_gr, old_width < width)
    # Select only the extension
    extension_gr <- GenomicRanges::resize(
        extended_gr,
        GenomicRanges::width(extended_gr) - extended_gr$old_width,
        fix = "end"
    )
    # Select the first base of extension
    first_base_extension_gr <-
        GenomicRanges::resize(
            extension_gr,
            width = 1
        )
    # Generate the original input
    original_gr <-
        GenomicRanges::resize(
            input_gr,
            width = input_gr$old_width
        )

    # Compute overlap between first base of extension and original
    # and get those from different gene_ids:
    pot_issues_first <- overlap_different_genes(first_base_extension_gr,
                                                original_gr)
    pot_issues_first$query_id <-
        first_base_extension_gr[pot_issues_first$queryHits]$id
    # These should not be extended
    # Restore them
    to_restore <- which(input_gr$id %in% pot_issues_first$query_id)
    input_gr[to_restore] <-
        original_gr[to_restore]
    # Store this into a data frame
    pot_issues_df <- as.data.frame(input_gr[to_restore])
    if (nrow(pot_issues_df) > 0) {
        pot_issues_df$issue <- "first_base_extension_overlaps_different_gene"
        pot_issues_df$solution <- "no extension"
        pot_issues_df$collision_base <- NA
    }
    # Select entries which have been extended:
    extended_gr <- subset(input_gr, old_width < width)
    # Select only the extension
    extension_gr <- GenomicRanges::resize(
        extended_gr,
        GenomicRanges::width(extended_gr) - extended_gr$old_width,
        fix = "end"
    )
    # To deal with tail-to-head extensions we
    # compute overlap between extension and five prime of input_gr
    # and get those from different gene_ids:
    pot_issues <- overlap_different_genes(extension_gr,
                                          GenomicRanges::resize(
                                              input_gr,
                                              width = 1
                                          ))

    if (nrow(pot_issues) > 0) {
        overlapped_range <- unlist(range(
            GenomicRanges::split(input_gr[pot_issues$subjectHits],
                                 pot_issues$queryHits)
        ))
        # Generate a GRanges with issues:
        extended_gr_issues <- extended_gr[as.numeric(names(overlapped_range))]
        extended_gr_issues$collision_base <- five_prime_pos(overlapped_range)
        # Compute the new_widths to avoid tail-to-head
        new_widths <- abs(extended_gr_issues$collision_base -
                              five_prime_pos(extended_gr_issues))
        names(new_widths) <- extended_gr_issues$id
        # Modify the input_gr
        input_gr[match(names(new_widths), input_gr$id)] <-
            GenomicRanges::resize(
                input_gr[match(names(new_widths), input_gr$id)],
                new_widths
            )
        # Store solving into a data.frame
        pot_issues_df_temp <- as.data.frame(extended_gr_issues)
        pot_issues_df_temp$issue <- "tail_to_head_collision"
        pot_issues_df_temp$solution <- "smaller extension"
        pot_issues_df <- rbind(pot_issues_df, pot_issues_df_temp)
    }
    # Now we compute overlap between extension and input_gr
    # We can restrict to exons with issues:
    new_extended_gr <- subset(input_gr, old_width < width)
    # Select only the extension
    new_extension_gr <- GenomicRanges::resize(
        new_extended_gr,
        GenomicRanges::width(new_extended_gr) - new_extended_gr$old_width,
        fix = "end"
    )
    new_pot_issues <- overlap_different_genes(
        new_extension_gr,
        input_gr
    )

    if (nrow(new_pot_issues) > 0) {
        # The only case that is remaining:
        # old
        # ----->
        #  ---->
        # new
        # ------------->
        #  ------------>
        new_pot_issues$query_id <-
            new_extension_gr$id[new_pot_issues$queryHits]
        to_restore <- which(input_gr$id %in% new_pot_issues$query_id)
        input_gr[to_restore] <-
            original_gr[to_restore]
        pot_issues_df_temp <-
            as.data.frame(new_extension_gr[new_pot_issues$queryHits])
        pot_issues_df_temp$issue <- "both_equal_extension"
        pot_issues_df_temp$solution <- "no extension"
        pot_issues_df_temp$collision_base <- NA
        pot_issues_df <- rbind(pot_issues_df, pot_issues_df_temp)
    }
    # Remove the id:
    input_gr$id <- NULL
    return(list(pot_issues = pot_issues_df,
                new_gr = input_gr))
}

#' Filter new exons for collision
#'
#' A function that from 2 GRanges filter exons from the first one
#' so they do not go three prime to the first collision with the second one.
#' @param all_exons_interesting A GRanges with exons to trim and filter
#' @param input_gr_to_extend A GRanges to overlap
#' @return A GRanges subset of `all_exons_interesting`
filter_new_exons <- function(all_exons_interesting, input_gr_to_extend) {
    progression_msg(
        "Check added exons which would overlap existing annotations.\n"
    )
    # Check that the new exons do not overlap existing annotations
    # ov_new_exons is the overlap between:
    # queries = all_exons_interesting
    # subjects = input_gr_to_extend
    # We arbitrarily decide to not include new exons if they overlap with any
    # existing exons (including those with the same gene_id).
    ov_new_exons <- suppressWarnings(GenomicRanges::findOverlaps(
        all_exons_interesting,
        input_gr_to_extend
    ))

    input_gr_to_extend_intersected <-
        input_gr_to_extend[S4Vectors::subjectHits(ov_new_exons)]

    # Get the five prime end of the intersected by new transcript_id
    input_gr_to_extend_intersected_range <-
        unlist(range(
            GenomicRanges::split(
                input_gr_to_extend_intersected,
                all_exons_interesting$transcript_id[
                    S4Vectors::queryHits(ov_new_exons)
                ]
            )
        ))
    collision_base_per_transcript <- five_prime_pos(
        input_gr_to_extend_intersected_range
    )
    names(collision_base_per_transcript) <-
        names(input_gr_to_extend_intersected_range)
    all_exons_interesting$collision_base <-
        collision_base_per_transcript[all_exons_interesting$transcript_id]
    progression_msg("Remove new candidate exons",
                    " three prime of overlap with existing annotations.\n")
    # Filter exons after collision:
    all_exons_interesting <-
        subset(all_exons_interesting,
               is.na(collision_base) |
                   (strand == "+" & start < collision_base) |
                   (strand == "-" & end > collision_base)
        )
    rlang::inform(
        paste("Stay", length(unique(all_exons_interesting$id)), "candidate",
              "exons that may be included into",
              length(unique(all_exons_interesting$transcript_id)),
              "transcripts.\n")
    )
    progression_msg(
        paste("Adjust start/end of candidate exons to not",
              "overlap with existing annotations.\n")
    )

    if (length(all_exons_interesting) == 0) {
        return()
    }

    to_adjust <- which(!is.na(all_exons_interesting$collision_base) &
                           (GenomicRanges::start(all_exons_interesting) <=
                                all_exons_interesting$collision_base) &
                           (GenomicRanges::end(all_exons_interesting) >=
                                all_exons_interesting$collision_base))

    # Compute the new_widths to avoid tail-to-head
    new_widths <- abs(all_exons_interesting[to_adjust]$collision_base -
                          five_prime_pos(all_exons_interesting[to_adjust]))
    # Apply to all_exons_interesting
    all_exons_interesting[to_adjust] <- GenomicRanges::resize(
        all_exons_interesting[to_adjust],
        new_widths
    )
    return(all_exons_interesting)
}

#' Add new exons
#'
#' A function that from 2 GRanges add exons from the second one
#' to the first one
#' if the 3p of the last exon of the transcript in the first GRanges
#' matches the 3p of an exon in the second one
#' @param input_gr_to_extend A GRanges to be complemented
#' @param input_gr_with_new_exons A GRanges with exons to be added
#' to the first one (exons with strand '*' are excluded)
#' @importFrom GenomicRanges strand
#' @return A GRanges identical to `input_gr_to_extend` with new exons
#'         whose `exon_id` contains BREW3R. `exon_number` may have changed.
#' @details
#' Potential new exons will be filtered for collision with exons present
#' in the first GRanges even if they belong to the same gene_id.
#' For the moment all potential exons extensions are added to the
#' same existing transcript_id so introns maybe artificial introns.
#'
add_new_exons <- function(input_gr_to_extend, input_gr_with_new_exons) {
    # Subset input_gr_with_new_exons
    input_gr_with_new_exons <- subset(
        input_gr_with_new_exons,
        strand(input_gr_with_new_exons) != "*"
    )
    # Add an id:
    input_gr_with_new_exons$id <- seq_along(input_gr_with_new_exons)
    progression_msg("Extract last exon of first GRanges.\n")
    # Get the last exons of input_gr_to_extend
    input_gr_last <- extract_last_exon(input_gr_to_extend)
    # First find exons which 'ends' at the same position
    progression_msg("Compute overlap with second GRanges.\n")
    last_base_input_gr_last <-
        GenomicRanges::resize(
            input_gr_last, width = 1, fix = "end"
        )
    last_base_input_gr_with_new_exons <-
        GenomicRanges::resize(
            input_gr_with_new_exons, width = 1, fix = "end"
        )
    ov <- suppressWarnings(GenomicRanges::findOverlaps(
        last_base_input_gr_last,
        last_base_input_gr_with_new_exons
    ))
    ov_df <- as.data.frame(ov)
    if (nrow(ov_df) == 0) {
        return(input_gr_to_extend)
    }
    progression_msg("Only keep those with exons to be added.\n")
    # First get transcript_id from exon overlapped
    ov_df$subject_transcript_id <-
        input_gr_with_new_exons$transcript_id[ov_df$subjectHits]

    # Get transcript_extremity and transcript_id from input_gr_with_new_exons
    input_gr_with_new_exons_subset <- subset(
        input_gr_with_new_exons,
        transcript_id %in% ov_df$subject_transcript_id
    )
    input_gr_with_new_exons_subset_range <-
        unlist(range(
            GenomicRanges::split(input_gr_with_new_exons_subset,
                                 input_gr_with_new_exons_subset$transcript_id)
        ))
    # Remove overlaps where it is the transcript subject extremity
    progression_msg(
        "Remove overlaps where transcript 3' matches between both\n"
    )
    ov_df <- subset(
        ov_df,
        GenomicRanges::start(last_base_input_gr_last[ov_df$queryHits]) !=
            three_prime_pos(
                input_gr_with_new_exons_subset_range[
                    ov_df$subject_transcript_id
                ]
            )
    )
    if (nrow(ov_df) == 0) {
        return(input_gr_to_extend)
    }
    progression_msg("Attribute exons to be added to the good transcript.\n")
    # Get all exons from subject_transcript_id
    # which are after the query_extremity
    # First get all exons of subject_transcript_id
    # (once per queryHits)
    input_gr_with_new_exons_subset <- subset(
        input_gr_with_new_exons,
        transcript_id %in% ov_df$subject_transcript_id
    )
    input_gr_with_new_exons_subset_split <-
        GenomicRanges::split(input_gr_with_new_exons_subset,
                             input_gr_with_new_exons_subset$transcript_id)
    nb_exons_per_transcript <-
        table(input_gr_with_new_exons_subset$transcript_id)
    all_pot_exons <-
        unlist(
            input_gr_with_new_exons_subset_split[ov_df$subject_transcript_id]
        )
    all_pot_exons$queryHits <- rep(
        ov_df$queryHits, nb_exons_per_transcript[ov_df$subject_transcript_id]
    )
    all_pot_exons$three_prime_query <- GenomicRanges::start(
        last_base_input_gr_last[all_pot_exons$queryHits]
    )
    # Filter, I don't know how to do better
    all_exons_interesting <- subset(
        all_pot_exons,
        ifelse(strand == "+",
               start > three_prime_query,
               end < three_prime_query)
    )
    rlang::inform(
        paste("Found", length(unique(all_exons_interesting$id)), " exons",
              "that may be included into",
              length(unique(all_exons_interesting$queryHits)),
              "transcripts.\n")
    )
    # Attribute the new transcript_id to new exons
    all_exons_interesting$transcript_id <-
        last_base_input_gr_last[all_exons_interesting$queryHits]$transcript_id
    # Filter
    all_exons_interesting <-
        filter_new_exons(all_exons_interesting, input_gr_to_extend)
    if (is.null(all_exons_interesting)) {
        return(input_gr_to_extend)
    }
    # Here I do something which is not ideal:
    # If multiple transcripts could contribute to extension
    # they will be merged into a single one:
    # For example:
    # gene_1:           --->
    # to_overlap_tr_A:  ----/    \-->
    # to_overlap_tr_B:  ----/          \--->
    # Will give
    # gene_1:           ----/    \----/\--->
    # While the second intron is not a real intron...

    progression_msg("Combine exons extending the same transcript.\n")
    all_exons_to_add_gr <- unlist(
        GenomicRanges::reduce(
            GenomicRanges::split(
                all_exons_interesting,
                all_exons_interesting$transcript_id
            )
        )
    )
    all_exons_to_add_gr$transcript_id <- names(all_exons_to_add_gr)
    names(all_exons_to_add_gr) <- NULL

    rlang::inform(
        paste("Finally", length(all_exons_to_add_gr), "combined exons",
              "will be included into",
              length(unique(all_exons_to_add_gr$transcript_id)),
              "transcripts.\n")
    )
    progression_msg("Annotate new exons with transcript info.\n")
    # Add transcript_id infos:
    mcols_input <- as.data.frame(GenomicRanges::mcols(
        input_gr_last[unique(all_exons_interesting$queryHits)]
    ))
    transcript_annotation <-
        unique(
            mcols_input[,
                        c(grep("^source$", colnames(mcols_input), value = TRUE),
                          grep("^gene_", colnames(mcols_input), value = TRUE),
                          grep("^transcript_", colnames(mcols_input),
                               value = TRUE))
                        ]
        )

    GenomicRanges::mcols(all_exons_to_add_gr) <-
        transcript_annotation[match(all_exons_to_add_gr$transcript_id,
                                    transcript_annotation$transcript_id), ]
    all_exons_to_add_gr$type <- "exon"

    progression_msg("Annotate new exons with exon_id.\n")
    # I add exon_id
    if (any(grepl("BREW3R", input_gr_to_extend$exon_id))) {
        last.id.used <- max(
            as.numeric(
                gsub(
                    gsub(
                        grep("BREW3R",
                             input_gr_to_extend$exon_id,
                             value = TRUE),
                        pattern = "BREW3R", replacement = ""
                    ),
                    pattern = "[^0-9]", replacement = ""
                )
            )
        )
    } else {
        last.id.used <- 0
    }
    all_exons_to_add_gr$exon_id <-
        paste0("BREW3R",
               sprintf("%010d",
                       last.id.used + seq_along(all_exons_to_add_gr)
               )
        )


    progression_msg("Update exon_number.\n")
    # Update exon_number
    # There are 2 exon_number modes:
    # 1. from 5p to 3p (I call it transcript)
    # 2. from start to end whatever the strand is (I call it coordinate)
    # To get the mode I select a second exon with strand -
    if (!"exon_number" %in%
        colnames(GenomicRanges::mcols(input_gr_to_extend))) {
        my_mode <- "transcript"
        gr_to_annotate <- suppressWarnings(
            c(all_exons_to_add_gr, input_gr_to_extend)
        )
        gr_to_add <- NULL
    } else {
        gr_to_annotate <- suppressWarnings(c(
            all_exons_to_add_gr,
            subset(
                input_gr_to_extend,
                transcript_id %in% all_exons_to_add_gr$transcript_id
            )
        ))
        gr_to_add <- subset(
            input_gr_to_extend,
            !transcript_id %in% all_exons_to_add_gr$transcript_id
        )
        potential_selected_exon <- subset(input_gr_to_extend,
                                          exon_number == 2 & strand == "-")
        if (length(potential_selected_exon) == 0) {
            my_mode <- "transcript"
        } else {
            selected_exon <- potential_selected_exon[1]
            first_corresponding_exon <-
                subset(input_gr_to_extend,
                       transcript_id == selected_exon$transcript_id &
                           exon_number == 1)
            if (GenomicRanges::start(selected_exon) <
                GenomicRanges::start(first_corresponding_exon)) {
                my_mode <- "transcript"
            } else {
                my_mode <- "coordinate"
            }
        }
    }
    # Order by transcript_id first:
    gr_to_annotate <- gr_to_annotate[
        order(gr_to_annotate$transcript_id,
              GenomicRanges::start(gr_to_annotate),
              GenomicRanges::end(gr_to_annotate))
    ]
    if (my_mode == "transcript") {
        gr_to_annotate_plus <-
            subset(gr_to_annotate,
                   strand(gr_to_annotate) == "+")
        gr_to_annotate_minus <-
            subset(gr_to_annotate,
                   strand(gr_to_annotate) == "-")
        gr_to_annotate <- c(
            gr_to_annotate_plus,
            rev(gr_to_annotate_minus)
        )
    }
    exons_per_transcript <- table(gr_to_annotate$transcript_id)
    exon_numbers_list <- lapply(exons_per_transcript, seq_len)
    exon_numbers <-
        unlist(exon_numbers_list[unique(gr_to_annotate$transcript_id)])
    gr_to_annotate$exon_number <- unname(exon_numbers)
    new_gr <- c(gr_to_annotate, gr_to_add)

    return(new_gr)
}
