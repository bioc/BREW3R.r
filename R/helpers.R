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
#' @importFrom dplyr %>%
extract_last_exon <- function(
        input_gr, groupping_variable = "transcript_id", invert = FALSE
) {
    # I only work with exons
    exons_gr <- base::subset(input_gr, type == "exon")

    # Get the last exons row numbers
    indices_last_exons <- as.data.frame(exons_gr) %>%
        dplyr::mutate(row_number = base::seq_len(dplyr::n())) %>%
        dplyr::group_by(!!!dplyr::syms(groupping_variable)) %>%
        dplyr::mutate(is_last = ifelse(strand == "+",
                                       end == max(end),
                                       start == min(start))) %>%
        dplyr::ungroup() %>%
        dplyr::filter(is_last) %>%
        dplyr::select(row_number)
    # Return the corresponding GRanges object
    if (invert) {
        return(exons_gr[-indices_last_exons$row_number])
    } else {
        return(exons_gr[indices_last_exons$row_number])
    }
}

#' Overlap exons and extend three prime end
#'
#' A function that from 2 GRanges returns a subset of the first
#' GRanges which have been extended using the second GRanges
#' @param input_gr_to_extend A GRanges with exons to extend
#'                           (strand * are excluded)
#' @param input_gr_to_overlap A GRanges with intervals to overlap
#' @param verbose An integer that indicates the level of verbosity
#'                0 = silent, 1 = statistics, 2 = progression.
#' @return A GRanges which is a subset of `input_gr_to_extend` where
#'         3' end have been modified to match the 3' end of
#'         `input_gr_to_overlap` if they overlap
#'         (initial start and end have been stored into old_start and old_end)
#' @importFrom GenomicRanges strand
extend_using_overlap <- function(input_gr_to_extend, input_gr_to_overlap,
                                  verbose = 1) {
    # Remove strands which are not in + - and non exonic features
    input_gr_to_extend <- subset(
        input_gr_to_extend,
        type == "exon" &
            as.character(strand(input_gr_to_extend)) != "*"
    )
    # Remove non-exonic features from input_gr_to_overlap
    input_gr_to_overlap <- subset(input_gr_to_overlap, type == "exon")
    if (verbose > 1) {
        message("Compute overlap. Between ", length(input_gr_to_extend),
                " exons and ", length(input_gr_to_overlap), " exons.")
    }
    # Compute overlap:
    full_overlap <- suppressWarnings(GenomicRanges::findOverlaps(
        input_gr_to_extend,
        input_gr_to_overlap
    ))
    if (verbose > 1) {
        message("Split by overlap.")
    }
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
    # Store old coordinates
    input_gr_to_extend_intersect_selected$old_start <-
        GenomicRanges::start(input_gr_to_extend_intersect_selected)
    input_gr_to_extend_intersect_selected$old_end <-
        GenomicRanges::end(input_gr_to_extend_intersect_selected)
    # Extend
    input_gr_to_extend_intersect_selected_extended <-
        GenomicRanges::resize(input_gr_to_extend_intersect_selected,
               width =
                   GenomicRanges::width(input_gr_to_extend_intersect_selected) +
                   potential_extention_selected)
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
#' A function that from a GRanges with 'old_start' and 'old_end'
#' Change the starts and ends to prevent collisions larger than
#' with old coordinates
#' @param input_gr A GRanges with 2 metas: 'old_start' and 'old_end'
#' @importFrom dplyr %>%
#' @return A list with:
#'         - 'pot_issues': A dataframe with overlaps between `input_gr`
#'                         and itself while gene_ids are different
#'         - 'new_gr': A GRanges identical to `input_gr` except that start/end
#'                     have been adjusted to prevent collisions.
adjust_for_collision <- function(input_gr) {
    # Use an 'id':
    input_gr$id <- paste0(input_gr$exon_id, "_", input_gr$transcript_id)
    # Select entries which have been extended:
    extended_gr <- subset(input_gr, old_start > start | old_end < end)
    # Compute overlap between extended and full
    # and get those from different gene_ids:
    pot_issues <- overlap_different_genes(extended_gr, input_gr)

    if (nrow(pot_issues) > 0) {
        # Add annotations from Subject
        pot_issues$subject_start <-
            GenomicRanges::start(input_gr)[pot_issues$subjectHits]
        pot_issues$subject_end <-
            GenomicRanges::end(input_gr)[pot_issues$subjectHits]
        pot_issues$subject_id <- input_gr$id[pot_issues$subjectHits]
        # Add annotations from Query
        pot_issues$query_old_start <-
            extended_gr$old_start[pot_issues$queryHits]
        pot_issues$query_old_end <-
            extended_gr$old_end[pot_issues$queryHits]
        # Summarize from tail-to-head collisions
        pot_issues_summary <- pot_issues %>%
            dplyr::group_by(queryHits) %>%
            dplyr::summarise(
                min_subject_start_collision =
                    tryCatch(
                        min(subject_start[query_old_end < subject_start]),
                        warning = function(w) {
                            NA
                        }
                    ),
                max_subject_end_collision =
                    tryCatch(
                        max(subject_end[query_old_start > subject_end]),
                        warning = function(w) {
                            NA
                        }
                    )
            )
        # Add to the corresponding GRanges
        pot_issues_df <-
            as.data.frame(extended_gr[pot_issues_summary$queryHits])
        pot_issues_df$queryHits <- pot_issues_summary$queryHits
        pot_issues_df <- merge(pot_issues_df, pot_issues_summary)
        # add fixed_start_collision and fixed_end_collision
        # Which takes into account collision (tail-to-head) with another gene_id
        pot_issues_df <- pot_issues_df %>%
            dplyr::group_by(queryHits) %>%
            dplyr::mutate(
                fixed_start_collision = ifelse(
                    strand == "+" |  is.na(max_subject_end_collision),
                    start,
                    max_subject_end_collision + 1
                ),
                fixed_end_collision = ifelse(
                    strand == "+" & !is.na(min_subject_start_collision),
                    min_subject_start_collision - 1,
                    end
                )
            ) %>%
            as.data.frame
        # Apply these modifications to input_gr
        new_start_values <- GenomicRanges::start(input_gr)
        new_start_values[match(pot_issues_df$id, input_gr$id)] <-
            pot_issues_df$fixed_start_collision
        input_gr <- GenomicRanges::`start<-`(input_gr, value = new_start_values)
        new_end_values <- GenomicRanges::end(input_gr)
        new_end_values[match(pot_issues_df$id, input_gr$id)] <-
            pot_issues_df$fixed_end_collision
        input_gr <- GenomicRanges::`end<-`(input_gr, value = new_end_values)

        # We need to recompute overlap
        # We can restrict to exons with issues:
        new_extended_gr <- subset(input_gr, old_start > start | old_end < end)
        input_gr_to_check <- subset(
            input_gr,
            id %in% c(new_extended_gr$id, pot_issues$subject_id)
        )

        new_pot_issues <- overlap_different_genes(
            new_extended_gr,
            input_gr_to_check
        )

        if (nrow(new_pot_issues) > 0) {
            # Annotate new_pot_issues with query:
            new_extended_gr$queryHits <- seq_along(new_extended_gr)
            new_pot_issues <- merge(new_pot_issues,
                                    as.data.frame(new_extended_gr))

            new_pot_issues$extension_width <- with(
                new_pot_issues,
                old_start - start + end - old_end
            )
            new_pot_issues$query_extension_limit <- with(
                new_pot_issues,
                ifelse(
                    strand == "+",
                    end,
                    start
                )
            )
            # Annotate new_pot_issues with subject:
            new_pot_issues$subject_extension_limit <- ifelse(
                new_pot_issues$strand == "+",
                GenomicRanges::end(input_gr_to_check)[new_pot_issues$subjectHits],
                GenomicRanges::start(input_gr_to_check)[new_pot_issues$subjectHits]
            )
            new_pot_issues$subject_id <-
                input_gr_to_check$id[new_pot_issues$subjectHits]

            # We should have now only genes that previously overlapped:
            # And:
            # - only one has been extended:
            # old:
            # ---------------------->
            #     ---->
            # new:
            # ---------------------->
            #     --------->
            # or
            # old:
            # ---------------------->
            #     ---->
            # new:
            # ---------------------->
            #     ------------------>
            # or
            # old:
            # ---------/     \------------->
            #     ---->
            # new:
            # ---------/     \------------->
            #     --------->
            # or
            # old:
            # ---------/     \------------->
            #     -->
            # new:
            # ---------/     \------------->
            #     --------->
            # - are both extended up to the same base,
            #   only one should be extended:
            # Like
            # ------------------->
            # ------>
            # or
            #    ---------------->
            # ------>
            # or
            # ------------------->
            #   ---->

            # Set a variable for both cases:
            new_pot_issues$keep <- FALSE
            # Deal with case 1:
            new_pot_issues_case1 <- subset(
                new_pot_issues,
                !subject_id %in% new_pot_issues$id
            )
            if (nrow(new_pot_issues_case1) > 0) {
                new_pot_issues_case1$subject_start <-
                    GenomicRanges::start(input_gr_to_check)[new_pot_issues_case1$subjectHits]
                new_pot_issues_case1$subject_end <-
                    GenomicRanges::end(input_gr_to_check)[new_pot_issues_case1$subjectHits]

                # Cases where extended part overlap subject
                # should not be extended
                new_pot_issues_case1 <- new_pot_issues_case1 %>%
                    dplyr::mutate(
                        case = "case1",
                        keep = ifelse(
                            strand == "+",
                            old_end == subject_end,
                            old_start == subject_start
                        )
                    )
            }
            # Deal with case 2:
            new_pot_issues_case2 <- subset(
                new_pot_issues,
                subject_id %in% new_pot_issues$id
            )
            if (nrow(new_pot_issues_case2) > 0) {
                # Then only the gene_id
                # with the smallest extension_width will be extended.
                # The others will be put back to before extension
                new_pot_issues_case2 <- new_pot_issues_case2 %>%
                    dplyr::group_by(strand, query_extension_limit) %>%
                    dplyr::mutate(
                        case = "case2",
                        min_extension = min(extension_width),
                        keep = query_gene_id %in%
                            query_gene_id[extension_width == min_extension]
                    ) %>%
                    dplyr::ungroup()
            }

            new_pot_issues_summary <- rbind(
                new_pot_issues_case1 %>%
                    dplyr::select(id, strand, keep,
                                  start, old_start, end, old_end),
                new_pot_issues_case2 %>%
                    dplyr::select(id, strand, keep,
                                  start, old_start, end, old_end)
            ) %>%
                dplyr::mutate(
                    fixed_start_same_extension = ifelse(
                        strand == "+" | keep,
                        start,
                        old_start
                    ),
                    fixed_end_same_extension = ifelse(
                        strand == "+" & !keep,
                        old_end,
                        end
                    )
                ) %>%
                dplyr::group_by(id) %>%
                dplyr::summarise(
                    max_fixed_start_same_extension = max(fixed_start_same_extension),
                    min_fixed_end_same_extension = min(fixed_end_same_extension)
                ) %>%
                as.data.frame()

            # Apply these modifications to input_gr
            new_start_values <- GenomicRanges::start(input_gr)
            new_start_values[match(new_pot_issues_summary$id, input_gr$id)] <-
                new_pot_issues_summary$max_fixed_start_same_extension
            input_gr <-
                GenomicRanges::`start<-`(input_gr, value = new_start_values)
            new_end_values <- GenomicRanges::end(input_gr)
            new_end_values[match(new_pot_issues_summary$id, input_gr$id)] <-
                new_pot_issues_summary$min_fixed_end_same_extension
            input_gr <- GenomicRanges::`end<-`(input_gr, value = new_end_values)

            # Combine both pot_issues
            pot_issues_df <- merge(pot_issues_df, new_pot_issues_summary)
        }

    } else {
        pot_issues_df <- data.frame()
    }
    # Remove the id:
    input_gr$id <- NULL
    return(list(pot_issues = pot_issues_df,
                new_gr = input_gr))
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
#' @param verbose An integer that indicates the level of verbosity
#'                0 = silent, 1 = statistics, 2 = progression.
#' @importFrom dplyr %>%
#' @importFrom GenomicRanges strand
#' @return A GRanges identical to `input_gr_to_extend` with new exons
#'         whose `exon_id` contains BREW3R. `exon_number` may have changed.
#' @details
#' Potential new exons will be filtered for collision with exons present
#' in the first GRanges even if they belong to the same gene_id.
#' For the moment all potential exons extensions are added to the
#' same existing transcript_id so introns maybe artificial introns.
#'
add_new_exons <- function(input_gr_to_extend, input_gr_with_new_exons,
                          verbose = 1) {
    # Subset input_gr_with_new_exons
    input_gr_with_new_exons <- subset(
        input_gr_with_new_exons,
        strand(input_gr_with_new_exons) != "*"
    )
    # Add an id:
    input_gr_with_new_exons$id <- seq_along(input_gr_with_new_exons)
    input_gr_to_extend$id <- paste0(input_gr_to_extend$exon_id, "_",
                                    input_gr_to_extend$transcript_id)
    if (verbose > 1) {
        message("Extract last exon of first GRanges.")
    }
    # Get the last exons of input_gr_to_extend
    input_gr_last <- extract_last_exon(input_gr_to_extend)
    # First find exons which 'ends' at the same position
    if (verbose > 1) {
        message("Compute overlap with second GRanges.")
    }
    ov <- suppressWarnings(GenomicRanges::findOverlaps(input_gr_last, input_gr_with_new_exons))
    ov_df <- as.data.frame(ov)
    if (verbose > 1) {
        message("Only keep those with identical three prime.")
    }
    # Get extremities
    ov_df$query_extremity <- ifelse(
        strand(input_gr_last[ov_df$queryHits]) == "+",
        GenomicRanges::end(input_gr_last[ov_df$queryHits]),
        GenomicRanges::start(input_gr_last[ov_df$queryHits])
    )
    ov_df$subject_extremity <- ifelse(
        strand(input_gr_with_new_exons[ov_df$subjectHits]) == "+",
        GenomicRanges::end(input_gr_with_new_exons[ov_df$subjectHits]),
        GenomicRanges::start(input_gr_with_new_exons[ov_df$subjectHits])
    )
    ov_df_potential <- subset(
        ov_df,
        query_extremity == subject_extremity
    )
    if (nrow(ov_df_potential) == 0) {
        input_gr_to_extend$id <- NULL
        return(input_gr_to_extend)
    }
    if (verbose > 1) {
        message("Only keep those with exons to be added.")
    }
    # We need to check if it does not matches transcript subject extremity
    # First get transcript_id from exon overlapped
    ov_df_potential$subject_transcript_id <-
        input_gr_with_new_exons$transcript_id[ov_df_potential$subjectHits]

    # Get transcript_extremity and transcript_id from input_gr_with_new_exons
    input_gr_with_new_exons_subset <- subset(
        input_gr_with_new_exons,
        transcript_id %in% ov_df_potential$subject_transcript_id
    )
    input_gr_with_new_exons_subset_range <-
        unlist(range(
            GenomicRanges::split(input_gr_with_new_exons_subset,
                                 input_gr_with_new_exons_subset$transcript_id)
        ))
    # Here seems slow
    subject_transcript_extremity_t <-
        GenomicRanges::as.data.frame(input_gr_with_new_exons_subset_range)
    rownames(subject_transcript_extremity_t) <-
        names(input_gr_with_new_exons_subset_range)
    subject_transcript_extremity_t$subject_transcript_extremity <-
        with(subject_transcript_extremity_t,
             ifelse(strand == "+",
                    end,
                    start)
        )
    ov_df_potential$subject_transcript_extremity <-
        subject_transcript_extremity_t[ov_df_potential$subject_transcript_id,
                                       "subject_transcript_extremity"]

    # Remove overlaps where it is the transcript subject extremity
    if (verbose > 1) {
        message("Remove overlaps where transcript 3' matches between both")
    }
    ov_df_potential <- subset(
        ov_df_potential,
        subject_extremity != subject_transcript_extremity
    )

    if (nrow(ov_df_potential) == 0) {
        input_gr_to_extend$id <- NULL
        return(input_gr_to_extend)
    }
    if (verbose > 1) {
        message("Attribute exons to be added to the good transcript (this is slow).")
    }

    # Get all exons from subject_transcript_id
    # which are after the query_extremity
    # First get all exons of subject_transcript_id
    # (once per queryHits)
    # This is slow could be improved:
    all_exons_interesting <- merge(
        data.frame(
            queryHits = ov_df_potential$queryHits,
            transcript_id = ov_df_potential$subject_transcript_id
        ),
        as.data.frame(subset(
            input_gr_with_new_exons,
            transcript_id %in% ov_df_potential$subject_transcript_id
        ))
    )
    # Annotate with query start/end/transcript_id
    all_exons_interesting$query_start <-
        GenomicRanges::start(input_gr_last[all_exons_interesting$queryHits])
    all_exons_interesting$query_end <-
        GenomicRanges::end(input_gr_last[all_exons_interesting$queryHits])
    all_exons_interesting$new_transcript_id <-
        input_gr_last$transcript_id[all_exons_interesting$queryHits]
    # Filter
    all_exons_interesting <- subset(
        all_exons_interesting,
        ifelse(strand == "+",
               start > query_end,
               end < query_start)
    )
    if (verbose > 0) {
        message("Found ", length(unique(all_exons_interesting$id)), " exons",
                " that may be included into ", length(unique(all_exons_interesting$queryHits)),
                " transcripts.")
        if (verbose > 1) {
            message("Check added exons which would ",
                    "overlap existing annotations.")
        }
    }


    # Check that the new exons do not overlap existing annotations
    # ov_new_exons is the overlap between:
    # queries = granges corresponding to all_exons_interesting
    # subjects = input_gr_to_extend
    # What could be improved is to give as gene_id the future gene_id
    # And select overlaps only for different gene_id
    ov_new_exons <- as.data.frame(
        suppressWarnings(GenomicRanges::findOverlaps(
            input_gr_with_new_exons[match(all_exons_interesting$id,
                                          input_gr_with_new_exons$id)],
            input_gr_to_extend
        ))
    )

    input_gr_to_extend_intersected <-
        input_gr_to_extend[ov_new_exons$subjectHits]

    ov_new_exons$transcript_id <-
        all_exons_interesting$transcript_id[ov_new_exons$queryHits]
    ov_new_exons$transcript_id_to_extend <-
        all_exons_interesting$new_transcript_id[ov_new_exons$queryHits]
    ov_new_exons$strand <- all_exons_interesting$strand[ov_new_exons$queryHits]
    ov_new_exons$to_extend_start <-
        GenomicRanges::start(input_gr_to_extend_intersected)
    ov_new_exons$to_extend_end <-
        GenomicRanges::end(input_gr_to_extend_intersected)
    ov_new_exons$to_extend_strand <-
        as.character(
            strand(input_gr_to_extend_intersected)
        )

    # Get the 5' end of the exon (from input_gr_to_extend)
    # that overlapped with the exons we wanted to add
    ov_new_exons$to_extend_5p <- with(ov_new_exons,
                                      ifelse(to_extend_strand == "+",
                                             to_extend_start,
                                             to_extend_end)
    )

    # Get of each transcript_id that wanted to be used for extension
    # And for each transcript_id to be extended
    # The first colliding base
    ov_new_exons_summary <- ov_new_exons %>%
        dplyr::mutate(transcript_pair =
                          paste0(transcript_id, "_", transcript_id_to_extend)
        ) %>%
        dplyr::group_by(transcript_pair) %>%
        dplyr::summarise(first_colliding_base = ifelse(strand[1] == "+",
                                                       min(to_extend_5p),
                                                       max(to_extend_5p))) %>%
        as.data.frame()
    rownames(ov_new_exons_summary) <-
        ov_new_exons_summary$transcript_pair
    # Add all info to all_exons_interesting
    all_exons_interesting$transcript_pair <-
        paste0(all_exons_interesting$transcript_id,
               "_", all_exons_interesting$new_transcript_id)
    # This is slow I don't know why:
    all_exons_interesting$first_colliding_base <-
        ov_new_exons_summary[
            all_exons_interesting$transcript_pair,
            "first_colliding_base"]
    if (verbose > 1) {
        message("Remove new candidate exons",
                " three prime of overlap with existing annotations.")
    }

    # Filter exons after collision:
    all_exons_interesting <-
        subset(all_exons_interesting,
               is.na(first_colliding_base) |
                   (strand == "+" & start < first_colliding_base) |
                   (strand == "-" & end > first_colliding_base)
        )
    if (verbose > 0) {
        message("Stay ", length(unique(all_exons_interesting$id)), " candidate",
                " exons that may be included into ",
                length(unique(all_exons_interesting$queryHits)),
                " transcripts.")
        if (verbose > 1) {
            message("Adjust start/end of candidate exons to not ",
                    "overlap with existing annotations.")
        }
    }

    if (nrow(all_exons_interesting) == 0) {
        input_gr_to_extend$id <- NULL
        return(input_gr_to_extend)
    }
    # Adjust start/end to avoid collision:
    all_exons_interesting$start <-
        with(all_exons_interesting,
             ifelse(strand == "+" | is.na(first_colliding_base) | start > first_colliding_base,
                    start,
                    first_colliding_base + 1)
        )
    all_exons_interesting$end <-
        with(all_exons_interesting,
             ifelse(strand == "-" | is.na(first_colliding_base) | end < first_colliding_base,
                    end,
                    first_colliding_base - 1)
        )

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

    if (verbose > 1) {
        message("Combine exons extending the same transcript.")
    }
    all_exons_to_add_gr <- unlist(
        GenomicRanges::reduce(
            GenomicRanges::split(
                GenomicRanges::makeGRangesFromDataFrame(
                    all_exons_interesting
                ),
                all_exons_interesting$new_transcript_id
            )
        )
    )
    all_exons_to_add_gr$transcript_id <- names(all_exons_to_add_gr)
    names(all_exons_to_add_gr) <- NULL

    if (verbose > 0) {
        message("Finally ", length(all_exons_to_add_gr), " combined exons",
                " will be included into ",
                length(unique(all_exons_to_add_gr$transcript_id)),
                " transcripts.")
        if (verbose > 1) {
            message("Annotate new exons with transcript info.")
        }
    }
    # Add transcript_id infos:
    transcript_annotation <-
        as.data.frame(GenomicRanges::mcols(input_gr_last[unique(all_exons_interesting$queryHits)])) %>%
        dplyr::select(tidyselect::starts_with("gene_") |
                          tidyselect::starts_with("transcript_") |
                          tidyselect::matches("source")) %>%
        unique() %>%
        as.data.frame()

    GenomicRanges::mcols(all_exons_to_add_gr) <-
        transcript_annotation[match(all_exons_to_add_gr$transcript_id,
                                    transcript_annotation$transcript_id), ]
    all_exons_to_add_gr$type <- "exon"

    if (verbose > 1) {
        message("Annotate new exons with exon_id.")
    }
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


    if (verbose > 1) {
        message("Update exon_number.")
    }
    # Update exon_number
    # There are 2 exon_number modes:
    # 1. from 5p to 3p (I call it transcript)
    # 2. from start to end whatever the strand is (I call it coordinate)
    # To get the mode I select a second exon with strand -
    if (!"exon_number" %in% colnames(GenomicRanges::mcols(input_gr_to_extend))) {
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
    exon_numbers <- unlist(exon_numbers_list[unique(gr_to_annotate$transcript_id)])
    gr_to_annotate$exon_number <- unname(exon_numbers)
    new_gr <- c(gr_to_annotate, gr_to_add)
    # Remove the id:
    new_gr$id <- NULL

    return(new_gr)
}
