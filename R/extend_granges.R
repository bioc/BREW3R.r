#' Extend GRanges
#'
#' A function that from a GRanges from gtf
#' will extend the 3' of transcripts using
#' another GRanges from gtf as a template
#' @param input_gr_to_extend A GRanges to extend
#'                           (only exons are kept and strand * are excluded)
#' @param input_gr_to_overlap A GRanges with intervals to overlap
#' @param extend_existing_exons A boolean that indicates
#'                              if existing exons should be extended
#' @param add_new_exons A boolean that indicates
#'                      if new exons with compatible splicing event
#'                      should be added
#' @param overlap_resolution_fn A file path where the dataframe giving details
#'                              on the collision resolution is written
#' @return A GRanges based on `input_gr_to_extend` where exons are extended
#'         and new exons can be added.
#'         Exons extended will have a '.ext' suffix to the original exon_id.
#'         Exons added will have a exon_id starting with 'BREW3R'.
#' @details
#' During the extension process a special care is taking to prevent extension
#' which would lead to overlap between different gene_ids.
#' @importFrom GenomicRanges strand
#' @export
#' @examples
#' # Very simple case
#' # input_gr:      ------->           ----->
#' # to_overlap:  ------------------>
#'
#' # output:        ---------------->  ----->
#'
#' input_gr <- GenomicRanges::GRanges(
#'     seqnames = "chr1",
#'     ranges = IRanges::IRanges(start = c(5, 20),
#'                               end = c(10, 30)),
#'     strand = "+",
#'     gene_id = c("gene1", "gene2"),
#'     transcript_id = c("transcript1", "transcript2"),
#'     type = "exon",
#'     exon_id = c("exon1", "exon2")
#' )
#'
#' input_gr_to_overlap <- GenomicRanges::GRanges(
#'     seqnames = "chr1",
#'     ranges = IRanges::IRanges(start = 3,
#'                               end = 15),
#'     strand = "+",
#'     gene_id = "geneA",
#'     transcript_id = "transcriptA",
#'     type = "exon",
#'     exon_id = "exonA"
#' )
#'
#' extend_granges(input_gr, input_gr_to_overlap)
#'
extend_granges <- function(input_gr_to_extend, input_gr_to_overlap,
                           extend_existing_exons = TRUE, add_new_exons = TRUE,
                           overlap_resolution_fn = NULL) {
    # Apply the filters
    input_gr_to_extend <- subset(
        input_gr_to_extend,
        type == "exon" &
            as.character(strand(input_gr_to_extend)) != "*"
    )

    if (extend_existing_exons) {
        progression_msg("Extend last exons.\nGetting last exons.\n")
        # First get the last exons
        last_exons_gr <- extract_last_exons(input_gr_to_extend)
        rlang::inform(
            paste("Found", length(last_exons_gr),
                  "last exons to potentially extend.\n")
        )
        debug_msg("last_exons_gr")
        # Then, these exons are extended ignoring potential collisions
        last_exons_gr_extended <-
            extend_using_overlap(last_exons_gr, input_gr_to_overlap)
        rlang::inform(
            paste(length(last_exons_gr_extended), " exons may be extended.\n")
        )
        progression_msg("Checking for collision with other genes.\n")
        debug_msg("last_exons_gr_extended")
        # Prepare the non-extended exons:
        last_exons_gr_extended$id <-
            paste0(last_exons_gr_extended$transcript_id, "_",
                   last_exons_gr_extended$exon_id)
        input_gr_to_extend$id <- paste0(input_gr_to_extend$transcript_id, "_",
                                        input_gr_to_extend$exon_id)
        non_last_exons_extended_gr <- subset(
            input_gr_to_extend,
            !id %in% last_exons_gr_extended$id
        )
        # Remove the id:
        last_exons_gr_extended$id <- NULL
        non_last_exons_extended_gr$id <- NULL
        # Add needed metadata:
        non_last_exons_extended_gr$old_width <-
            GenomicRanges::width(non_last_exons_extended_gr)
        debug_msg("non_last_exons_extended_gr")
        # Then overlaps are resolved and exons extensions are removed
        # or shortened:
        extension_resolved <- adjust_for_collision(
            c(
                non_last_exons_extended_gr,
                last_exons_gr_extended
            )
        )
        if (!is.null(overlap_resolution_fn)) {
            tryCatch(
                utils::write.table(extension_resolved[["pot_issues"]],
                                   overlap_resolution_fn,
                                   quote = FALSE, sep = "\t",
                                   row.names = FALSE),
                error = function(e){
                    message("Could not save table.\n Continuing merge.\n")
                }
            )
        }
        debug_msg("extension_resolved$pot_issues")
        extension_resolved_gr <- extension_resolved[["new_gr"]]
        debug_msg("extension_resolved_gr")
        # exon_id of exons changed by BREW3R are modified
        modified <-
            GenomicRanges::width(extension_resolved_gr) !=
            extension_resolved_gr$old_width
        extension_resolved_gr$exon_id[modified] <- paste0(
            extension_resolved_gr$exon_id[modified],
            ".ext"
        )
        rlang::inform(
            paste(sum(modified), "exons have been extended",
                  "while preventing collision with other genes.\n")
        )
        # We remove old_width
        extension_resolved_gr$old_width <- NULL
        debug_msg("extension_resolved_gr")
    } else {
        extension_resolved_gr <- input_gr_to_extend
    }
    if (add_new_exons) {
        progression_msg("Adding exons after existing ones.\n")
        extension_resolved_gr_new_exons <-
            add_new_exons(extension_resolved_gr,
                          input_gr_to_overlap)
        progression_msg(
            paste("Added",
                  length(grep("BREW3R",
                              extension_resolved_gr_new_exons$exon_id)) -
                      length(grep("BREW3R",
                                  extension_resolved_gr$exon_id)),
                  " exons.\n")
        )
    } else {
        extension_resolved_gr_new_exons <- extension_resolved_gr
    }
    return(extension_resolved_gr_new_exons)
}
