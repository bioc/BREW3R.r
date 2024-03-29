
input_to_overlap_case1_2_3_4_6_7_8 <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 3,
                              end = 25),
    strand = "+",
    gene_id = "geneA",
    transcript_id = "transcriptA",
    type = "exon",
    exon_id = "exonA"
)
input_to_overlap_case5_9 <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(1, 33, 45, 72),
                              end = c(25, 40, 60, 75)),
    strand = "+",
    gene_id = "geneA",
    transcript_id = "transcriptA",
    type = "exon",
    exon_id = c("exonA", "exonB", "exonC", "exonD")
)


# Case1:
# input_gene1:      ------->
# input_gene2:                    -------->
# to_overlap:  ---------------------->

test_that("case1 works", {
    expect_equal(
        sort(extend_granges(
            GenomicRanges::GRanges(
                seqnames = "chr1",
                ranges = IRanges::IRanges(start = c(5, 20),
                                          end = c(10, 30)),
                strand = "+",
                gene_id = c("gene1", "gene2"),
                transcript_id = c("transcript1", "transcript2"),
                type = "exon",
                exon_id = c("exon1", "exon2")
            ),
            input_to_overlap_case1_2_3_4_6_7_8,
            verbose = 0
        )),
        GenomicRanges::GRanges(
            seqnames = "chr1",
            ranges = IRanges::IRanges(start = c(5, 20),
                                      end = c(19, 30)),
            strand = "+",
            gene_id = c("gene1", "gene2"),
            transcript_id = c("transcript1", "transcript2"),
            type = "exon",
            exon_id = c("exon1.ext", "exon2")
        )
    )
})

# Case 2:
# input_gene1:      ------->
# input_gene2:                    -------->
# to_overlap:           ------------------>

test_that("case2 works", {
    expect_equal(
        sort(extend_granges(
            GenomicRanges::GRanges(
                seqnames = "chr1",
                ranges = IRanges::IRanges(start = c(5, 20),
                                          end = c(10, 25)),
                strand = "+",
                gene_id = c("gene1", "gene2"),
                transcript_id = c("transcript1", "transcript2"),
                type = "exon",
                exon_id = c("exon1", "exon2")
            ),
            input_to_overlap_case1_2_3_4_6_7_8,
            verbose = 0
        )),
        GenomicRanges::GRanges(
            seqnames = "chr1",
            ranges = IRanges::IRanges(start = c(5, 20),
                                      end = c(19, 25)),
            strand = "+",
            gene_id = c("gene1", "gene2"),
            transcript_id = c("transcript1", "transcript2"),
            type = "exon",
            exon_id = c("exon1.ext", "exon2")
        )
    )
})

# Case 3:
# input_gene1:      ------->
# input_gene2:                    ---->
# to_overlap:       ---------------------->

test_that("case3 works", {
    expect_equal(
        sort(extend_granges(
            GenomicRanges::GRanges(
                seqnames = "chr1",
                ranges = IRanges::IRanges(start = c(5, 20),
                                          end = c(10, 22)),
                strand = "+",
                gene_id = c("gene1", "gene2"),
                transcript_id = c("transcript1", "transcript2"),
                type = "exon",
                exon_id = c("exon1", "exon2")
            ),
            input_to_overlap_case1_2_3_4_6_7_8,
            verbose = 0
        )),
        GenomicRanges::GRanges(
            seqnames = "chr1",
            ranges = IRanges::IRanges(start = c(5, 20),
                                      end = c(19, 25)),
            strand = "+",
            gene_id = c("gene1", "gene2"),
            transcript_id = c("transcript1", "transcript2"),
            type = "exon",
            exon_id = c("exon1.ext", "exon2.ext")
        )
    )
})

# Case 4:
# input_gene1:      ------->
# input_gene2:      ------------->
# to_overlap:       ---------------------->

test_that("case4 works", {
    expect_equal(
        sort(extend_granges(
            GenomicRanges::GRanges(
                seqnames = "chr1",
                ranges = IRanges::IRanges(start = c(5, 5),
                                          end = c(10, 22)),
                strand = "+",
                gene_id = c("gene1", "gene2"),
                transcript_id = c("transcript1", "transcript2"),
                type = "exon",
                exon_id = c("exon1", "exon2")
            ),
            input_to_overlap_case1_2_3_4_6_7_8,
            verbose = 0
        )),
        GenomicRanges::GRanges(
            seqnames = "chr1",
            ranges = IRanges::IRanges(start = c(5, 5),
                                      end = c(10, 25)),
            strand = "+",
            gene_id = c("gene1", "gene2"),
            transcript_id = c("transcript1", "transcript2"),
            type = "exon",
            exon_id = c("exon1", "exon2.ext")
        )
    )
})

# Case 5:
# input_gene1:      ------->
# input_gene2:                    -------->
# to_overlap:           ----------------------/\----/\-/\->

test_that("case5 works", {
    expect_equal(
        sort(extend_granges(
            GenomicRanges::GRanges(
                seqnames = "chr1",
                ranges = IRanges::IRanges(start = c(5, 20),
                                          end = c(10, 22)),
                strand = "+",
                gene_id = c("gene1", "gene2"),
                transcript_id = c("transcript1", "transcript2"),
                type = "exon",
                exon_id = c("exon1", "exon2")
            ),
            input_to_overlap_case5_9,
            verbose = 0
        )),
        GenomicRanges::GRanges(
            seqnames = "chr1",
            ranges = IRanges::IRanges(start = c(5, 20, 33, 45, 72),
                                      end = c(19, 25, 40, 60, 75)),
            strand = "+",
            gene_id = rep(c("gene1", "gene2"), c(1, 4)),
            transcript_id = rep(c("transcript1", "transcript2"), c(1, 4)),
            type = "exon",
            exon_id = c("exon1.ext", "exon2.ext", "BREW3R0000000001",
                        "BREW3R0000000002", "BREW3R0000000003"),
            exon_number = c(1, 1:4)
        )
    )
})

# Case 6:
# input_transcript1:     --------->
# input_transcript2:     ----------/        \-------->
# to_overlap:               ----------->

test_that("case6 works", {
    expect_equal(
        sort(extend_granges(
            GenomicRanges::GRanges(
                seqnames = "chr1",
                ranges = IRanges::IRanges(start = c(1, 1, 30),
                                          end = c(10, 10, 40)),
                strand = "+",
                gene_id = "gene1",
                transcript_id = c("transcript1", "transcript2", "transcript2"),
                type = "exon",
                exon_id = c("exon1", "exon1", "exon2")
            ),
            input_to_overlap_case1_2_3_4_6_7_8,
            verbose = 0
        )),
        sort(GenomicRanges::GRanges(
            seqnames = "chr1",
            ranges = IRanges::IRanges(start = c(1, 1, 30),
                                      end = c(25, 10, 40)),
            strand = "+",
            gene_id = "gene1",
            transcript_id = c("transcript1", "transcript2", "transcript2"),
            type = "exon",
            exon_id = c("exon1.ext", "exon1", "exon2")
        ))
    )
})

# Case 6bis:
# input_gene1:     --------->
# input_gene2:     ----------/        \-------->
# to_overlap:          ----------->

test_that("case6bis works", {
    expect_equal(
        sort(extend_granges(
            GenomicRanges::GRanges(
                seqnames = "chr1",
                ranges = IRanges::IRanges(start = c(1, 1, 30),
                                          end = c(10, 10, 40)),
                strand = "+",
                gene_id = c("gene1", "gene2", "gene2"),
                transcript_id = c("transcript1", "transcript2", "transcript2"),
                type = "exon",
                exon_id = c("exon1", "exon2", "exon3")
            ),
            input_to_overlap_case1_2_3_4_6_7_8,
            verbose = 0
        )),
        sort(GenomicRanges::GRanges(
            seqnames = "chr1",
            ranges = IRanges::IRanges(start = c(1, 1, 30),
                                      end = c(25, 10, 40)),
            strand = "+",
            gene_id = c("gene1", "gene2", "gene2"),
            transcript_id = c("transcript1", "transcript2", "transcript2"),
            type = "exon",
            exon_id = c("exon1.ext", "exon2", "exon3")
        ))
    )
})
# Case 7:
# input_transcript1:     ------>
# input_transcript2:     ----------/        \-------->
# to_overlap:               ----------->

test_that("case7 works", {
    expect_equal(
        sort(extend_granges(
            GenomicRanges::GRanges(
                seqnames = "chr1",
                ranges = IRanges::IRanges(start = c(1, 1, 30),
                                          end = c(8, 10, 40)),
                strand = "+",
                gene_id = "gene1",
                transcript_id = c("transcript1", "transcript2", "transcript2"),
                type = "exon",
                exon_id = c("exon1", "exon2", "exon3")
            ),
            input_to_overlap_case1_2_3_4_6_7_8,
            verbose = 0
        )),
        sort(GenomicRanges::GRanges(
            seqnames = "chr1",
            ranges = IRanges::IRanges(start = c(1, 1, 30),
                                      end = c(25, 10, 40)),
            strand = "+",
            gene_id = "gene1",
            transcript_id = c("transcript1", "transcript2", "transcript2"),
            type = "exon",
            exon_id = c("exon1.ext", "exon2", "exon3")
        ))
    )
})

# Case 8 should not be extended:
# input_gene1:        ----->
# input_gene2:     ----------------/      \--------
# to_overlap:      ------------------>

test_that("case8 works", {
    expect_equal(
        sort(extend_granges(
            GenomicRanges::GRanges(
                seqnames = "chr1",
                ranges = IRanges::IRanges(start = c(1, 1, 30),
                                          end = c(8, 10, 40)),
                strand = "+",
                gene_id = c("gene1", "gene2", "gene2"),
                transcript_id = c("transcript1", "transcript2", "transcript2"),
                type = "exon",
                exon_id = c("exon1", "exon2", "exon3")
            ),
            input_to_overlap_case1_2_3_4_6_7_8,
            verbose = 0
        )),
        sort(GenomicRanges::GRanges(
            seqnames = "chr1",
            ranges = IRanges::IRanges(start = c(1, 1, 30),
                                      end = c(8, 10, 40)),
            strand = "+",
            gene_id = c("gene1", "gene2", "gene2"),
            transcript_id = c("transcript1", "transcript2", "transcript2"),
            type = "exon",
            exon_id = c("exon1", "exon2", "exon3")
        ))
    )
})

# Case 9:
# input_gene1:  --->
# input_gene2:                   --->
# to_overlap:   ----/\--------/\---/  \---->

test_that("case9 works with debug and verbose and overlapresolution", {
    overlap_resolution_fn <- tempfile()
    expect_equal(
        sort(extend_granges(
            GenomicRanges::GRanges(
                seqnames = "chr1",
                ranges = IRanges::IRanges(start = c(5, 55),
                                          end = c(10, 70)),
                strand = "+",
                gene_id = c("gene1", "gene2"),
                transcript_id = c("transcript1", "transcript2"),
                type = "exon",
                exon_id = c("exon1", "exon2")
            ),
            input_to_overlap_case5_9,
            verbose = 2, debug = TRUE,
            overlap_resolution_fn = overlap_resolution_fn
        )),
        GenomicRanges::GRanges(
            seqnames = "chr1",
            ranges = IRanges::IRanges(start = c(5, 33, 45, 55),
                                      end = c(25, 40, 54, 70)),
            strand = "+",
            gene_id = rep(c("gene1", "gene2"), c(3, 1)),
            transcript_id = rep(c("transcript1", "transcript2"), c(3, 1)),
            type = "exon",
            exon_id = c("exon1.ext", "BREW3R0000000001",
                        "BREW3R0000000002", "exon2"),
            exon_number = c(1:3, 1)
        )
    )
})

