# Define GR Used in different places
input_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(1, 10, 20, 30),
                              end = c(3, 13, 23, 33)),
    strand = rep(c("+", "-"), each = 2),
    gene_id = rep(c("gene1", "gene2"), each = 2),
    transcript_id = c("transcript1", "transcript2",
                      "transcript3", "transcript3"),
    type = "exon"
)

input_gr_to_overlap <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(1, 7, 17, 25, 28),
                              end = c(3, 15, 23, 30, 35)),
    strand = c("+", "+", "+", "-", "-"),
    gene_id = c(rep(c("geneA", "geneB"), each = 2),  "geneC"),
    transcript_id = c("transcriptA", "transcriptB", "transcriptC",
                      "transcriptC", "transcriptD"),
    type = "exon"
)

small_gr_pot_extended <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(10, 25),
                              end = c(15, 33)),
    strand = c("+", "-"),
    gene_id = c("gene1", "gene2"),
    transcript_id = c("transcript2", "transcript3"),
    type = "exon",
    old_width = c(4, 4)
)
# Case1:
# input_gene1:      ------->
# input_gene2:                    -------->
# to_overlap:  ---------------------->
# or Case 2:
# to_overlap:           ------------------>
# Case 3:
# input_gene1:      ------->                (will be extended by Case1)
# input_gene2:                    ---->
# to_overlap:       ---------------------->
# Case 4:
# input_gene1:      ------->
# input_gene2:      ------------->
# to_overlap:       ---------------------->
# Case 4bis:
# input_gene1:      ------->
# input_gene2:      ------------->
# to_overlap:       ------------->
# Case 5:
# input_gene1:      ------->
# input_gene2:                    -------->
# to_overlap:           ----------------------/\----/\-/\->
# Case 6:
# input_transcript1:     --------->
# input_transcript2:     ----------/        \-------->
# to_overlap:               ----------->
# Case 6bis:
# input_gene1:     --------->
# input_gene2:     ----------/        \-------->
# to_overlap:          ----------->
# Case 7:
# input_transcript1:     ------>
# input_transcript2:     ----------/        \-------->
# to_overlap:               ----------->
# Case 8 should not be extended:
# input_gene1:        ----->
# input_gene2:     ----------------/      \--------
# to_overlap:      ------------------>
# Case 9:
# input_gene1:  --->
# input_gene2:                   --->
# to_overlap:   ----/\--------/\---/  \---->


input_gr_case1 <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(5, 20),
                              end = c(10, 30)),
    strand = "+",
    gene_id = c("gene1", "gene2"),
    transcript_id = c("transcript1", "transcript2"),
    type = "exon",
    exon_id = c("exon1", "exon2")
)
input_gr_case2 <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(5, 20),
                              end = c(10, 25)),
    strand = "+",
    gene_id = c("gene1", "gene2"),
    transcript_id = c("transcript1", "transcript2"),
    type = "exon",
    exon_id = c("exon1", "exon2")
)
input_gr_case3_5 <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(5, 20),
                              end = c(10, 22)),
    strand = "+",
    gene_id = c("gene1", "gene2"),
    transcript_id = c("transcript1", "transcript2"),
    type = "exon",
    exon_id = c("exon1", "exon2")
)
input_gr_case4 <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(5, 5),
                              end = c(10, 22)),
    strand = "+",
    gene_id = c("gene1", "gene2"),
    transcript_id = c("transcript1", "transcript2"),
    type = "exon",
    exon_id = c("exon1", "exon2")
)
input_gr_case4bis <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(5, 5),
                              end = c(10, 25)),
    strand = "+",
    gene_id = c("gene1", "gene2"),
    transcript_id = c("transcript1", "transcript2"),
    type = "exon",
    exon_id = c("exon1", "exon2")
)
input_gr_case6 <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(1, 1, 30),
                              end = c(10, 10, 40)),
    strand = "+",
    gene_id = "gene1",
    transcript_id = c("transcript1", "transcript2", "transcript2"),
    type = "exon",
    exon_id = c("exon1", "exon1", "exon2")
)
input_gr_case6bis <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(1, 1, 30),
                              end = c(10, 10, 40)),
    strand = "+",
    gene_id = c("gene1", "gene2", "gene2"),
    transcript_id = c("transcript1", "transcript2", "transcript2"),
    type = "exon",
    exon_id = c("exon1", "exon2", "exon3")
)
input_gr_case7 <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(1, 1, 30),
                              end = c(8, 10, 40)),
    strand = "+",
    gene_id = "gene1",
    transcript_id = c("transcript1", "transcript2", "transcript2"),
    type = "exon",
    exon_id = c("exon1", "exon2", "exon3")
)
input_gr_case8 <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(1, 1, 30),
                              end = c(8, 10, 40)),
    strand = "+",
    gene_id = c("gene1", "gene2", "gene2"),
    transcript_id = c("transcript1", "transcript2", "transcript2"),
    type = "exon",
    exon_id = c("exon1", "exon2", "exon3")
)
input_gr_case9 <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(5, 55),
                              end = c(10, 70)),
    strand = "+",
    gene_id = c("gene1", "gene2"),
    transcript_id = c("transcript1", "transcript2"),
    type = "exon",
    exon_id = c("exon1", "exon2")
)
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

input_gr_case1_pot_extended <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(5, 20),
                              end = c(25, 30)),
    strand = "+",
    gene_id = c("gene1", "gene2"),
    transcript_id = c("transcript1", "transcript2"),
    type = "exon",
    exon_id = c("exon1", "exon2"),
    old_width = c(6, 11)
)

input_gr_case1_extended_fixed <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(5, 20),
                              end = c(19, 30)),
    strand = "+",
    gene_id = c("gene1", "gene2"),
    transcript_id = c("transcript1", "transcript2"),
    type = "exon",
    exon_id = c("exon1", "exon2"),
    old_width = c(6, 11)
)

input_gr_case2_pot_extended <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(5, 20),
                              end = c(25, 25)),
    strand = "+",
    gene_id = c("gene1", "gene2"),
    transcript_id = c("transcript1", "transcript2"),
    type = "exon",
    exon_id = c("exon1", "exon2"),
    old_width = c(6, 6)
)

input_gr_case2_extended_fixed <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(5, 20),
                              end = c(19, 25)),
    strand = "+",
    gene_id = c("gene1", "gene2"),
    transcript_id = c("transcript1", "transcript2"),
    type = "exon",
    exon_id = c("exon1", "exon2"),
    old_width = c(6, 6)
)

input_gr_case3_pot_extended <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(5, 20),
                              end = c(25, 25)),
    strand = "+",
    gene_id = c("gene1", "gene2"),
    transcript_id = c("transcript1", "transcript2"),
    type = "exon",
    exon_id = c("exon1", "exon2"),
    old_width = c(6, 3)
)

input_gr_case3_extended_fixed <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(5, 20),
                              end = c(19, 25)),
    strand = "+",
    gene_id = c("gene1", "gene2"),
    transcript_id = c("transcript1", "transcript2"),
    type = "exon",
    exon_id = c("exon1", "exon2"),
    old_width = c(6, 3)
)

input_gr_case4_pot_extended <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(5, 5),
                              end = c(25, 25)),
    strand = "+",
    gene_id = c("gene1", "gene2"),
    transcript_id = c("transcript1", "transcript2"),
    type = "exon",
    exon_id = c("exon1", "exon2"),
    old_width = c(6, 18)
)

input_gr_case4_extended_fixed <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(5, 5),
                              end = c(10, 25)),
    strand = "+",
    gene_id = c("gene1", "gene2"),
    transcript_id = c("transcript1", "transcript2"),
    type = "exon",
    exon_id = c("exon1", "exon2"),
    old_width = c(6, 18)
)

input_gr_case4bis_pot_extended <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(5, 5),
                              end = c(25, 25)),
    strand = "+",
    gene_id = c("gene1", "gene2"),
    transcript_id = c("transcript1", "transcript2"),
    type = "exon",
    exon_id = c("exon1", "exon2"),
    old_width = c(6, 21)
)

input_gr_case4bis_extended_fixed <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(5, 5),
                              end = c(10, 25)),
    strand = "+",
    gene_id = c("gene1", "gene2"),
    transcript_id = c("transcript1", "transcript2"),
    type = "exon",
    exon_id = c("exon1", "exon2"),
    old_width = c(6, 21)
)

input_gr_case5_pot_extended <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(5, 20),
                              end = c(25, 25)),
    strand = "+",
    gene_id = c("gene1", "gene2"),
    transcript_id = c("transcript1", "transcript2"),
    type = "exon",
    exon_id = c("exon1", "exon2"),
    old_width = c(6, 3)
)

input_gr_case5_extended_fixed <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(5, 20),
                              end = c(19, 25)),
    strand = "+",
    gene_id = c("gene1", "gene2"),
    transcript_id = c("transcript1", "transcript2"),
    type = "exon",
    exon_id = c("exon1", "exon2"),
    old_width = c(6, 3)
)

input_gr_case6_pot_extended <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(1, 1, 30),
                              end = c(25, 10, 40)),
    strand = "+",
    gene_id = "gene1",
    transcript_id = c("transcript1", "transcript2", "transcript2"),
    type = "exon",
    exon_id = c("exon1", "exon1", "exon2"),
    old_width = c(10, 10, 11)
)

input_gr_case6bis_pot_extended <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(1, 1, 30),
                              end = c(25, 10, 40)),
    strand = "+",
    gene_id = c("gene1", "gene2", "gene2"),
    transcript_id = c("transcript1", "transcript2", "transcript2"),
    type = "exon",
    exon_id = c("exon1", "exon2", "exon3"),
    old_width = c(10, 10, 11)
)

input_gr_case7_pot_extended <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(1, 1, 30),
                              end = c(25, 10, 40)),
    strand = "+",
    gene_id = "gene1",
    transcript_id = c("transcript1", "transcript2", "transcript2"),
    type = "exon",
    exon_id = c("exon1", "exon2", "exon3"),
    old_width = c(8, 10, 11)
)

input_gr_case8_pot_extended <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(
        start = c(1, 1, 30),
        end = c(25, 10, 40)
    ),
    strand = "+",
    gene_id = c("gene1", "gene2", "gene2"),
    transcript_id = c("transcript1", "transcript2", "transcript2"),
    type = "exon",
    exon_id = c("exon1", "exon1", "exon2"),
    old_width = c(8, 10, 11)
)

input_gr_case8_extended_fixed <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(
        start = c(1, 1, 30),
        end = c(8, 10, 40)
    ),
    strand = "+",
    gene_id = c("gene1", "gene2", "gene2"),
    transcript_id = c("transcript1", "transcript2", "transcript2"),
    type = "exon",
    exon_id = c("exon1", "exon1", "exon2"),
    old_width = c(8, 10, 11)
)

input_gr_case9_pot_extended <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(5, 55),
                              end = c(25, 70)),
    strand = "+",
    gene_id = c("gene1", "gene2"),
    transcript_id = c("transcript1", "transcript2"),
    type = "exon",
    exon_id = c("exon1", "exon2"),
    old_width = c(21, 16)
)

test_that("three_prime_pos works", {
    expect_equal(three_prime_pos(input_gr),
                 c(3, 13, 20, 30))
})

test_that("five_prime_pos works", {
    expect_equal(five_prime_pos(input_gr),
                 c(1, 10, 23, 33))
})

test_that("extract_last_exon works", {
    expect_equal(extract_last_exon(input_gr, "gene_id"), input_gr[2:3])
    expect_equal(extract_last_exon(input_gr, "transcript_id"), input_gr[1:3])
    expect_equal(extract_last_exon(input_gr, "gene_id", invert = TRUE),
                 input_gr[c(1, 4)])
    expect_equal(extract_last_exon(input_gr, "transcript_id", invert  = TRUE),
                 input_gr[4])
})

test_that("extend_using_overlap works", {
    expect_equal(
        extend_using_overlap(input_gr, input_gr_to_overlap),
        small_gr_pot_extended
    )
})

test_that("adjust_for_collision works on no issue", {
    res <- adjust_for_collision(small_gr_pot_extended)
    expect_equal(
        nrow(res$pot_issues),
        0
    )
    expect_equal(
        res$new_gr,
        small_gr_pot_extended
    )
})

test_that("case1 works", {
    expect_equal(
        extend_using_overlap(extract_last_exon(input_gr_case1),
                                  input_to_overlap_case1_2_3_4_6_7_8),
        subset(
            input_gr_case1_pot_extended,
            old_width != width
        )
    )
    expect_equal(
        adjust_for_collision(input_gr_case1_pot_extended)[["new_gr"]],
        input_gr_case1_extended_fixed
    )
})

test_that("case2 works", {
    expect_equal(
        extend_using_overlap(extract_last_exon(input_gr_case2),
                                  input_to_overlap_case1_2_3_4_6_7_8),
        subset(
            input_gr_case2_pot_extended,
            old_width != width
        )
    )
    expect_equal(
        adjust_for_collision(input_gr_case2_pot_extended)[["new_gr"]],
        input_gr_case2_extended_fixed
    )
})

test_that("case3 works", {
    expect_equal(
        extend_using_overlap(extract_last_exon(input_gr_case3_5),
                                  input_to_overlap_case1_2_3_4_6_7_8),
        subset(
            input_gr_case3_pot_extended,
            old_width != width
        )
    )
    expect_equal(
        adjust_for_collision(input_gr_case3_pot_extended)[["new_gr"]],
        input_gr_case3_extended_fixed
    )
})

test_that("case4 works", {
    expect_equal(
        extend_using_overlap(extract_last_exon(input_gr_case4),
                                  input_to_overlap_case1_2_3_4_6_7_8),
        subset(
            input_gr_case4_pot_extended,
            old_width != width
        )
    )
    expect_equal(
        adjust_for_collision(input_gr_case4_pot_extended)[["new_gr"]],
        input_gr_case4_extended_fixed
    )
})


test_that("case4bis works", {
    expect_equal(
        extend_using_overlap(extract_last_exon(input_gr_case4bis),
                                  input_to_overlap_case1_2_3_4_6_7_8),
        subset(
            input_gr_case4bis_pot_extended,
            old_width != width
        )
    )
    expect_equal(
        adjust_for_collision(input_gr_case4bis_pot_extended)[["new_gr"]],
        input_gr_case4bis_extended_fixed
    )
})

test_that("case5 works", {
    expect_equal(
        extend_using_overlap(extract_last_exon(input_gr_case3_5),
                                  input_to_overlap_case5_9),
        subset(
            input_gr_case5_pot_extended,
            old_width != width
        )
    )
    expect_equal(
        adjust_for_collision(input_gr_case5_pot_extended)[["new_gr"]],
        input_gr_case5_extended_fixed
    )
})

test_that("case6 works", {
    expect_equal(
        extend_using_overlap(extract_last_exon(input_gr_case6),
                                  input_to_overlap_case1_2_3_4_6_7_8),
        subset(
            input_gr_case6_pot_extended,
            old_width != width
        )
    )
    expect_equal(
        adjust_for_collision(input_gr_case6_pot_extended)[["new_gr"]],
        input_gr_case6_pot_extended
    )
})

test_that("case6bis works", {
    expect_equal(
        extend_using_overlap(extract_last_exon(input_gr_case6bis),
                                  input_to_overlap_case1_2_3_4_6_7_8),
        subset(
            input_gr_case6bis_pot_extended,
            old_width != width
        )
    )
    expect_equal(
        adjust_for_collision(input_gr_case6bis_pot_extended)[["new_gr"]],
        input_gr_case6bis_pot_extended
    )
    expect_equal(
        sort(adjust_for_collision(
            input_gr_case6bis_pot_extended[c(2,1,3)]
            )[["new_gr"]]),
        sort(input_gr_case6bis_pot_extended)
    )
})

test_that("case7 works", {
    expect_equal(
        extend_using_overlap(extract_last_exon(input_gr_case7),
                                  input_to_overlap_case1_2_3_4_6_7_8),
        subset(
            input_gr_case7_pot_extended,
            old_width != width
        )
    )
    expect_equal(
        adjust_for_collision(input_gr_case7_pot_extended)[["new_gr"]],
        input_gr_case7_pot_extended
    )
})

test_that("case8 works", {
    expect_equal(
        extend_using_overlap(extract_last_exon(input_gr_case8),
                                  input_to_overlap_case1_2_3_4_6_7_8),
        subset(
            input_gr_case8_pot_extended,
            old_width != width
        )
    )
    expect_equal(
        adjust_for_collision(input_gr_case8_pot_extended)[["new_gr"]],
        input_gr_case8_extended_fixed
    )
})

test_that("add_new_exons works on case 5", {
    expect_equal(
        add_new_exons(
            input_gr_case5_extended_fixed,
            input_to_overlap_case5_9,
            verbose = 0
            ),
        GenomicRanges::GRanges(
            seqnames = "chr1",
            ranges = IRanges::IRanges(start = c(5, 20, 33, 45, 72),
                                      end = c(19, 25, 40, 60, 75)),
            strand = "+",
            gene_id = rep(c("gene1", "gene2"), c(1, 4)),
            transcript_id = rep(c("transcript1", "transcript2"), c(1, 4)),
            type = "exon",
            exon_id = c("exon1", "exon2", "BREW3R0000000001",
                        "BREW3R0000000002", "BREW3R0000000003"),
            old_width = c(6, 3, NA, NA, NA),
            exon_number = c(1, 1:4)
        )
    )
})

test_that("add_new_exons works on case 9", {
    expect_equal(
        add_new_exons(
            input_gr_case9_pot_extended,
            input_to_overlap_case5_9,
            verbose = 0
        ),
        GenomicRanges::GRanges(
            seqnames = "chr1",
            ranges = IRanges::IRanges(start = c(5, 33, 45, 55),
                                      end = c(25, 40, 54, 70)),
            strand = "+",
            gene_id = rep(c("gene1", "gene2"), c(3, 1)),
            transcript_id = rep(c("transcript1", "transcript2"), c(3, 1)),
            type = "exon",
            exon_id = c("exon1", "BREW3R0000000001", "BREW3R0000000002",
                        "exon2"),
            old_width = c(21, NA, NA, 16),
            exon_number = c(1:3, 1)
        )
    )
})


test_that("add_new_exons works when run twice", {
    expect_equal(
        sort(add_new_exons(
            GenomicRanges::GRanges(
                seqnames = "chr1",
                ranges = IRanges::IRanges(start = c(5, 33, 45, 55),
                                          end = c(25, 40, 54, 70)),
                strand = "+",
                gene_id = rep(c("gene1", "gene2"), c(3, 1)),
                transcript_id = rep(c("transcript1", "transcript2"), c(3, 1)),
                type = "exon",
                exon_id = c("exon1", "BREW3R0000000001", "BREW3R0000000002",
                            "exon2"),
                exon_number = c(1:3, 1)
            ),
            GenomicRanges::GRanges(
                seqnames = "chr1",
                ranges = IRanges::IRanges(start = c(60, 92),
                                          end = c(70, 95)),
                strand = "+",
                gene_id = "geneA",
                transcript_id = "transcriptA",
                type = "exon",
                exon_id = c("exonA", "exonB")
            ),
            verbose = 0
        )),
        sort(GenomicRanges::GRanges(
            seqnames = "chr1",
            ranges = IRanges::IRanges(start = c(5, 33, 45, 55, 92),
                                      end = c(25, 40, 54, 70, 95)),
            strand = "+",
            gene_id = rep(c("gene1", "gene2"), c(3, 2)),
            transcript_id = rep(c("transcript1", "transcript2"), c(3, 2)),
            type = "exon",
            exon_id = c("exon1", "BREW3R0000000001", "BREW3R0000000002",
                        "exon2", "BREW3R0000000003"),
            exon_number = c(1:3, 1:2)
        ))
    )
})
