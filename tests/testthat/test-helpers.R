test_that("remove_singletons drops taxa with total abundance <= 1", {
  otu <- matrix(c(5, 0, 1, 0,
                  3, 2, 0, 0,
                  0, 0, 0, 1),
                nrow = 3, byrow = TRUE)
  rownames(otu) <- paste0("OTU", 1:3)
  colnames(otu) <- paste0("Sample", 1:4)
  ps <- phyloseq::phyloseq(phyloseq::otu_table(otu, taxa_are_rows = TRUE))

  ps_filtered <- app_env$remove_singletons(ps)

  expect_true(all(phyloseq::taxa_sums(ps_filtered) > 1))
  expect_false("OTU3" %in% phyloseq::taxa_names(ps_filtered))
})

test_that("safe_tax_glom collapses at a valid rank and no-ops on a missing rank", {
  otu <- matrix(c(10, 5,
                  0, 8),
                nrow = 2, byrow = TRUE)
  rownames(otu) <- c("OTU1", "OTU2")
  colnames(otu) <- c("SampleA", "SampleB")
  tax <- matrix(c("Firmicutes", "GenusA",
                  "Firmicutes", "GenusB"),
                nrow = 2, byrow = TRUE)
  rownames(tax) <- c("OTU1", "OTU2")
  colnames(tax) <- c("Phylum", "Genus")
  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(otu, taxa_are_rows = TRUE),
    phyloseq::tax_table(tax)
  )

  glommed <- app_env$safe_tax_glom(ps, "Phylum")
  expect_equal(phyloseq::ntaxa(glommed), 1)

  unchanged <- app_env$safe_tax_glom(ps, "Species")
  expect_equal(phyloseq::ntaxa(unchanged), phyloseq::ntaxa(ps))
})

test_that("compositional transform makes each sample sum to 1", {
  otu <- matrix(c(10, 30,
                  20, 20),
                nrow = 2, byrow = TRUE)
  rownames(otu) <- c("OTU1", "OTU2")
  colnames(otu) <- c("SampleA", "SampleB")
  ps <- phyloseq::phyloseq(phyloseq::otu_table(otu, taxa_are_rows = TRUE))

  ps_comp <- app_env$compositional(ps)
  sums <- phyloseq::sample_sums(ps_comp)

  expect_equal(unname(sums), c(1, 1), tolerance = 1e-8)
})
