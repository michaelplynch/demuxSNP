test_that("Only accepts SCE object", {
  expect_error(top_genes(snps))
})
