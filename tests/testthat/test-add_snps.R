test_that("Only accepts SCE object", {
  expect_error(add_snps(snps,snps))
})
