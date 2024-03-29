# Load the fixture
annotated_gene_lists_fixture <- readRDS('./fixtures/AnnotatedGeneLists.rds')
# background_fixture
background_fixture <- annotated_gene_lists_fixture@inbuilt[["background"]]@ENTREZID

test_that("GO_enrichment_analysis returns expected results", {
  skip_if_not_installed("clusterProfiler")
  skip_if_offline()

  # Run the function
  result <- GO_enrichment_analysis(
    annotated_gene_lists = annotated_gene_lists_fixture,
    background = background_fixture,
    organism = "org.Hs.eg.db",
    keyType = "ENTREZID",
    minGSSize = 10
  )


  # Test 1 Check if the result is a list
  expect_is(result, "list")

  # Test 2 Check if the list is not empty
  expect_true(length(result) > 0)

  expected_names <- c('Lasso', 'Univariate', 'RandomForest', 'boruta', 'DEG_rural', 'DEG_urban')
  expect_equal(names(result), expected_names)

})

test_that("GO_enrichment_analysis throws an error for invalid list_type", {

  skip_if_not_installed("clusterProfiler")
  skip_if_offline()

  # Test 3
  expect_error(
    GO_enrichment_analysis(
      annotated_gene_lists = annotated_gene_lists_fixture,
      list_type = "invalid_type",
      background = background_fixture,
      organism = "org.Hs.eg.db",
      keyType = "ENTREZID",
      minGSSize = 10
    ),
    "list_type should be 'inbuilt' or 'permutation'"
  )
})
