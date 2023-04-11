test_that("Annotation works correctly",{
  
  data("ribo_toy")
  
  ribo_toy <- rRMSAnalyzer::rename_rna(ribo_toy)
  
  data("human_methylated")
  ribo_toy <- rRMSAnalyzer::annotate_site(ribo_toy,human_methylated)
  ribom <- rRMSAnalyzer::extract_data(ribo_toy,only_annotated = T)
  
  expected_ribom <- readRDS(testthat::test_path("testdata","anno.rds"))
  
  testthat::expect_equal(ribom,expected_ribom)
  
  
})