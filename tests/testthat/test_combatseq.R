test_that("test if combatSeq does its job", {

    data(ribo_toy)
    ribo_corrected <- rRMSAnalyzer::adjust_bias(ribo_toy,"run")
    ribo_corrected_m <- extract_data(ribo_corrected)
    
    expected_correction <- readRDS(testthat::test_path("testdata","cseq.rds"))
    
    testthat::expect_equal(ribo_corrected_m,expected_correction)
    
})