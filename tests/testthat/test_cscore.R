test_that("compute_cscore outputs the expected RiboClass with default params",{
    data(ribo_toy)
    # Get the already computed c-score from ribo_toy
    ribo_expected_cscore <- rRMSAnalyzer::extract_data(ribo_toy)
    
    # Recompute the same C-score and check if it matches with the old one
    ribo_new_cscore <- rRMSAnalyzer::compute_cscore(ribo_toy)
    ribo_new_matrix <- rRMSAnalyzer::extract_data(ribo_new_cscore)
    expect_equal(ribo_new_matrix,ribo_expected_cscore)
})