test_that("test if plots creation do not fail", {
    data(ribo_toy)
    
    expect_no_error(rRMSAnalyzer::plot_pca(ribo_toy, "condition"))
    expect_no_error(rRMSAnalyzer::plot_coa(ribo_toy, "condition"))


})