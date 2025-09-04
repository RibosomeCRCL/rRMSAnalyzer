test_that("test if plots creation do not fail", {
    data(ribo_toy)
  
  custom_anno <- data.frame(
    rnapos = c(15,76,100,401),
    rna = c("5.8S","5.8S","18S","28S"),
    nomenclature = c("5.8S_Um14","5.8S_Gm75","18S_Am99","28S_Am391")
  )
  
  ribo_toy <- rename_rna(ribo_toy)
  ribo_toy <- annotate_site(ribo_toy,custom_anno)
    
  expect_no_error(rRMSAnalyzer::plot_pca(ribo_toy, "condition"))
  expect_no_error(rRMSAnalyzer::plot_coa(ribo_toy, "condition"))
  
  expect_no_error(rRMSAnalyzer::boxplot_count(ribo_toy))
  expect_no_error(rRMSAnalyzer::boxplot_cscores(ribo_toy))
  expect_no_error(rRMSAnalyzer::plot_rlc(ribo_toy))
  
  expect_no_error(rRMSAnalyzer::plot_heatmap(ribo_toy,
                                             color_col = c("run","condition")))
  expect_no_error(rRMSAnalyzer::plot_heatmap_corr(ribo_toy,"count","run"))
  
  expect_no_error(rRMSAnalyzer::plot_counts_env(ribo_toy,"5S",50))
  expect_no_error(rRMSAnalyzer::plot_counts_env(ribo_toy,"5S",50,c("S1","S2")))
  
  expect_no_error(rRMSAnalyzer::plot_diff_sites(ribo_toy,"condition"))
  
                  
  

})