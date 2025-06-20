#----------------------------------------------------------------------------------
#
#   rRMSAnalyzer: analysis script of RiboMethSeq | Human et Murine
#   by llyson Moureaux inspired by Virginie Marcel the 19/06/2025
#   thankes to rRMSAnalyzer & rRMSReports packages & rRMSFonctions
#
#   For more specific help: https://ribosomecrcl.github.io/rRMSAnalyzer/
#
#----------------------------------------------------------------------------------

#----------------------------------------------------------------------------------
# 0. Before starting analysis
#----------------------------------------------------------------------------------


#0.0. Create a new Project and prepare working environment

#go to /file/new project/ create a new folder
#In this folder, copy and paste files of directory "rRMSAnalyser_template":
#metadata.csv + rRMSAnalyzer_template.R
#if it's mouse analysis take "mouse.methylated.rda" and "mouse.suspected.rda"
#Don't forget to put your 5'end counts file in this directory.
#open metadata.csv in Excel, fill in the file and save it under .csv with ";" as separator.


#0.1. Package update
library(devtools)

#0.1.1. RMSAnalyzer and RMSReport
devtools::install_github("RibosomeCRCL/rRMSAnalyzer")

library(rRMSAnalyzer)
library(rRMSReports)

#----------------------------------------------------------------------------------
# 1. RiboClass creation
#----------------------------------------------------------------------------------

# 1.2. Personal data preparation

ribo <- load_ribodata(
  #data & metadata files path
  count_path = "~/Bureau/RMS/chemin/counts/",
  metadata = "~/Bureau/RMS/chemin/metadata_RMS2b.csv",# data & metadata files separator
  count_sep = "\t",
  metadata_sep = ",",
  # count data parameters :
  count_header = FALSE,
  count_value = 3,
  count_rnaid = 1,
  count_pos = 2,
  # Metadata parameters :
  metadata_key = "filename",
  metadata_id = "samplename",
  # c-score parameters :
  flanking = 6,
  method = "median",
  ncores = 1)

# Check importation with a PCA
plot_pca(ribo,"lib") #change column name of the library if necessary

#Plot coa
plot_coa(
  ribo,
  color_col = NULL,
  axes = c(1, 2),
  only_annotated = FALSE,
  title = "default",
  subtitle = "default",
  object_only = FALSE
)
#----------------------------------------------------------------------------------
# 2. QC step
#----------------------------------------------------------------------------------

# 2.1. QC_report for all data

report_qc(ribo = ribo, library_col = "lib")

# 2.2. adjustment ComBat-Seq if necessary

plot_pca(ribo,"library")                #ribo_toy library_col = "run"

ribo_adj <- adjust_bias(ribo,"library") 

plot_pca(ribo_adj,"library")            

# 2.3. QC_report for adjusted data

report_qc(ribo = ribo_adj, library_col = "library", project_name = "xxx", comments = "./commentaire_QC.Rmd") 

#----------------------------------------------------------------------------------
# 3. Remove problematic or not used samples
#----------------------------------------------------------------------------------

# 3.1. Keep a small number of samples
ribo_adj_small <- keep_ribo_samples(ribo_adj,c("sample1","sample2","sample3","sample4"))
                                 
# 3.1.bis  remove a small number of samples                        
ribo_adj_small <- remove_ribo_samples(ribo_adj,c("xxxx","xxxccc","RNA_ref-1","RNA_ref-2"))

# 3.3. Uniformisation of Riboclass name
ribo_adj <- ribo_adj_small 


#----------------------------------------------------------------------------------
# 4A. Annotation des sites de 2'Ome des ARNr HUMAINS
#----------------------------------------------------------------------------------

data("human_methylated")
cat("human_methylated's rna names: ", unique(human_methylated$rRNA),"\n")
cat("ribo's rna names: ", as.character(ribo$rna_names$current_name)) # change name of RiboClass if necessary

ribo_adj_name <- rename_rna(ribo, #change name of RiboClass if necessary for ribo_adj if you adjusted the data
                            new_names = c("5S", "5.8S", "18S", "28S"))

ribo_adj_annot <- annotate_site(ribo_adj_name, 
                                annot = human_methylated,
                                anno_rna = "rRNA",
                                anno_pos = "Position",
                                anno_value = "Nomenclature")

#----------------------------------------------------------------------------------
# 4B. redo all the analysis for adjusted data if necessary (replace step 4A then ABis by 4B then BBis)
#----------------------------------------------------------------------------------

load(file = "/home/bioinfo/Bureau/RMS/mouse.suspected.rda") 

ribo_adj_annot <- annotate_site(ribo_adj_annot, human.suspected) #vérifier si les noms des ARNr sont les mêmes ou non


ribo_adj_annot <- annotate_site(ribo_adj_annot,
                                annot = human.suspected,
                                anno_rna = "rRNA",
                                anno_pos = "Position",
                                anno_value = "Nomenclature")

#----------------------------------------------------------------------------------
# 4Abis. 2'Ome sites of MURINS rRNA annotation
#------------------------------------

load(file = "/home/bioinfo/Bureau/RMS/mouse.methylated.rda") 

ribo_adj_annot <- rename_rna(ribo,
                             new_names = c("5S", "5.8S", "18S", "28S")) 

ribo_adj_annot <- annotate_site(ribo_adj_annot,
                                annot = mouse.methylated,
                                anno_rna = "rRNA",
                                anno_pos = "Position",
                                anno_value = "Nomenclature")
#------------------------------------
# 4Bis. redo all the analysis for adjusted data if necessary (replace step 4A then ABis by 4B then BBis)
#------------------------------------

load(file = "/home/bioinfo/Bureau/RMS/mouse.suspected.rda") 

ribo_adj_annot <- annotate_site(ribo_adj_annot, mouse.suspected) #vérifier si les noms des ARNr sont les mêmes ou non


ribo_adj_annot <- annotate_site(ribo_adj_annot,
                                annot = mouse.suspected,
                                anno_rna = "rRNA",
                                anno_pos = "Position",
                                anno_value = "Nomenclature")

#----------------------------------------------------------------------------------
# 5. Reports
#----------------------------------------------------------------------------------

# To launch as many time as the number of comparisons asked(A, B, C ...). Ex of metadata:

# filename	              samplename	      variableA	  variableB	    variableC	  condition	  treatment	  lib
# S26_R1.5_counts.csv     variableB_N2_P2	  NA	        P2	          NA	        variableB	  P2	        L2
# S27_R1.5_counts.csv	    variableA_N2_P1	  P1	        NA	          NA	        variableA	  P1	        L2
# S28_R1.5_counts.csv	    variableC_D4_ctrl	NA	        NA	          ctrl	      variableC	  ctrl	      L2
# S30_R1.5_counts.csv	    variableB_N1_P2	  NA	        P2	          NA	        variableB	  P2	        L2
# S31_R1.5_counts.csv	    variableC_D2_P2	  NA	        NA	          P2	        variableC	  P2	        L2
# S32_R1.5_counts.csv	    variableA_N3_P1	  P1	        NA	          NA	        variableA	  P1        	L2
# S33_R1.5_counts.csv	    variableA_N3_ctrl	ctrl	      NA	          NA	        variableA	  ctrl	      L2
# S34_R1.5_counts.csv	    variableC_D1_ctrl	NA	        NA	          ctrl	      variableC	  ctrl       	L2
# S35_R1.5_counts.csv	    variableA_N1_P2	  P2	        NA	          NA        	variableA	  P2	        L2
# S36_R1.5_counts.csv	    variableA_N2_ctrl	ctrl	      NA	          NA	        variableA	  ctrl	      L2
# S37_R1.5_counts.csv	    variableC_D2_ctrl	NA	        NA	          ctrl	      variableC	  ctrl	      L2
# S38_R1.5_counts.csv	    variableA_N2_P2	  P2	        NA	          NA	        variableA	  P2	        L2
# S39_R1.5_counts.csv	    variableC_D1_P2	  NA	        NA	          P2	        variableC	  P2	        L2
# S40_R1.5_counts.csv	    variableB_N1_ctrl	NA	        ctrl	        NA	        variableB	  ctrl	      L2
# S41_R1.5_counts.csv	    variableA_N1_P1	  P1	        NA	          NA	        variableA	  P1	        L2
# S42_R1.5_counts.csv	    variableA_N1_ctrl	ctrl	      NA	          NA	        variableA	  ctrl	      L2
# S43_R1.5_counts.csv	    variableC_D3_P2	  NA	        NA	          P2	        variableC	  P2	        L2
# S44_R1.5_counts.csv	    variableB_N3_ctrl	NA	        ctrl	        NA	        variableB	  ctrl        L2
# S45_R1.5_counts.csv	    variableA_N3_P2	  P2	        NA	          NA	        variableA	  P2	        L2

# If only one column condition, it is still the same thing.
# --------------------------------------A------------------------------------------

# 5.1 filter with metadata to keep only samples present in the colum "variableA"
kept_samples <- ribo_adj_annot$metadata %>%
  dplyr::filter(!is.na(variableA)) %>% # keep lines that are not "NA"
  dplyr::pull(samplename)

ribo_adj_annot_variableA <- keep_ribo_samples(ribo_adj_annot,kept_samples)

# create the necessary comparison table for the diff_sites report
comparisons <- tibble::tibble(
  comp = c("comp1", "comp2"),
  ctrl = c("ctrl", "ctrl"),
  cases = c("P1", "P2")
)

# 5.2. 2'Ome of rRNAs analysis of the global profile  and the xx most variant sites
report_2ome_sites(ribo_adj_annot_variableA,"variableA", project_name = "nom", comments = "./commentaire_2ome_variableA.Rmd") 

# 5.3. Analysis of rRNAs 2'Ome level differences
report_diff_sites(ribo_adj_annot_variableA, condition_col = "variableA", project_name = "nom", comparisons = comparisons, comments = "./commentaire_diff_site_variableA.Rmd") 

# --------------------------------------B------------------------------------------

# 5.1 filter with metadata to keep only samples present in the colum "variableB"
kept_samples <- ribo_adj_annot$metadata %>%
  dplyr::filter(!is.na(variableB)) %>%
  dplyr::pull(samplename)

ribo_adj_annot_variableB <- keep_ribo_samples(ribo_adj_annot,kept_samples)

# create the necessary comparison table for the diff_sites report
comparisons <- tibble::tibble(
  comp = c("comp1"),
  ctrl = c("ctrl"),
  cases = c("P2")
)

# 5.2. 2'Ome of rRNAs analysis of the global profile  and the xx most variant sites
report_2ome_sites(ribo_adj_annot_variableB,"variableB", project_name = "nom", comments = "./commentaire_2ome_variableB.Rmd") 

# 5.3.  Analysis of rRNAs 2'Ome level differences
report_diff_sites(ribo_adj_annot_variableB, condition_col = "variableB", project_name = "nom", comparisons = comparisons, comments = "./commentaire_diff_site_variableB.Rmd") 

# --------------------------------------C------------------------------------------

# 5.1 filter with metadata to keep only samples present in the colum "variableC"
kept_samples <- ribo_adj_annot$metadata %>%
  dplyr::filter(!is.na(variableC)) %>%
  dplyr::pull(samplename)

ribo_adj_annot_variableC <- keep_ribo_samples(ribo_adj_annot,kept_samples)

# create the necessary comparison table for the diff_sites report
comparisons <- tibble::tibble(
  comp = c("comp1"),
  ctrl = c("ctrl"),
  cases = c("P2")
)


# 5.2. 2'Ome of rRNAs analysis of the global profile  and the xx most variant sites
report_2ome_sites(ribo_adj_annot_variableC,"variableC", project_name = "nom", comments = "./commentaire_2ome_variableC.Rmd") 

# 5.3.  Analysis of rRNAs 2'Ome level differences
report_diff_sites(ribo_adj_annot_variableC, condition_col = "variableC", project_name = "nom", comparisons = comparisons, comments = "./commentaire_diff_site_variableC.Rmd") 


#----------------------------------------------------------------------------------
# 6. Check local environments of sites of interest
#----------------------------------------------------------------------------------


# 6.1 For all samples
plot_counts_env(ribo_adj_annot, rna= "18S", pos = 1442, samples = "all", condition = "condition")

# 6.2 For some samples
plot_counts_env(ribo_adj_annot, rna = "18S", pos = 1442, samples = c("sample1", "sample2"))


#----------------------------------------------------------------------------------
# 7. data export
#----------------------------------------------------------------------------------

ribo_df <- extract_data(ribo_adj_annot,
                        col = "cscore",
                        only_annotated = TRUE)

write.csv(ribo_df, "/chemin/vers/dossier/cscore_df.csv") #give path

#----------------------------------------------------------------------------------
