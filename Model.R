#----------------------------------------------------------------------------------
#
#   rRMSAnalyzer: script d'analyse de RiboMethSeq | Humain et Murin
#   par VM le 22/06/2023
#   à l'aide des packages rRMSAnalyzer & rRMSReports & rRMSFonctions
#
#   Pour de l'aide plus spécifique: https://ribosomecrcl.github.io/rRMSAnalyzer/
#
#----------------------------------------------------------------------------------

#----------------------------------------------------------------------------------
# 0. Manipulation sous R préalable
#----------------------------------------------------------------------------------


#0.0. Créer un nouveau projet et préparation de l'environnement de travail

#Pour des premières prises en main: aller directment à #0.1 et utiliser le ribo_toy puis créer un riboClass à partir des données de 5FU fournies

#Pour analyser ses propres données:
#Aller dans /file/new projet/ et créer un nouveau dossier
#Dans ce nouveau dossier, copier/coller les fichiers du dossier "rRMSAnalyser_template":
#metadata.csv + rRMSAnalyzer_template.R
#si c'est une analyse souris, prendre aussi "mouse.methylated.rda" et "mouse.suspected.rda"
#De plus, n'oublier pas de
#copier/coller le dossier avec les 5'end counts
#ouvrir le metadata.csv sous Excel,  remplir le fichier puis l'enregistrer en .csv avec séparation en ";" dans le dossier de comptage


#0.1. Mise à jour des packages
library(devtools)

#0.1.1. RMSAnalyzer
devtools::install_github("RibosomeCRCL/rRMSAnalyzer")

library(rRMSAnalyzer)
library(rRMSReports)

#----------------------------------------------------------------------------------
# 1. Création du RiboClass
#----------------------------------------------------------------------------------


# 1.2. Préparation d'un jeu de données personnelles.

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

# Vérification de l'importation de données avec PCA
plot_pca(ribo,"lib") #changer le nom de la colonne en fonction du message d'erreur

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
# 2. Etape de QC
#----------------------------------------------------------------------------------

# 2.1. Réalisation du QC_report pour toutes les données

report_qc(ribo = ribo, library_col = "lib")

# 2.2. Ajustement par ComBat-Seq si nécessaire

plot_pca(ribo,"library")                #ribo_toy library_col = "run"

ribo_adj <- adjust_bias(ribo,"library") 

plot_pca(ribo_adj,"library")            

# 2.3. Réalisation du QC_report pour toutes les données ajustées

report_qc(ribo = ribo_adj, library_col = "library", project_name = "xxx", comments = "./commentaire_QC.Rmd") 

#----------------------------------------------------------------------------------
# 3. Retirer les échantillons problématiques ou qui ne font pas partis de l'étude
#----------------------------------------------------------------------------------

# 3.1. Garder un faible nombre d'échantillons.
ribo_adj_small <- keep_ribo_samples(ribo_adj,c("sample1","sample2","sample3","sample4"))
                                 
# 3.1.bis Enlever un faible nombre d'échantillons                                
ribo_adj_small <- remove_ribo_samples(ribo_adj,c("xxxx","xxxccc","RNA_ref-1","RNA_ref-2"))

# 3.3. Uniformiser le nom du riboClass de cette étape
ribo_adj <- ribo_adj_small 


#----------------------------------------------------------------------------------
# 4A. Annotation des sites de 2'Ome des ARNr HUMAINS
#----------------------------------------------------------------------------------

data("human_methylated")
cat("human_methylated's rna names: ", unique(human_methylated$rRNA),"\n")
cat("ribo's rna names: ", as.character(ribo$rna_names$current_name)) # changer le nom du RiboClass

ribo_adj_name <- rename_rna(ribo,
                            new_names = c("5S", "5.8S", "18S", "28S"))

ribo_adj_annot <- annotate_site(ribo_adj_name,
                                annot = human_methylated,
                                anno_rna = "rRNA",
                                anno_pos = "Position",
                                anno_value = "Nomenclature")

#----------------------------------------------------------------------------------
# 4B. dans un deuxième temps, refaire toute l'analyse avec les sites suspectés si besoin
#----------------------------------------------------------------------------------

load(file = "/home/bioinfo/Bureau/RMS/mouse.suspected.rda") 

ribo_adj_annot <- annotate_site(ribo_adj_annot, human.suspected) #vérifier si les noms des ARNr sont les mêmes ou non


ribo_adj_annot <- annotate_site(ribo_adj_annot,
                                annot = human.suspected,
                                anno_rna = "rRNA",
                                anno_pos = "Position",
                                anno_value = "Nomenclature")

#----------------------------------------------------------------------------------
# 4Abis. Annotation des sites de 2'Ome des ARNr MURINS
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
# 4Bis. dans un deuxième temps, refaire toute l'analyse avec les sites suspectés si besoin
#------------------------------------

load(file = "/home/bioinfo/Bureau/RMS/mouse.suspected.rda") 

ribo_adj_annot <- annotate_site(ribo_adj_annot, mouse.suspected) #vérifier si les noms des ARNr sont les mêmes ou non


ribo_adj_annot <- annotate_site(ribo_adj_annot,
                                annot = mouse.suspected,
                                anno_rna = "rRNA",
                                anno_pos = "Position",
                                anno_value = "Nomenclature")

#----------------------------------------------------------------------------------
# 5. Rapports automatiques
#----------------------------------------------------------------------------------

# A refaire autant de fois qu'il y a de variables (A, B, C ...). Ex de metadata :

# filename	              samplename	      variableA	  variableB	    variableC	  condition	  traitement	lib
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

# Si vous n'avez qu'une seule colonne "condition" le fonctionnement est le même. 
# --------------------------------------A------------------------------------------

# 5.1 filtrer avec les métadonnées pour ne garder que les échantillons présents dans la colonne "variableA"
kept_samples <- ribo_adj_annot$metadata %>%
  dplyr::filter(!is.na(variableA)) %>% # ne garde que les lignes qui ne sont pas rempli avec "NA"
  dplyr::pull(samplename)

ribo_adj_annot_variableA <- keep_ribo_samples(ribo_adj_annot,kept_samples)

# créer le tableau de comparaison nécessaire pour le rapport diff_sites
comparisons <- tibble::tibble(
  comp = c("comp1", "comp2"),
  ctrl = c("ctrl", "ctrl"),
  cases = c("P1", "P2")
)

# 5.2. Analyse du profil global de 2'Ome des ARNr et des xx sites les plus variants
report_2ome_sites(ribo_adj_annot_variableA,"variableA", project_name = "nom", comments = "./commentaire_2ome_variableA.Rmd") 

# 5.3. Analyse des différences de niveau de 2'Ome des ARNr
report_diff_sites(ribo_adj_annot_variableA, condition_col = "variableA", project_name = "nom", comparisons = comparisons, comments = "./commentaire_diff_site_variableA.Rmd") 

# --------------------------------------B------------------------------------------

# 5.1 filtrer avec les métadonnées pour ne garder que les échantillons présents dans la colonne "variableB"
kept_samples <- ribo_adj_annot$metadata %>%
  dplyr::filter(!is.na(variableB)) %>%
  dplyr::pull(samplename)

ribo_adj_annot_variableB <- keep_ribo_samples(ribo_adj_annot,kept_samples)

# créer le tableau de comparaison nécessaire pour le rapport diff_sites
comparisons <- tibble::tibble(
  comp = c("comp1"),
  ctrl = c("ctrl"),
  cases = c("P2")
)

# 5.2. Analyse du profil global de 2'Ome des ARNr et des xx sites les plus variants
report_2ome_sites(ribo_adj_annot_variableB,"variableB", project_name = "nom", comments = "./commentaire_2ome_variableB.Rmd") 

# 5.3. Analyse des différences de niveau de 2'Ome des ARNr
report_diff_sites(ribo_adj_annot_variableB, condition_col = "variableB", project_name = "nom", comparisons = comparisons, comments = "./commentaire_diff_site_variableB.Rmd") 

# --------------------------------------C------------------------------------------

# 5.1 filtrer avec les métadonnées pour ne garder que les échantillons présents dans la colonne "variableC"
kept_samples <- ribo_adj_annot$metadata %>%
  dplyr::filter(!is.na(variableC)) %>%
  dplyr::pull(samplename)

ribo_adj_annot_variableC <- keep_ribo_samples(ribo_adj_annot,kept_samples)

# créer le tableau de comparaison nécessaire pour le rapport diff_sites
comparisons <- tibble::tibble(
  comp = c("comp1"),
  ctrl = c("ctrl"),
  cases = c("P2")
)


# 5.2. Analyse du profil global de 2'Ome des ARNr et des xx sites les plus variants
report_2ome_sites(ribo_adj_annot_variableC,"variableC", project_name = "nom", comments = "./commentaire_2ome_variableC.Rmd") 

# 5.3. Analyse des différences de niveau de 2'Ome des ARNr
report_diff_sites(ribo_adj_annot_variableC, condition_col = "variableC", project_name = "nom", comparisons = comparisons, comments = "./commentaire_diff_site_variableC.Rmd") 


#----------------------------------------------------------------------------------
# 6. Vérifier les environnements locaux des sites d'intérêts
#----------------------------------------------------------------------------------


# 6.1 Pour tous les échantillons
plot_counts_env(ribo_adj_annot, rna= "18S", pos = 1442, samples = "all", condition = "condition")

# 6.2 Pour certain échantillons
plot_counts_env(ribo_adj_annot, rna = "18S", pos = 1442, samples = c("sample1", "sample2"))


#----------------------------------------------------------------------------------
# 7. Export des données
#----------------------------------------------------------------------------------

ribo_df <- extract_data(ribo_adj_annot,
                        col = "cscore",
                        only_annotated = TRUE)

write.csv(ribo_df, "/chemin/vers/dossier/ribo_df.csv") #indiquer le chemin

#----------------------------------------------------------------------------------
