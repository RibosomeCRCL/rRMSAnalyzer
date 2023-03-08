test_that("RiboClass creation",{

    path <- system.file("extdata", package="rRMSAnalyzer")
    
    # Metadata loaded from file
    expect_no_error(ribo <- load_ribodata(
                        #data & metadata files path
                        count_path = file.path(path,"miniglioma/"),
                        metadata = file.path(path,"metadata.csv"),
                        # data & metadata files separator
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
                        ncores = 1))

    # Metadata previously loaded in a dataframe

    met <- read.csv(paste0(path,"/metadata.csv"),sep = ",")
    ribo <- load_ribodata(
                        #data & metadata files path
                        count_path = paste0(path,"/miniglioma"),
                        metadata = met,
                        # data & metadata files separator
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

    expect_s3_class(ribo,"RiboClass")


})