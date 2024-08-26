##################################################
# install and load dependencies
if (!require("devtools", quietly = TRUE))
  install.packages("devtools", repos = "http://cran.us.r-project.org")
if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr", repos = "http://cran.us.r-project.org")
if (!require("ggfortify", quietly = TRUE))
  install.packages("ggfortify", repos = "http://cran.us.r-project.org")
# install BgeeDB from the devel branch of the github repository as Bgee 15.2 data are not yet available on Bioconductor.
# The master branch or last version of the Bioconductor release should be used once the new Bioconductor release (3.20) is available
if (!require("BgeeDB", quietly = TRUE))
  devtools::install_github(repo = "https://github.com/BgeeDB/BgeeDB_R", ref = "devel")

library(BgeeDB)
library(dplyr)
library(ggfortify)

# source the file containing functions used to generate the
# PCA input file
source("/workspace/Biodiversity-Bioinformatics/scripts/retrieveExpressionOrthologsFromBgee.R")


##############################################
## Retrieve one-tot-one orthologs from Bgee ##
##############################################
orthologs_one_to_one <- listOrthologs(mandatorySpecies = c(9031, 9986, 10090, 10116, 9544, 9606, 13616), onlyOneToOne = TRUE)
# export the ono-to-one orthologs as a tsv file
write.table(x = orthologs_one_to_one, file = "/workspace/Biodiversity-Bioinformatics/data/orthologs_9031_9544_9606_9986_10090_10116.tsv",
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

##########################################################################
## load processed expression values from kaessmann 2019 paper from Bgee ##
##########################################################################
# This step download processed expression values from 1890 bulk RNA-Seq
# libraries. It is long to process
expression_kaessmann_2019_bgee <- loadExpressionKaessmann2019()
# Save processed expression for all gene families of selected species
# in order not to have to redownload all the data
saveRDS(object = expression_kaessmann_2019_bgee, file = "/workspace/Biodiversity-Bioinformatics/data/expression_kaessmann_2019_bgee.Rds")

# Load data from a file if you already donwloaded it
# To do so comment the 2 lines of code loading data from Bgee and saving it to a file
# and uncomment the line below
#expression_kaessmann_2019_bgee <- readRDS(file = "/workspace/Biodiversity-Bioinformatics/data/expression_kaessmann_2019_bgee.Rds")

###############################################
## Remap dev. stages and anatomical entities ##
###############################################
# remap dev. stages to embryo and post embry
kaessmann_remapped <- remapDevStages(expression_kaessmann_2019_bgee)
# remap anatomical entities from Bgee to higher level terms when too precise
kaessmann_remapped <- remapAnatEntities(kaessmann_remapped)

##########################################################################################
## Transform expression data to a format compatible to the input of the prcomp function ##
##########################################################################################
transformed_expression <- transfrom_pca(orthologs_one_to_one, kaessmann_remapped)

#remove opossum data because they contained a lot of NAs
transformed_expression <- transformed_expression[transformed_expression$speciesId != 13616,]

#remove gene family for which some species have NA TPM value
transformed_expression <- transformed_expression %>% select_if(~ !any(is.na(.)))

# save data transformed and ready to use as input of PCA in order to allow student to directly
# import those data
write.table(x = transformed_expression, file = "/workspace/Biodiversity-Bioinformatics/data/pca_dataset.tsv", sep = "\t",
 row.names = FALSE, col.names = TRUE)
