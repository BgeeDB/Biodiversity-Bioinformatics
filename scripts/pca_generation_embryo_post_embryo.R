if (!require('ggfortify', quietly = TRUE)) {
  install.packages("ggfortify", repos = "http://cran.us.r-project.org")
}
library(ggfortify)

pca_file_name <- "/workspace/Biodiversity-Bioinformatics/data/pca_dataset.tsv"

expression_dataset <- NULL
if(file.exists(pca_file_name)) {
  expression_dataset = read.table(file = pca_file_name, header = TRUE, sep = "\t")
} else {
  expression_dataset = read.table(file = gzfile(paste0(pca_file_name,".gz")), header = TRUE, sep = "\t")
}

#replace type of speciesId from numeric to character
expression_dataset$speciesId <- as.character(expression_dataset$speciesId)

#filter on embryo / post-embryo
expression_post_embryo <- 
  expression_dataset[expression_dataset$devStageName == 'post-embryonic stage',]
expression_embryo <- 
  expression_dataset[expression_dataset$devStageName == 'embryo stage',]

#run the pca
pca_res_post_embryo <- prcomp(expression_post_embryo[5:length(expression_post_embryo)], scale. = TRUE)
pca_res_embryo <- prcomp(expression_embryo[5:length(expression_embryo)], scale. = TRUE)

#plot results per species
jpeg(file="/workspace/Biodiversity-Bioinformatics/figures/Species_PCA_post_embryo.jpeg")
autoplot(pca_res_post_embryo, data = expression_post_embryo, colour = 'speciesId')
dev.off()
jpeg(file="/workspace/Biodiversity-Bioinformatics/figures/Species_PCA_embryo.jpeg")
autoplot(pca_res_embryo, data = expression_embryo, colour = 'speciesId')
dev.off()

#plot results per organ
jpeg(file="/workspace/Biodiversity-Bioinformatics/figures/Organ_PCA_post_embryo.jpeg")
autoplot(pca_res_post_embryo, data = expression_post_embryo, colour = 'anatEntityName')
dev.off()
jpeg(file="/workspace/Biodiversity-Bioinformatics/figures/Organ_PCA_embryo.jpeg")
autoplot(pca_res_embryo, data = expression_embryo, colour = 'anatEntityName')
dev.off()