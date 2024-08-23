if (!require('ggfortify', quietly = TRUE)) {
  install.packages("ggfortify", repos = "http://cran.us.r-project.org")
}
library(ggfortify)

expression_dataset = read.table(file = gzfile("/workspace/Biodiversity-Bioinformatics/data/pca_dataset.tsv.gz"), header = TRUE, sep = "\t")

#replace type of speciesId from numeric to character
expression_dataset$speciesId <- as.character(expression_dataset$speciesId)

#run the pca
pca_res <- prcomp(expression_dataset[5:length(expression_dataset)], scale. = TRUE)

#plot result per species
jpeg(file="/workspace/Biodiversity-Bioinformatics/figures/Species_PCA.jpeg")
autoplot(pca_res, data = expression_dataset, colour = 'speciesId')
dev.off()

#plot results per organ
jpeg(file="/workspace/Biodiversity-Bioinformatics/figures/Organ_PCA.jpeg")
autoplot(pca_res, data = expression_dataset, colour = 'anatEntityName')
dev.off()