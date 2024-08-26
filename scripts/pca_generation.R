if (!require('ggfortify', quietly = TRUE)) {
  install.packages("ggfortify", repos = "http://cran.us.r-project.org")
}
library(ggfortify)

####################
# Read data
pca_file_name <- "/workspace/Biodiversity-Bioinformatics/data/pca_dataset.tsv"

expression_dataset <- NULL
# Read tsv or tsv.gz file in case the gzipped file has been uncompressed by the student
if(file.exists(pca_file_name)) {
  expression_dataset = read.table(file = pca_file_name, header = TRUE, sep = "\t")
} else {
  expression_dataset = read.table(file = gzfile(paste0(pca_file_name,".gz")), header = TRUE, sep = "\t")
}
# replace type of speciesId from numeric to character
expression_dataset$speciesId <- as.character(expression_dataset$speciesId)

####################
# Run PCAs
# run the pca for post-embryonic stages
pca_nonemb <- prcomp(expression_dataset[which(expression_dataset$devStageName=="post-embryonic stage"),5:length(expression_dataset)], scale. = TRUE)
#run the pca for embryonic stages
pca_emb <- prcomp(expression_dataset[which(expression_dataset$devStageName=="embryo stage"),5:length(expression_dataset)], scale. = TRUE)

####################
# Plot result post-embryonic stages
# coloring per species
jpeg(file="/workspace/Biodiversity-Bioinformatics/figures/PCA_post-embryonic_ColoredbySpecies.jpeg")
a<-autoplot(pca_nonemb, data = expression_dataset[which(expression_dataset$devStageName=="post-embryonic stage"),], colour = 'speciesId')
a + scale_color_manual(values = c("black","grey50","orange","forestgreen","red","purple")) 
dev.off()

# coloring per organ
jpeg(file="/workspace/Biodiversity-Bioinformatics/figures/PCA_post-embryonic_ColoredbyOrgan.jpeg.jpeg")
a<-autoplot(pca_nonemb, data = expression_dataset[which(expression_dataset$devStageName=="post-embryonic stage"),], colour = 'anatEntityName')
a + scale_color_manual(values = c("black","grey50","orange","forestgreen","red","purple","blue")) 
dev.off()

####################
# Plot result embryonic stages
# coloring per species
jpeg(file="/workspace/Biodiversity-Bioinformatics/figures/PCA_embryonic_ColoredbySpecies.jpeg.jpeg")
a<-autoplot(pca_emb, data = expression_dataset[which(expression_dataset$devStageName=="embryo stage"),], colour = 'speciesId')
a + scale_color_manual(values = c("black","grey50","orange","forestgreen","red","purple")) 
dev.off()

# coloring per organ
jpeg(file="/workspace/Biodiversity-Bioinformatics/figures/PCA_embryonic_ColoredbyOrgan.jpeg.jpeg")
a<-autoplot(pca_emb, data = expression_dataset[which(expression_dataset$devStageName=="embryo stage"),], colour = 'anatEntityName')
a + scale_color_manual(values = c("black","grey50","orange","forestgreen","red","purple","blue")) 
dev.off()
