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
# Mapping speciesId to species names
species_names <- c("9031" = "Chicken",
                   "9544" = "Macaque", 
                   "9606" = "Human",
                   "9986" = "Rabbit",
                   "10090" = "Mouse",
                   "10116" = "Rat")

# Create a new column with species names
expression_dataset$speciesName <- species_names[expression_dataset$speciesId]
####################
# Plot using species names for coloring
jpeg(file="PCA_post-embryonic_ColoredbySpecies.jpeg")
a <- autoplot(pca_nonemb, 
              data = expression_dataset[which(expression_dataset$devStageName=="post-embryonic stage"),], 
              colour = 'speciesName')
a + scale_color_manual(values = c("black","grey50","orange","forestgreen","red","purple"))
dev.off()
# coloring per organ
jpeg(file="PCA_post-embryonic_ColoredbyOrgan.jpeg")
a<-autoplot(pca_nonemb, data = expression_dataset[which(expression_dataset$devStageName=="post-embryonic stage"),], colour = 'anatEntityName')
a + scale_color_manual(values = c("black","grey50","orange","forestgreen","red","purple","blue"))
dev.off()
####################
# Plot result embryonic stages
# coloring per species
jpeg(file="PCA_embryonic_ColoredbySpecies.jpeg")
a<-autoplot(pca_emb, data = expression_dataset[which(expression_dataset$devStageName=="embryo stage"),], colour = 'speciesName')
a + scale_color_manual(values = c("black","grey50","orange","forestgreen","red","purple"))
dev.off()
# coloring per organ
jpeg(file="PCA_embryonic_ColoredbyOrgan.jpeg")
a<-autoplot(pca_emb, data = expression_dataset[which(expression_dataset$devStageName=="embryo stage"),], colour = 'anatEntityName')
a + scale_color_manual(values = c("black","grey50","orange","forestgreen","red","purple","blue"))
dev.off()
