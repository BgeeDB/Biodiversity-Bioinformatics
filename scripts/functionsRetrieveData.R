retrieveOrthologsExpression <- function(bgeeRelease = "current", downloadPath = getwd(),
                                        experimentId = NULL, referenceSpeciesId = NULL, speciesIds = NULL, anatEntityIds = NULL, removeDownloadFiles = FALSE,
                                        onlyOneToOneOrthologs = FALSE, onlyOneToManyOrthologs = FALSE, createFile = FALSE, orthologs = NULL) {
  # we consider that BgeeDB package is already installed... lazy

  # first retrieve ortholog genes for the specified species if not provided
  if (is.null(orthologs)) {
    orthologs <- listOrthologs(mandatorySpecies = speciesIds, referenceSpecies = referenceSpeciesId, onlyOneToOne = onlyOneToOneOrthologs,
                               onlyOneToMany = onlyOneToManyOrthologs, bgeeRelease = bgeeRelease)
  }

  # ugly quick and dirty implementation
  gene_families <- unique(orthologs[1])

  # now retrieve expression data for all libraries of one species
  all_expression <- NULL
  allSpeciesIds <- c(referenceSpeciesId, speciesIds)
  for (speciesId in allSpeciesIds) {
    # if 1-to-N gene family have to be considered we need to retrieve all genes corresponding to this
    # species. They are in different columns in the orthologs data.frame. We first merge those columns in one
    # single list called species_orthologous_genes
    species_orthologous_genes <- as.data.frame(matrix(nrow = 0, ncol = 2, dimnames = list(NULL, c("geneId", "geneFamily"))))
    for (i in grep(paste0(speciesId, "|", speciesId, "_\\d+"), as.character(colnames(orthologs)))) {
      column_orthologous_genes <- as.data.frame(list(orthologs[,i], seq_len(nrow(orthologs))))
      colnames(column_orthologous_genes) <- c("geneId", "geneFamily")
      species_orthologous_genes <- rbind(species_orthologous_genes, column_orthologous_genes)
    }
    bgee_object <- Bgee$new(species = as.character(speciesId), dataType = "rna_seq")
    library_metadata <- getAnnotation(bgee_object)[[1]]
    #only keep metadata from provided experiment
    if (! is.null(experimentId)) {
      library_metadata <- library_metadata[library_metadata$Experiment.ID == experimentId,]
    }
    # now filter on anatomical entities
    if (!is.null(anatEntityIds)) {
      library_metadata <- library_metadata[library_metadata$Anatomical.entity.ID %in% anatEntityIds,]
    }
    if (nrow(library_metadata) == 0) {
      stop("No library available for species ", speciesId, ". Please remove that species or update the experimentId, ",
           "anatEntityIds or bgeeRelease")
    }
    expression_data <- getSampleProcessedData(sampleId = unique(library_metadata$Library.ID), myBgeeObject = bgee_object)
    expression_data$Species.ID <- speciesId
    orthologs_expression_data <- expression_data[expression_data$Gene.ID %in% species_orthologous_genes$geneId,]
    #implementation with lapply (slow so commented it)
    # which_f <- function(x, orthologs, speciesId, gene_families) {return(which(gene_families == orthologs[which(orthologs[,as.character(speciesId)] == x),1]))}
    # orthologs_expression_data$Gene.Family <- lapply(X = orthologs_expression_data$Gene.ID, FUN = which_f, orthologs = orthologs, speciesId = speciesId, gene_families = gene_families)
    # naive implementation with for loop... faster than using lapply
    orthologs_expression_data$Gene.Family <- NA
    for (row_number in seq_len(nrow(orthologs_expression_data))) {
      orthologs_expression_data$Gene.Family[row_number] <- species_orthologous_genes$geneFamily[species_orthologous_genes$geneId == orthologs_expression_data$Gene.ID[row_number]]

    }
    all_expression <- rbind(all_expression, orthologs_expression_data)
  }

  if (createFile) {
    output_file <- file.path(downloadPath, paste0("expression_", paste(allSpeciesIds, collapse = "_"), ".tsv"))
    message("writing expression data in the file ", output_file)
    write.table(x = all_expression, file = output_file,
                sep = "\t", quote = F, row.names = F, col.names = T)
  } else {
    message("properly retrieved Bgee expression data")
    return(all_expression)
  }

}

listOrthologs <- function(bgeeRelease = "current", downloadPath = getwd(),
                          mandatorySpecies = NULL, optionalSpecies = NULL, referenceSpecies = NULL, onlyOneToOne = FALSE,
                          onlyOneToMany = FALSE, removeDownloadedFiles = FALSE, createFile = FALSE) {

  # TODO check all arguments combination and throw error if not an expected one (e.g both onlyOneToOne and onlyOneToMany)
  if (!is.null(referenceSpecies) && onlyOneToOne) {
    stop("please provide a reference species only when 1-to-many orthologs have to be retrieved")
  }
  if (onlyOneToMany && onlyOneToOne) {
    stop("please select only one argument out of onlyOneToOne and onlyOneToMany")
  }
  if (is.null(referenceSpecies) && !onlyOneToOne) {
    stop("please provide a reference species if 1-to-many orthologs have to be retrieved")
  }

  # first draft of the implementation lots of future functionalities to implement....
  if (!is.null(optionalSpecies)) {stop("optionalSpecies not yet implemented")}

  ftp_url <- "https://www.bgee.org/ftp/RELEASE_VERSION/homologous_genes/OMA_orthologs.zip"
  orthologs_dir <- "bgeeOrthologs"
  archive_file <- basename(ftp_url)
  download_dir <- file.path(downloadPath, orthologs_dir)
  if (!dir.exists(download_dir)) {
    message("dir does not exist")
    dir.create(download_dir)
    downloaded_archive <- file.path(download_dir, archive_file)
    ftp_url <- gsub(pattern = "RELEASE_VERSION", replacement = bgeeRelease, x = ftp_url)
    # for now we only provide one zip archive containing all orthologs files
    file_path <- download.file(url = ftp_url, destfile = downloaded_archive)
    unzip_files <- unzip(zipfile = downloaded_archive, exdir = download_dir)
  }
  unzipped_files <- list.files(path = download_dir, include.dirs = FALSE, full.names = TRUE, recursive = TRUE,
                               pattern = "orthologs_.*.csv")

  # now files are downloaded and orthologs retrieval can starts
  #for now we consider retrieval of all orthologs from mandatory species
  retrieved_orthologs_files <- NULL
  retrieved_species <- NULL
  all_orthologs <- NULL
  for (file_path in unzipped_files) {
    speciesIds <- unlist(regmatches(x = basename(file_path), m = gregexpr(pattern = "[0-9]+", basename(file_path))))
    if (!is.null(referenceSpecies) && referenceSpecies %in% as.numeric(speciesIds) &&
        (as.numeric(speciesIds[1]) %in% mandatorySpecies || as.numeric(speciesIds[2]) %in% mandatorySpecies)) {
      orthologs <- read.table(file = file_path, header = TRUE, sep = ",", quote = "")
      if (speciesIds[1] == as.character(referenceSpecies)) {
        orthologs <- orthologs[,c(1,2)]
        colnames(orthologs) <- c(referenceSpecies, speciesIds[2])
      } else {
        orthologs <- orthologs[,c(2,1)]
        colnames(orthologs) <- c(referenceSpecies, speciesIds[1])
      }
      duplicated_species2 <- unique(orthologs[duplicated(orthologs[2]),][2])
      orthologs <- orthologs[! orthologs[,2] %in% duplicated_species2[[1]],]
      if (onlyOneToMany) {
        duplicated_species1 <- unique(orthologs[duplicated(orthologs[1]),][1])
        orthologs <- orthologs[ orthologs[,1] %in% duplicated_species1[[1]],]
      }
      if (as.character(referenceSpecies) %in% retrieved_species) {
        all_orthologs <- merge(x = all_orthologs, y = orthologs, by = as.character(referenceSpecies))
      } else {
        all_orthologs <- orthologs
      }
      retrieved_species <- unique(c(retrieved_species, speciesIds))
    } else if (is.null(referenceSpecies) && as.numeric(speciesIds[1]) %in% mandatorySpecies &&
               as.numeric(speciesIds[2] %in% mandatorySpecies)) {
      retrieved_orthologs_files <- c(retrieved_orthologs_files, file_path)
      orthologs <- read.table(file = file_path, header = TRUE, sep = ",", quote = "")[,1:2]
      # only one-to-one requested so we remove all duplicated genes
      duplicated_species2 <- unique(orthologs[duplicated(orthologs[2]),][2])
      duplicated_species1 <- unique(orthologs[duplicated(orthologs[1]),][1])
      orthologs <- orthologs[! orthologs[,1] %in% duplicated_species1[[1]],]
      orthologs <- orthologs[! orthologs[,2] %in% duplicated_species2[[1]],]
      colnames(orthologs) <- c(speciesIds[1], speciesIds[2])
      if ((speciesIds[1] %in% retrieved_species) & (speciesIds[2] %in% retrieved_species)) {
        all_orthologs <- merge(x = all_orthologs, y = orthologs, by = speciesIds)
      } else if (speciesIds[1] %in% retrieved_species) {
        all_orthologs <- merge(x = all_orthologs, y = orthologs, by = speciesIds[1])
      } else if (speciesIds[2] %in% retrieved_species) {
        all_orthologs <- merge(x = all_orthologs, y = orthologs, by = speciesIds[2])
      } else {
        all_orthologs <- orthologs
      }
      retrieved_species <- unique(c(retrieved_species, speciesIds))
    }
  }
  if (removeDownloadedFiles) {
    file.remove(download_dir, recursive = TRUE)
  }
  if (createFile) {
    output_file <- file.path(downloadPath, paste0("orthologs_", paste(retrieved_species, collapse = "_"), ".tsv"))
    message("writing orthologous genes in the file ", output_file)
    write.table(x = all_orthologs, file = output_file,
                sep = "\t", quote = F, row.names = F, col.names = T)
  } else {
    message("properly retrieved Bgee orthologs")
    return(all_orthologs)
  }
}

loadExpressionKaessmann2019 <- function() {
  expression_9031 <- retrieveOrthologsExpression(experimentId = "ERP108704", speciesIds = 9031, orthologs = orthologs_one_to_one)
  expression_9986 <- retrieveOrthologsExpression(experimentId = "ERP108867", speciesIds = 9986, orthologs = orthologs_one_to_one)
  expression_10090 <- retrieveOrthologsExpression(experimentId = "ERP108893", speciesIds = 10090, orthologs = orthologs_one_to_one)
  expression_10116 <- retrieveOrthologsExpression(experimentId = "ERP108920", speciesIds = 10116, orthologs = orthologs_one_to_one)
  expression_9544 <- retrieveOrthologsExpression(experimentId = "ERP108988", speciesIds = 9544, orthologs = orthologs_one_to_one)
  expression_9606 <- retrieveOrthologsExpression(experimentId = "ERP109002", speciesIds = 9606, orthologs = orthologs_one_to_one)
  expression_13616 <- retrieveOrthologsExpression(experimentId = "ERP109071", speciesIds = 13616, orthologs = orthologs_one_to_one)
  # And then merge the results
  expression_kaessmann_2019_bgee <- rbind(expression_9031, expression_9986, expression_10090,
                                          expression_10116, expression_9544, expression_9606, expression_13616)
  return(expression_kaessmann_2019_bgee)
}

transfrom_pca <- function(orthologs_one_to_one, processed_expression_data) {
  kaessmann_remapped <- processed_expression_data
  all_exression_transformed <- as.data.frame(matrix(nrow = 0, ncol = nrow(orthologs_one_to_one) + 4))
  for (speciesId in unique(kaessmann_remapped$Species.ID)) {
    for (libraryId in unique(kaessmann_remapped$Library.ID[kaessmann_remapped$Species.ID == speciesId])) {
      
      transformed <- as.data.frame(t(kaessmann_remapped[kaessmann_remapped$Species.ID == speciesId & kaessmann_remapped$Library.ID== libraryId, c("TPM", "Gene.Family")]))
      colnames(transformed) <- as.integer(transformed["Gene.Family",])
      transformed <- transformed[!row.names(transformed) %in% "Gene.Family",]
      transformed$libraryId <- libraryId
      transformed$speciesId <- speciesId
      transformed$anatEntityName <- unique(kaessmann_remapped$Anatomical.entity.name[kaessmann_remapped$Species.ID == speciesId
                                                                                     & kaessmann_remapped$Library.ID== libraryId])
      transformed$devStageName <- unique(kaessmann_remapped$Stage.name[kaessmann_remapped$Species.ID == speciesId
                                                                       & kaessmann_remapped$Library.ID== libraryId])
      if(nrow(all_exression_transformed) == 0) {
        all_exression_transformed <- transformed
      } else {
        all_exression_transformed <- dplyr::bind_rows(all_exression_transformed, transformed)
      }
    }
  }
  # reorder columns of transformed data
  all_exression_transformed <- all_exression_transformed[, c("libraryId", "speciesId", "anatEntityName", "devStageName", seq_len(nrow(orthologs_one_to_one)))]
  
  # remove quotes from the anat. etntiy and dev. stage names
  all_exression_transformed$anatEntityName <- gsub('^.|.$', '', all_exression_transformed$anatEntityName)
  all_exression_transformed$devStageName <- gsub('^.|.$', '', all_exression_transformed$devStageName)
  
  return(all_exression_transformed)
}

remapAnatEntities <- function(dataframe) {
  
  kidney_ids <- c("UBERON:0000080", "UBERON:0000082", "UBERON:0002113")
  kidney_names <- c("\"mesonephros\"", "\"adult mammalian kidney\"", "\"kidney\"")
  brain_ids <- c("UBERON:0001890", "UBERON:0002028", "UBERON:0000955")
  brain_names <- c("\"forebrain\"", "\"hindbrain\"", "\"brain\"")
  
  dataframe$Anatomical.entity.ID[dataframe$Anatomical.entity.ID %in% kidney_ids] <- "UBERON:0002113"
  dataframe$Anatomical.entity.ID[dataframe$Anatomical.entity.ID %in% brain_ids] <- "UBERON:0000955"
  dataframe$Anatomical.entity.name[dataframe$Anatomical.entity.name %in% kidney_names] <- "\"kidney\""
  dataframe$Anatomical.entity.name[dataframe$Anatomical.entity.name %in% brain_names] <- "\"brain\""
  
  return(dataframe)
}

remapDevStages <- function(dataframe) {
  #9031
  bgee <- Bgee$new(species = "9031", dataType = "rna_seq")
  embryo_9031 <-getSampleProcessedData(myBgeeObject = bgee, experimentId = "ERP108704", stageId = "UBERON:0000068",withDescendantStages = TRUE)
  post_embryo_9031 <- getSampleProcessedData(myBgeeObject = bgee, experimentId = "ERP108704", stageId = "UBERON:0000092",withDescendantStages = TRUE)
  embryo <- unique(embryo_9031$Stage.ID)
  post_embryo <- unique(post_embryo_9031$Stage.ID)
  #9986
  bgee <- Bgee$new(species = "9986", dataType = "rna_seq")
  embryo_9986 <-getSampleProcessedData(myBgeeObject = bgee, experimentId = "ERP108867", stageId = "UBERON:0000068",withDescendantStages = TRUE)
  post_embryo_9986 <- getSampleProcessedData(myBgeeObject = bgee, experimentId = "ERP108867", stageId = "UBERON:0000092",withDescendantStages = TRUE)
  embryo <- c(embryo, unique(embryo_9986$Stage.ID))
  post_embryo <- c(post_embryo, unique(post_embryo_9986$Stage.ID))
  #10090
  bgee <- Bgee$new(species = "10090", dataType = "rna_seq")
  embryo_10090 <-getSampleProcessedData(myBgeeObject = bgee, experimentId = "ERP108893", stageId = "UBERON:0000068",withDescendantStages = TRUE)
  post_embryo_10090 <- getSampleProcessedData(myBgeeObject = bgee, experimentId = "ERP108893", stageId = "UBERON:0000092",withDescendantStages = TRUE)
  embryo <- c(embryo, unique(embryo_10090$Stage.ID))
  post_embryo <- c(post_embryo, unique(post_embryo_10090$Stage.ID))
  #10116
  bgee <- Bgee$new(species = "10116", dataType = "rna_seq")
  embryo_10116 <-getSampleProcessedData(myBgeeObject = bgee, experimentId = "ERP108920", stageId = "UBERON:0000068",withDescendantStages = TRUE)
  post_embryo_10116 <- getSampleProcessedData(myBgeeObject = bgee, experimentId = "ERP108920", stageId = "UBERON:0000092",withDescendantStages = TRUE)
  embryo <- c(embryo, unique(embryo_10116$Stage.ID))
  post_embryo <- c(post_embryo, unique(post_embryo_10116$Stage.ID))
  #9544
  bgee <- Bgee$new(species = "9544", dataType = "rna_seq")
  embryo_9544 <-getSampleProcessedData(myBgeeObject = bgee, experimentId = "ERP108988", stageId = "UBERON:0000068",withDescendantStages = TRUE)
  post_embryo_9544 <- getSampleProcessedData(myBgeeObject = bgee, experimentId = "ERP108988", stageId = "UBERON:0000092",withDescendantStages = TRUE)
  embryo <- c(embryo, unique(embryo_9544$Stage.ID))
  post_embryo <- c(post_embryo, unique(post_embryo_9544$Stage.ID))
  #9606
  bgee <- Bgee$new(species = "9606", dataType = "rna_seq")
  embryo_9606 <-getSampleProcessedData(myBgeeObject = bgee, experimentId = "ERP109002", stageId = "UBERON:0000068",withDescendantStages = TRUE)
  post_embryo_9606 <- getSampleProcessedData(myBgeeObject = bgee, experimentId = "ERP109002", stageId = "UBERON:0000092",withDescendantStages = TRUE)
  embryo <- c(embryo, unique(embryo_9606$Stage.ID))
  post_embryo <- c(post_embryo, unique(post_embryo_9606$Stage.ID))
  #13616
  bgee <- Bgee$new(species = "13616", dataType = "rna_seq")
  embryo_13616 <-getSampleProcessedData(myBgeeObject = bgee, experimentId = "ERP109071", stageId = "UBERON:0000068",withDescendantStages = TRUE)
  post_embryo_13616 <- getSampleProcessedData(myBgeeObject = bgee, experimentId = "ERP109071", stageId = "UBERON:0000092",withDescendantStages = TRUE)
  embryo <- c(embryo, unique(embryo_13616$Stage.ID))
  post_embryo <- c(post_embryo, unique(post_embryo_13616$Stage.ID))
  
  dataframe$Stage.ID[dataframe$Stage.ID %in% embryo] <- "UBERON:0000068"
  dataframe$Stage.ID[dataframe$Stage.ID %in% post_embryo] <- "UBERON:0000092"
  dataframe$Stage.name[dataframe$Stage.ID == "UBERON:0000092"] <- "\"post-embryonic stage\""
  dataframe$Stage.name[dataframe$Stage.ID == "UBERON:0000068"] <- "\"embryo stage\""
  
  return(dataframe)
}