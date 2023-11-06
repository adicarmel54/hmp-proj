library(MASS)
library(ape)
library(dplyr)
library(tidyverse)
library(stringr)

# Load data
setwd("/Users/adicarmel/Desktop/HARVARD/S1/BST\ 210/microbiome-metabolome-curated-data/data/processed_data/iHMP_IBDMDB_2019")
metadata <- read.table("metadata.tsv", sep = "\t", header = TRUE)
genera <- read.table("genera.tsv", sep = "\t", header = TRUE)
mtb <- read.table("mtb.map.tsv", sep = "\t", header = TRUE)
dim(mtb)
# Define outcomes
Y <- metadata$Study.Group
Y_bin <- as.integer(Y %in% c("UC", "CD"))              # 1 for IBD (UC or CD), 0 otherwise
Y_cat <- ifelse(Y == "UC", 1, ifelse(Y == "CD", 2, 0)) # 1 for UC, 2 for CD, 0 otherwise
# Remove outcome from metadata
metadata <- metadata[, !colnames(metadata) %in% "Study.Group"]

# Filter metadata to only include columns without missing data
metadata_complete_cols <- metadata[ , colSums(is.na(metadata)) == 0]
# Remove non-features from metadata
metadata_feats <- metadata_complete_cols[, c("Sample", "Subject", "Gender", "week_num", "interval_days", "visit_num", "Antibiotics", "race")]
metadata_feats$Antibiotics <- as.numeric(metadata_feats$Antibiotics == "Yes")
metadata_feats$Gender <- as.numeric(metadata_feats$Gender == "Female")
metadata_feats <- cbind(metadata_feats, model.matrix(~ race - 1, data = metadata_feats))
metadata_feats <- metadata_feats[, -which(names(metadata_feats) %in% c("race"))]

# Remove genera columns that contain all 0's
genera_non_zero <- genera[, !apply(genera, 2, function(x) all(x == 0))]

# PCA -- Just genera
numeric_data <- genera_non_zero[, sapply(genera_non_zero, is.numeric)]
pca_result <- prcomp(numeric_data, scale = TRUE)
summary(pca_result)
plot(pca_result)

barplot(pca_result$sdev^2 / sum(pca_result$sdev^2), 
        names.arg = 1:length(pca_result$sdev), 
        xlab = "Principal Component", 
        ylab = "Proportion of Variance Explained",
        main = "Variance Explained by Principal Components")

pc_scores <- pca_result$x
top_50_pc_scores <- pc_scores[, 1:50]
top_50_pc_dataframe <- as.data.frame(top_50_pc_scores)

eigenvalues <- pca_result$sdev^2
variance_explained <- sum(eigenvalues[1:50])
total_variance <- sum(eigenvalues)
proportion_explained <- variance_explained / total_variance
print(paste("Proportion of variance explained by top 50 principal components:", proportion_explained))

cumulative_var <- cumsum(eigenvalues) / sum(eigenvalues)
num_pcs <- length(eigenvalues)
pcs <- 1:num_pcs
plot(pcs, cumulative_var, type = "b", xlab = "Number of Principal Components", ylab = "Cumulative Proportion of Variance Explained", main = "Scree Plot")

# So, the useful things from above are:
# Outcomes: Y_bin (binary), Y_cat (categorical)
# metadata_feats: the metadata features we may use
# genera_non_zero: genera where zero-valued columns are filtered out
# top_50_pc_dataframe: the vectors of the top 50 principal components for the genera data





# -------- THIS DOESNT YET WORK ------
# Combining by taxonomy 

taxonomic_levels <- strsplit(names(genera_non_zero[, !colnames(genera_non_zero) %in% c("Sample")]), "\\.")
taxonomic_levels <- lapply(taxonomic_levels, function(x) x[1:6])
taxonomic_df <- do.call(rbind, taxonomic_levels)
colnames(taxonomic_df) <- c("kingdom", "phylum", "class", "order", "genus", "species")
dim(taxonomic_df)
length(colnames(genera_non_zero[, !colnames(genera_non_zero) %in% c("Sample")]))
taxonomic_df$full <- colnames(genera_non_zero[, !colnames(genera_non_zero) %in% c("Sample")])
genera_taxonomy <- cbind(genera_non_zero, taxonomic_df)

# Merge genera_non_zero with taxonomic_levels_df by the new Taxonomic_Level column
merged_data <- merge(genera_non_zero, taxonomic_levels_df, by="Taxonomic_Level")

merged_data <- merge(genera_non_zero, taxonomic_df, by="Sample")

# Sum counts by kingdom
kingdom_df <- merged_data %>%
  group_by(kingdom) %>%
  summarise_all(sum) %>%
  select(-phylum, -class, -order, -genus, -species)  # Remove unnecessary columns

# Sum counts by phylum
phylum_df <- merged_data %>%
  group_by(kingdom, phylum) %>%
  summarise_all(sum) %>%
  select(-class, -order, -genus, -species)  # Remove unnecessary columns

# Sum counts by class
class_df <- merged_data %>%
  group_by(kingdom, phylum, class) %>%
  summarise_all(sum) %>%
  select(-order, -genus, -species)  # Remove unnecessary columns

numeric_taxonomic_df <- taxonomic_df[, sapply(taxonomic_df, is.numeric)]

colnames(taxonomic_df) <- c("kingdom", "phylum", "class", "order", "genus", "species")
install.packages("ape")  # Install the ape package if you haven't already

# Create a phylogenetic tree using taxonomic levels
tree <- hclust(dist(as.matrix(numeric_taxonomic_df)), method = "complete")
phylo_tree <- as.phylo(tree)

# Plot the phylogenetic tree
plot(phylo_tree, cex = 0.5, label.offset = 0.5, main = "Taxonomic Tree Visualization")


#LDA
merged_data <- merge(genera, metadata, by = "Sample")
merged_data$class <- factor(merged_data$Study.Group, levels = unique(merged_data$Study.Group))
merged_data$class <- factor(as.numeric(merged_data$class), levels = c(1, 2, 3))
lda_model <- lda(class ~ ., data = merged_data)
print(lda_model)