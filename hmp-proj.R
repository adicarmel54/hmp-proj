# setwd("/Users/mayalightfoot/Desktop/bst210/project")
#load("RData")

setwd("/Users/adicarmel/Desktop/HARVARD/S1/BST\ 210/microbiome-metabolome-curated-data/data/processed_data/iHMP_IBDMDB_2019")

metadata <- read.table("metadata.tsv", sep = "\t", header = TRUE)
genera <- read.table("genera.tsv", sep = "\t", header = TRUE)
mtb <- read.table("mtb.map.tsv", sep = "\t", header = TRUE)


# EDA
## metadata
dim(metadata)
print(head(metadata))
summary(metadata)

## microbiome
dim(genera)
print(head(genera))
summary(genera)


## metabolome
dim(mtb)
print(head(mtb))
summary

# plot densities
library(ggplot2)
## age
metadata |>
  filter(!is.na(consent_age)) |>
  ggplot(aes(x = consent_age, color = Study.Group), data=metadata) +
  geom_density() +
  ggtitle("Age Distribution by Study Group")

## gender
metadata |> 
  filter(!is.na(Gender)) |>
  ggplot(aes(x = Gender, fill = Study.Group), data=metadata) +
  geom_bar() +
  ggtitle("Gender by Study Group")

## antibiotic use
metadata |> 
  filter(!is.na(Antibiotics)) |>
  ggplot(aes(x = Antibiotics, fill = Study.Group), data=metadata) +
  geom_bar() +
  ggtitle("Antibiotic Use by Study Group")


# Missing data analysis
install.packages("naniar")
library(naniar)
library(dplyr)

vis_miss(metadata) # useful result, 8.2% missing
vis_miss(t(genera)) # doesn't work, too many columns
vis_miss(mtb) # no missing data
missing_proportion <- colMeans(is.na(genera))
zero_proportion <- colMeans(genera == 0)
table(missing_proportion)
table(zero_proportion)

ggplot(data = data.frame(zero_proportion), aes(x = zero_proportion)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "Frequency Distribution of Zeroes in Genera Columns",
       x = "Proportion of Zero Values",
       y = "Frequency (# of columns/microbes)")

sum(zero_proportion == 1)/length(zero_proportion)
sum(zero_proportion > 0.995)/length(zero_proportion)
sum(zero_proportion > 0.99)/length(zero_proportion)
sum(zero_proportion > 0.95)/length(zero_proportion)


# Reduce Dimensionality of Data
## Identify and remove outliers
## Remove strains with zero abundance
## Aggregation 
## Regularization 
## 
## Statistical testing 
## PCA


# Identify potential metadata covariates
## age
## gender
## antibiotic use 


# Select model
## Step-wise
## LASSO/Elastic Net
## maaslin?






# Re-configure Data Frames
library(Maaslin2)

## Genera
genera_data <- as.data.frame(genera)
row.names(genera_data) <- genera_data$Sample  # Replace 'SampleID' with the actual column name
genera_data <- genera_data[, -which(names(genera_data) == "Sample")] # Remove the SampleID column from the data frame (if needed)

## Metabolome 
mtb_data <- as.data.frame(mtb)
row.names(mtb_data) <- mtb_data$Sample  # Replace 'SampleID' with the actual column name
mtb_data <- mtb_data[, -which(names(mtb_data) == "Sample")] # Remove the SampleID column from the data frame (if needed)

## Metadata
metadata_data <- as.data.frame(metadata)
row.names(metadata_data) <- metadata_data$Sample  # Replace 'SampleID' with the actual column name
metadata_data <- metadata_data[, -which(names(metadata_data) == "Sample")] # Remove the SampleID column from the data frame (if needed)



# Maaslin2 for microbiome analysis
## Crohn's disease and ulcerative colitis
maaslin2_micro_CDUC <- Maaslin2(
  input_data = genera_data,
  input_metadata = metadata_data,
  min_prevalence = 0,
  normalization = "NONE",
  fixed_effects = "Study.Group",
  output = "microbiome_output_CDUC",
  reference = c("Study.Group", "nonIBD"))

## IBD
metadata_IBD <- metadata_data # Combine CD and UC outcomes in metadata 
metadata_IBD$diagnosis[metadata_IBD$Study.Group == "UC"] <- "IBD"
metadata_IBD$diagnosis[metadata_IBD$Study.Group == "CD"] <- "IBD"
metadata_IBD$diagnosis[metadata_IBD$Study.Group == "nonIBD"] <- "nonIBD"

maaslin2_micro_IBD <- Maaslin2(
  input_data = genera_data,
  input_metadata = metadata_IBD,
  min_prevalence = 0,
  normalization = "NONE",
  fixed_effects = "diagnosis",
  output = "microbiome_output_IBD",
  reference = c("diagnosis", "nonIBD"))



# Maaslin2 for metabolome analysis
## Crohn's disease and ulcerative Ccolitis
maaslin2_mtb_CDUC <- Maaslin2(
  input_data = mtb_data,
  input_metadata = metadata_data,
  min_prevalence = 0,
  normalization = "NONE",
  fixed_effects = "Study.Group",
  output = "metabolome_output_CDUC",
  reference = c("Study.Group", "nonIBD"))

## IBD
maaslin2_mtb_IBD <- Maaslin2(
  input_data = mtb_data,
  input_metadata = metadata_IBD,
  min_prevalence = 0,
  normalization = "NONE",
  fixed_effects = "diagnosis",
  output = "metabolome_output_IBD",
  reference = c("diagnosis", "nonIBD"))





# Hosmer Lemenshow to test calibrations (more variables than samples)

# HALLA to compare microbiome/metabolome association w outcome


  