

setwd("/Users/mayalightfoot/Desktop/bst210/project")
load("RData")

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
## age
metadata |> 
  filter(!is.na(consent_age)) |>
  ggplot(aes(x = consent_age, color = Study.Group)) +
  geom_density() +
  ggtitle("Age Distribution by Study Group")

## gender
metadata |> 
  filter(!is.na(Gender)) |>
  ggplot(aes(x = Gender, fill = Study.Group)) +
  geom_bar() +
  ggtitle("Gender by Study Group")

## antibiotic use
metadata |> 
  filter(!is.na(Antibiotics)) |>
  ggplot(aes(x = Antibiotics, fill = Study.Group)) +
  geom_bar() +
  ggtitle("Antibiotic Use by Study Group")





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


  