
setwd("/Users/mayalightfoot/Desktop/bst210/project")
load("RData")



# Retreive these data:
# Y_bin (binary), Y_cat (categorical)
# metadata_feats: the metadata features we may use
# genera_non_zero: genera where zero-valued columns are filtered out
# top_50_pc_dataframe: the vectors of the top 50 principal components for the genera data


# Convert outcome vectors into data frames

Y <- data.frame(Y)
Y_bin <- data.frame(Y_bin)
Y_cat <- data.frame(Y_cat)

Y_bin_data <- cbind(top_50_pc_dataframe, Y_bin)
Y_cat_data <- cbind(top_50_pc_dataframe, Y_cat)



## Linear Regression

linear_mod_full <- lm(Y_bin ~ ., data = Y_bin_data)
linear_mod_reduced <- lm(Y_bin ~ 1, data = Y_bin_data)
summary(linear_mod_full)
summary(linear_mod_reduced)

step_mod <- step(linear_mod_full, direction = c("both"))

step_mod_forward <- step(linear_mod_reduced, direction = "forward", scope = (~ .), data = Y_bin_data)

step_mod_backward <- step(linear_mod_full, direction = "backward")



## Logistic Regression

glm_bin_full <- glm(Y_bin ~ ., data = Y_bin_data)
summary(glm_bin_full)


x <- as.matrix(Y_bin_data[1:50])
y <- Y_bin_data[[51]]

ridge.mod = glmnet(x, y, alpha = 0, family = "binomial")
cv.ridge <- cv.glmnet(x, y, alpha = 0, family = "binomial")
lambda_min_ridge <- cv.ridge$lambda.min
lambda_1se_ridge <- cv.ridge$lambda.1se
coef(cv.ridge, s = lambda_min_ridge)

lasso.mod = glmnet(x, y, alpha = 1, family = "binomial")
cv.lasso <- cv.glmnet(x, y, alpha = 1, family = "binomial")
lambda_min_lasso <- cv.lasso$lambda.min
lambda_1se_lasso <- cv.lasso$lambda.1se
coef(cv.lasso, s = lambda_min_lasso)

elasticnet = glmnet(x, y, alpha = 0.5, family = "binomial")
cv.EN <- cv.glmnet(x, y, alpha = 0.5, family = "binomial")
lambda_min_EN <- cv.lasso$lambda.min
lambda_1se_EN <- cv.lasso$lambda.1se
coef(cv.EN, s = lambda_1se_EN)


out <- cbind(coef(cv.ridge, s = lambda_min_ridge),
             coef(cv.lasso, s = lambda_min_lasso),
             coef(cv.EN, s = lambda_min_EN))
colnames(out) <- c("Ridge", "LASSO", "EN")
out



## Multinomial Regression

glm_cat_full <- glm(Y_cat ~ ., data = Y_cat_data)
summary(glm_cat_full)


x <- as.matrix(Y_cat_data[1:50])
y <- Y_cat_data[[51]]

ridge.mod = glmnet(x, y, alpha = 0, family = "multinomial")
cv.ridge <- cv.glmnet(x, y, alpha = 0, family = "multinomial")
lambda_min_ridge <- cv.ridge$lambda.min
lambda_1se_ridge <- cv.ridge$lambda.1se
coef(cv.ridge, s = lambda_min_ridge)

lasso.mod = glmnet(x, y, alpha = 1, family = "multinomial")
cv.lasso <- cv.glmnet(x, y, alpha = 1, family = "multinomial")
lambda_min_lasso <- cv.lasso$lambda.min
lambda_1se_lasso <- cv.lasso$lambda.1se
coef(cv.lasso, s = lambda_min_lasso)

elasticnet = glmnet(x, y, alpha = 0.5, family = "multinomial")
cv.EN <- cv.glmnet(x, y, alpha = 0.5, family = "multinomial")
lambda_min_EN <- cv.lasso$lambda.min
lambda_1se_EN <- cv.lasso$lambda.1se
coef(cv.EN, s = lambda_1se_EN)




## Poisson Regression
## Binomial

poisson_bin_full <- glm(Y_bin ~ ., data = Y_bin_data, family = poisson)
summary(poisson_bin_full)

## Categorical

poisson_cat_full <- glm(Y_cat ~ ., data = Y_cat_data, family = poisson)
summary(poisson_cat_full)











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