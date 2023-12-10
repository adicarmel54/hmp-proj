
model <- function(log_x, log_y) {
  
  log_ridge.mod = glmnet(log_x, log_y, alpha = 0, family = "binomial")
  log_cv.ridge <- cv.glmnet(log_x, log_y, alpha = 0, family = "binomial")
  log_lambda_min_ridge <- log_cv.ridge$lambda.min
  log_lambda_1se_ridge <- log_cv.ridge$lambda.1se
  
  log_lasso.mod = glmnet(log_x, log_y, alpha = 1, family = "binomial")
  log_cv.lasso <- cv.glmnet(log_x, log_y, alpha = 1, family = "binomial")
  log_lambda_min_lasso <- log_cv.lasso$lambda.min
  log_lambda_1se_lasso <- log_cv.lasso$lambda.1se
  
  log_EN = glmnet(log_x, log_y, alpha = 0.5, family = "binomial")
  log_cv.EN <- cv.glmnet(log_x, log_y, alpha = 0.5, family = "binomial")
  log_lambda_min_EN <- log_cv.lasso$lambda.min
  log_lambda_1se_EN <- log_cv.lasso$lambda.1se
  
  log_out <- cbind(coef(log_cv.ridge, s = log_lambda_1se_ridge),
                   coef(log_cv.lasso, s = log_lambda_1se_lasso),
                   coef(log_cv.EN, s = log_lambda_1se_EN))
  colnames(log_out) <- c("Ridge", "LASSO", "EN")
  
  return(log_out)
}


# phylum 

phylum_bin <- cbind(phylum, Y_bin)
phylum_bin <- phylum_bin[-1]

log_full <- glm(Y_bin ~ ., data = phylum_bin)

log_x <- as.matrix(phylum_bin[1:112])
log_y <- phylum_bin[[113]]

model(log_x, log_y)
print(log_out)



# order 

order_bin <- cbind(order, Y_bin)
order_bin <- order_bin[-1]

log_full <- glm(Y_bin ~ ., data = order_bin)

log_x <- as.matrix(order_bin[1:825])
log_y <- order_bin[[826]]

model(log_x, log_y)
print(log_out)



# PCA_50

# from github
numeric_data <- genera_non_zero[, sapply(genera_non_zero, is.numeric)]
pca_result <- prcomp(numeric_data, scale = TRUE)
summary(pca_result)
plot(pca_result)
pc_scores <- pca_result$x
top_50_pc_scores <- pc_scores[, 1:50]
top_50_pc_dataframe <- as.data.frame(top_50_pc_scores)

PCA_50 <- top_50_pc_dataframe

PCA_50_bin <- cbind(PCA_50, Y_bin)

log_full <- glm(Y_bin ~ ., data = PCA_50_bin)

log_x <- as.matrix(PCA_50_bin[1:50])
log_y <- PCA_50_bin[[51]]

model(log_x, log_y)
print(log_out)



# PCA_100

PCA_100_bin <- cbind(PCA_100, Y_bin)
PCA_100_bin <- PCA_100_bin[-1]

log_full <- glm(Y_bin ~ ., data = PCA_100_bin)

log_x <- as.matrix(PCA_100_bin[1:100])
log_y <- PCA_100_bin[[101]]

model(log_x, log_y)
print(log_out)



#encode_50

encode_50_bin <- cbind(encode_50, Y_bin)
encode_50_bin <- encode_50_bin[-1]

log_full <- glm(Y_bin ~ ., data = encode_100_bin)

log_x <- as.matrix(encode_50_bin[1:50])
log_y <- encode_50_bin[[51]]

model(log_x, log_y)
print(log_out)



#encode_100

encode_100_bin <- cbind(encode_100, Y_bin)
encode_100_bin <- encode_100_bin[-1]

log_full <- glm(Y_bin ~ ., data = encode_100_bin)

log_x <- as.matrix(encode_100_bin[1:100])
log_y <- encode_100_bin[[101]]

model(log_x, log_y)
print(log_out)
