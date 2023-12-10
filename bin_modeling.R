
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


# metadata

metadata_feats_bin <- cbind(metadata_feats, Y_bin)
metadata_feats_bin <- metadata_feats_bin[-c(1,2)]

log_full <- glm(Y_bin ~ ., data = metadata_feats_bin)

log_x <- as.matrix(metadata_feats_bin[1:9])
log_y <- metadata_feats_bin[[10]]

model(log_x, log_y)
summary(log_full)

log_EN <- glm(Y_bin ~ Gender + Antibiotics + race_Black.or.African.American + race_Other + race_White, data = metadata_feats_bin)
summary(log_EN)

# phylum 

phylum_bin <- cbind(phylum, Y_bin)
phylum_bin <- phylum_bin[-1]

log_full <- glm(Y_bin ~ ., data = phylum_bin)

log_x <- as.matrix(phylum_bin[1:112])
log_y <- phylum_bin[[113]]

model(log_x, log_y)
summary(log_full)

log_EN <- glm(Y_bin ~ p__Aerophobota + p__Aquificota + p__CSSED10.310 + p__Caldisericota + p__Calditrichota + p__Coprothermobacterota + p__Dependentiae + p__Desulfobacterota_G + p__JACPQY01 + p__Krumholzibacteriota + p__Margulisbacteria + p__Mcinerneyibacteriota + p__Nitrospirota + p__QNDG01 + p__SAR324 + p__UBA10199 + p__Verrucomicrobiota, data = phylum_bin)
summary(log_EN)




# order 

order_bin <- cbind(order, Y_bin)
order_bin <- order_bin[-1]

log_full <- glm(Y_bin ~ ., data = order_bin)

log_x <- as.matrix(order_bin[1:825])
log_y <- order_bin[[826]]

model(log_x, log_y)
summary(log_full)

log_EN <- glm(Y_bin ~ o__Aerophobales +
              o__B22.G15 +
              o__Babeliales +   
              o__Bacillales_D +
              o__CACIAM.22H2 +
              o__CAIYCZ01 +
              o__CPR2 +  
              o__CSSED10.310 +
              o__Calditrichales +
              o__Carboxydothermales +
              o__Ch66 +
              o__DTPP01 +
              o__DTU022 +
              o__Daviesbacterales +
              o__Elsterales +
              o__Fen.727 +
              o__Fibrobacterales +
              o__GCA.002796325 +
              o__GCA.2747355 +
              o__Ga0077536 +
              o__Hydrogenothermales +
              o__Immundisolibacterales +
              o__JABMRG01 +
              o__JACPQY01 +
              o__JAHIPI01 +
              o__JGIOTU.2 +
              o__Krumholzibacteriales +
              o__ML615J.28 +
              o__Minwuiales +
              o__O2.12.FULL.45.9 +
              o__Oligoflexales +
              o__Reyranellales +
              o__SM23.28.2 +
              o__SMXP01 +
              o__SP197 +
              o__SS1.B.03.39 +
              o__SZUA.146 +
              o__Silvanigrellales +
              o__Syntrophomonadales +
              o__Thiohalorhabdales +
              o__UBA10575 +
              o__UBA1280 +
              o__UBA1400 +
              o__UBA2968 +
              o__UBA6002 +
              o__UBA6077 +
              o__UBA7916 +
              o__UBA796 +
              o__UBA9942 +
              o__UBA9983_A +
              o__Vicinamibacterales +
              o__XYD2.FULL.50.16
              , data = order_bin)
summary(log_EN)



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
summary(log_full)

log_EN <- glm(Y_bin ~ PC1 + PC6 + PC7 + PC8 + PC11 + PC13 + PC14 + PC18 + PC21 + PC22 + PC30 + PC31 + PC34 + PC36 + PC38 + PC42 + PC43 + PC49, data = PCA_50_bin)
summary(log_EN)



# PCA_100

PCA_100_bin <- cbind(PCA_100, Y_bin)
PCA_100_bin <- PCA_100_bin[-1]

log_full <- glm(Y_bin ~ ., data = PCA_100_bin)

log_x <- as.matrix(PCA_100_bin[1:100])
log_y <- PCA_100_bin[[101]]

model(log_x, log_y)
summary(log_full)

log_EN <- glm(Y_bin ~ PC6 + PC7 + PC8 + PC11 + PC14 + PC18 + PC30 + PC36 + PC38 + PC42 + PC43 + PC52 + PC58 + PC70 + PC92 + PC95 + PC97, data = PCA_100_bin)
summary(log_EN)



# encode_10

encode_10_bin <- cbind(encode_10, Y_bin)
encode_10_bin <- encode_10_bin[-1]

log_full <- glm(Y_bin ~ ., data = encode_10_bin)

log_x <- as.matrix(encode_10_bin[1:10])
log_y <- encode_10_bin[[11]]

model(log_x, log_y)
summary(log_full)



# encode_50

encode_50_bin <- cbind(encode_50, Y_bin)
encode_50_bin <- encode_50_bin[-1]

log_full <- glm(Y_bin ~ ., data = encode_50_bin)

log_x <- as.matrix(encode_50_bin[1:50])
log_y <- encode_50_bin[[51]]

model(log_x, log_y)
summary(log_full)


log_lasso<- glm(Y_bin ~ encoded_feature_0 + encoded_feature_1 + encoded_feature_4 + encoded_feature_7, data = encode_50_bin)
summary(log_lasso)

log_EN <- glm(Y_bin ~ . - encoded_feature_33 - encoded_feature_37 - encoded_feature_38 - encoded_feature_40 - encoded_feature_42 - encoded_feature_44 - encoded_feature_46, data = encode_50_bin)
summary(log_EN)



#encode_100

encode_100_bin <- cbind(encode_100, Y_bin)
encode_100_bin <- encode_100_bin[-1]

log_full <- glm(Y_bin ~ ., data = encode_100_bin)

log_x <- as.matrix(encode_100_bin[1:100])
log_y <- encode_100_bin[[101]]

model(log_x, log_y)
summary(log_full)


log_lasso<- glm(Y_bin ~ encoded_feature_0 + encoded_feature_3 + encoded_feature_10 + encoded_feature_14, data = encode_100_bin)
summary(log_lasso)

log_EN <- glm(Y_bin ~ . - encoded_feature_18 - 
                encoded_feature_29 - 
                encoded_feature_32 - 
                encoded_feature_36 - 
                encoded_feature_38 - 
                encoded_feature_44 - 
                encoded_feature_45 - 
                encoded_feature_46 - 
                encoded_feature_47 - 
                encoded_feature_48 - 
                encoded_feature_50 - 
                encoded_feature_51 -
                encoded_feature_52 -
                encoded_feature_60 -
                encoded_feature_63 -
                encoded_feature_67 -
                encoded_feature_69 -
                encoded_feature_70 -
                encoded_feature_72 -
                encoded_feature_73 -
                encoded_feature_74 -
                encoded_feature_75 -
                encoded_feature_76 -
                encoded_feature_77 -
                encoded_feature_79 -
                encoded_feature_83 -
                encoded_feature_84 -
                encoded_feature_85 -
                encoded_feature_86 -
                encoded_feature_90 -
                encoded_feature_91 -
                encoded_feature_92 -
                encoded_feature_94 -
                encoded_feature_96 -
                encoded_feature_99
              , data = encode_100_bin)
summary(log_EN)

