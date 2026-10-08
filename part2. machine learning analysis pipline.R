library(Seurat)
library(dplyr)
library(ggplot2)
library(caret)
library(glmnet)
library(randomForest)
library(xgboost)
library(e1071)
library(Boruta)
library(pROC)
library(ROCR)
library(doParallel)
library(ggrepel)
library(tidyr)

results_dir <- "results"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

### 1. reading data 
set.seed(42)
Macro <- readRDS("macro.rds")
DefaultAssay(Macro) <- "SCT"
Macro$donor <- as.factor(Macro$orig.ident)
print(table(Macro$donor, Macro$disease))    
expr_matrix <- GetAssayData(Macro, layer = "data")
expr_matrix <- as.matrix(expr_matrix)
y_binary <- ifelse(Macro$disease == "normal", 0, 1)
table(y_binary)

### 2. obtaining DEGs between disease and normal 
Idents(Macro) <- "disease"
deg <- FindMarkers(Macro,
                   ident.1 = "atrial fibrillation || mitral valve insufficiency",
                   ident.2 = "normal",
                   min.pct = 0.1,
                   logfc.threshold = 0.25,    
                   only.pos = FALSE)
deg$gene <- rownames(deg)
deg_sig <- deg[deg$p_val_adj < 0.05, ]
cat("total DEGs:", nrow(deg_sig), "\n")
cat("up:", sum(deg_sig$avg_log2FC > 0), "\n")
cat("down:", sum(deg_sig$avg_log2FC < 0), "\n")
candidate_genes <- rownames(deg_sig)
X <- t(expr_matrix[candidate_genes, ])
y <- factor(y_binary, levels = c(0, 1), labels = c("Normal", "Disease"))

### 3. train-test split based on donors
set.seed(42)
donors <- levels(Macro$donor)
train_donors <- sample(donors, size = round(0.7 * length(donors)))
test_donors  <- setdiff(donors, train_donors)

train_cells <- which(Macro$donor %in% train_donors)
test_cells  <- which(Macro$donor %in% test_donors)

X_train <- X[train_cells, , drop = FALSE]
y_train <- y[train_cells]
X_test  <- X[test_cells, , drop = FALSE]
y_test  <- y[test_cells]

train_donor_vec <- Macro$donor[train_cells]
group_folds <- groupKFold(train_donor_vec, k = 5)
donor_tbl <- data.frame(donor = train_donor_vec, class = y_train) |>
  dplyr::distinct()
table(donor_tbl$class)
n_minor_donors <- min(table(donor_tbl$class))
if (n_minor_donors < 2) {
  stop("Donor-level CV unavailable: not enough donors in the minority class. Consider using a different splitting strategy.")
}
k_folds <- min(5, n_minor_donors)
set.seed(42)
donor_tbl$fold <- NA_integer_
for (cl in levels(y_train)) {
  idx <- which(donor_tbl$class == cl)
  donor_tbl$fold[idx] <- sample(rep(seq_len(k_folds), length.out = length(idx)))
}
foldid <- donor_tbl$fold[match(train_donor_vec, donor_tbl$donor)]
n_cores <- max(1, detectCores() - 1)
registerDoParallel(cores = n_cores)

### 4. multi-feature selection
selected_features <- list()

## 4.1 LASSO ----------------------------------------------------------------
set.seed(42)
cv_lasso <- cv.glmnet(
  x = X_train,
  y = y_train,
  family = "binomial",
  alpha = 1,
  foldid = foldid,          
  parallel = TRUE
)
lasso_coef <- coef(cv_lasso, s = "lambda.1se")
lasso_genes <- rownames(lasso_coef)[which(lasso_coef[, 1] != 0)]
lasso_genes <- setdiff(lasso_genes, "(Intercept)")
lasso_coef_full <- as.matrix(coef(cv_lasso, s = "lambda.min"))
lasso_coef_full <- lasso_coef_full[rownames(lasso_coef_full) != "(Intercept)", , drop = FALSE]
lasso_imp <- data.frame(
  gene = rownames(lasso_coef_full),
  coef = as.numeric(lasso_coef_full),
  abs_coef = abs(as.numeric(lasso_coef_full))
)
lasso_imp <- data.frame(gene = rownames(lasso_coef_full),
                          coef = as.numeric(lasso_coef_full),
                          abs_coef = abs(as.numeric(lasso_coef_full))) |>
    arrange(desc(abs_coef))
## Figure 2B
top20_lasso <- lasso_imp |>
  filter(coef != 0) |>
  slice_head(n = 20)
selected_features[["LASSO"]] <- top20_lasso$gene

ggplot(top20_lasso, aes(x = reorder(gene, coef), y = coef, fill = coef > 0)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"),
                    labels = c("down", "up")) +
  labs(x = NULL, y = "Coefficient (lambda.min)", fill = "black",
       title = "LASSO Top Genes") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.text = element_text(size = 13)
  )

## Figure 2C
top20_genes_lasso <- intersect(top20_lasso$gene, rownames(cv_lasso$glmnet.fit$beta))
beta_all <- as.matrix(cv_lasso$glmnet.fit$beta)          # genes x lambda
beta_top20 <- beta_all[top20_genes_lasso, , drop = FALSE] # 20 x lambda
matplot(
  x = log(cv_lasso$lambda),
  y = t(beta_top20),
  type = "l",
  lty = 1,
  col = rainbow(length(top20_genes_lasso)),
  lwd = 1.5,
  xaxs = "i",
  xlab = expression(log(lambda)),
  ylab = "Coefficient",
  main = "LASSO Coefficient Path",
  cex.main = 1.6, font.main = 1,
  cex.lab  = 1.4,  
  cex.axis = 1.2, 
)

## 4.2 Random Forest --------------------------------------------------------
ntree_val <- ifelse(ncol(X_train) > 500, 300, 500)
set.seed(42)
rf_imp <- randomForest(
  x = X_train,
  y = y_train,
  ntree = ntree_val,
  importance = TRUE,
  mtry = max(1, floor(sqrt(ncol(X_train))))
)
rf_importance <- as.data.frame(rf_imp$importance)
rf_importance$gene <- rownames(rf_importance)
rf_importance <- rf_importance[order(rf_importance$MeanDecreaseGini, decreasing = TRUE), ]
rf_cumsum <- cumsum(rf_importance$MeanDecreaseGini) / sum(rf_importance$MeanDecreaseGini)
n_rf_select <- min(max(which(rf_cumsum <= 0.90)), 200)
n_rf_select <- max(n_rf_select, 20)
selected_features[["RF"]] <- rf_importance$gene[1:n_rf_select]

## Figure 2D
top20_rf <- rf_importance |>
  arrange(desc(MeanDecreaseGini)) |>
  slice_head(n = 20) |>
  mutate(gene = factor(gene, levels = gene))

p1 <- ggplot(top20_rf, aes(x = MeanDecreaseGini, y = gene, fill = MeanDecreaseGini)) +
  geom_col(width = 0.7) +
  scale_fill_gradient(low = "steelblue", high = "#E41A1C") +
  labs(x = "MeanDecreaseGini", y = NULL, fill = NULL,
       title = "RF Top Genes") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.y = element_text(size = 13),        
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    plot.title = element_text(size = 18, hjust = 0.5)
  )
p1
## Figure 2E
rf_cum_df <- data.frame(n_genes = 1:length(rf_cumsum), cum_imp = rf_cumsum)
n90 <- which(rf_cumsum >= 0.9)[1]

p2 <- ggplot(rf_cum_df, aes(x = n_genes, y = cum_imp)) +
  geom_line(color = "#377EB8", linewidth = 1) +
  geom_vline(xintercept = n90, linetype = "dashed", color = "grey50") +
  annotate("text", x = n90 + 80, y = 0.5,
           label = paste0(n90, " genes\n(90% importance)"),
           size = 4.5, color = "grey30") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Number of Genes", y = "Cumulative Importance",
       title = "RF Cumulative Importance") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    plot.title = element_text(size = 16, face = "plain", hjust = 0.5),
    axis.title = element_text(size = 14, face = "plain"),
    axis.text  = element_text(size = 12, face = "plain", color = "black")
  )
p2

## 4.3 XGBoost
dtrain <- xgb.DMatrix(data = X_train, label = ifelse(y_train == "Disease", 1, 0))
set.seed(42)
xgb_model <- xgb.train(
  params = list(
    objective      = "binary:logistic",
    eval_metric    = "auc",
    max_depth      = 4,
    eta            = 0.1,
    subsample      = 0.8,
    colsample_bytree = 0.8
  ),
  data   = dtrain,
  nrounds = 100,
  verbose = 0
)
xgb_importance <- xgb.importance(model = xgb_model)
xgb_cumsum <- cumsum(xgb_importance$Gain) / sum(xgb_importance$Gain)
n_xgb_select <- min(max(which(xgb_cumsum <= 0.90)), 200)
n_xgb_select <- max(n_xgb_select, 20)
selected_features[["XGBoost"]] <- xgb_importance$Feature[1:n_xgb_select]

## Figure 2F
top20_xgb <- xgb_importance |>
  arrange(desc(Gain)) |>
  slice_head(n = 20) |>
  mutate(Feature = factor(Feature, levels = Feature))

p3 <- ggplot(top20_xgb, aes(x = Gain, y = Feature, fill = Gain)) +
  geom_col(width = 0.7) +
  scale_fill_gradient(low = "steelblue", high = "#E41A1C") +
  labs(x = "Gain", y = NULL, fill = NULL, title = "XGBoost Top Genes") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.y = element_text(size = 13),  
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    plot.title = element_text(size = 18, hjust = 0.5)
  )
p3
## Figure 2G
xgb_cum_df <- data.frame(n_genes = 1:length(xgb_cumsum), cum_imp = xgb_cumsum)
n90_xgb <- which(xgb_cumsum >= 0.9)[1]
if (is.na(n90_xgb)) n90_xgb <- length(xgb_cumsum)

p4 <- ggplot(xgb_cum_df, aes(x = n_genes, y = cum_imp)) +
  geom_line(color = "#4DAF4A", linewidth = 1) +
  geom_vline(xintercept = n90_xgb, linetype = "dashed", color = "grey50") +
  annotate("text", x = n90_xgb + 80, y = 0.5,
           label = paste0(n90_xgb, " genes\n(90% importance)"),
           size = 4.5, color = "black") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Number of Genes", y = "Cumulative Importance",
       title = "XGBoost Cumulative Importance") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    plot.title = element_text(size = 16, face = "plain", hjust = 0.5),
    axis.title = element_text(size = 14, face = "plain"),
    axis.text  = element_text(size = 12, face = "plain", color = "black")
  )
p4

## 4.4 SVM
set.seed(42)
svm_model <- svm(
  x = X_train,
  y = y_train,
  kernel = "linear",
  cost = 1,
  scale = TRUE
)
svm_weights <- t(svm_model$coefs) %*% svm_model$SV
svm_imp <- data.frame(
  gene = colnames(X_train),
  weight = as.numeric(svm_weights)
) |>
  mutate(abs_weight = abs(weight)) |>
  arrange(desc(abs_weight))

svm_cumsum <- cumsum(svm_imp$abs_weight) / sum(svm_imp$abs_weight)
n_svm_select <- min(max(which(svm_cumsum <= 0.90)), 200)
n_svm_select <- max(n_svm_select, 20)
selected_features[["SVM"]] <- svm_imp$gene[1:n_svm_select]

## Figure 2H
top20_svm <- svm_imp |>
  slice_head(n = 20) |>
  mutate(gene = factor(gene, levels = gene))
p5 <- ggplot(top20_svm, aes(x = abs_weight, y = gene, fill = abs_weight)) +
  geom_col(width = 0.7) +
  scale_fill_gradient(low = "steelblue", high = "#E41A1C") +
  labs(x = "|Weight|", y = NULL, fill = NULL, title = "SVM Top Genes") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.y = element_text(size = 13),     
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    plot.title = element_text(size = 18, hjust = 0.5)
  )
p5

svm_cum_df <- data.frame(n_genes = 1:length(svm_cumsum), cum_imp = svm_cumsum)
n90_svm <- which(svm_cumsum >= 0.9)[1]
if (is.na(n90_svm)) n90_svm <- length(svm_cumsum)

p6 <- ggplot(svm_cum_df, aes(x = n_genes, y = cum_imp)) +
  geom_line(color = "#984EA3", linewidth = 1) +
  geom_vline(xintercept = n90_svm, linetype = "dashed", color = "black") +
  annotate("text", x = n90_svm + 80, y = 0.5,
           label = paste0(n90_svm, " genes\n(90% importance)"),
           size = 4.5, color = "black") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Number of Genes", y = "Cumulative Importance",
       title = "SVM Cumulative Importance") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    plot.title = element_text(size = 16, face = "plain", hjust = 0.5),
    axis.title = element_text(size = 14, face = "plain"),
    axis.text  = element_text(size = 12, face = "plain", color = "black")
  )
p6

## 4.5 Boruta ---------------------------------------------------------------
boruta_genes <- candidate_genes
if (length(boruta_genes) > 500) {
  rf_genes_ranked <- rownames(rf_importance)
  boruta_genes <- intersect(rf_genes_ranked[1:500], candidate_genes)
}
boruta_X <- X_train[, boruta_genes, drop = FALSE]
set.seed(42)
boruta_result <- Boruta(
  x = boruta_X,
  y = y_train,
  maxRuns = 100,
  doTrace = 0,
  ntree = 100
)
boruta_confirmed <- names(boruta_result$finalDecision[boruta_result$finalDecision == "Confirmed"])
selected_features[["Boruta"]] <- boruta_confirmed

## Figure 2I
imp_history <- as.data.frame(boruta_result$ImpHistory)
confirmed_names <- intersect(boruta_confirmed, colnames(imp_history))
med_imp <- sort(sapply(imp_history[confirmed_names], median, na.rm = TRUE), decreasing = TRUE)
top_b <- names(head(med_imp, 20))

top_features <- imp_history[, top_b, drop = FALSE] |>
  tidyr::pivot_longer(everything(), names_to = "Feature", values_to = "Importance")
top_features$Feature <- factor(top_features$Feature, levels = top_b)

p1_boruta <- ggplot(top_features, aes(x = Feature, y = Importance, fill = "Confirmed")) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3) +
  scale_fill_manual(values = c("Confirmed" = "#1E90FF"), name = "Decision") +
  geom_hline(yintercept = max(imp_history[, "shadowMax"], na.rm = TRUE),
             linetype = "dashed", color = "#666666", linewidth = 0.6) +
  labs(x = NULL, y = "Importance (Z-score)",
       title = "Boruta Feature Importance (Confirmed)") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                               size = 10, face = "plain")
  )
p1_boruta


for (nm in names(selected_features)) {
  cat("  ", nm, ":", length(selected_features[[nm]]), "\n")
}
all_genes <- unique(unlist(selected_features))
gene_votes <- data.frame(gene = all_genes, count = 0)
for (nm in names(selected_features)) {
  gene_votes$count <- gene_votes$count + (all_genes %in% selected_features[[nm]])
}
consensus_2 <- gene_votes$gene[gene_votes$count >= 2]
consensus_3 <- gene_votes$gene[gene_votes$count >= 3]
cat("selected by at least two methods:", length(consensus_2), "\n")
cat("selected by at least three methods:", length(consensus_3), "\n")
gene_score <- data.frame(gene = all_genes, score = 0)
for (nm in names(selected_features)) {
  n_genes <- length(selected_features[[nm]])
  ranks <- seq_len(n_genes)
  weights <- (n_genes - ranks + 1) / n_genes
  matched <- match(all_genes, selected_features[[nm]])
  gene_score$score <- gene_score$score + ifelse(is.na(matched), 0, weights[matched])
}
gene_score <- gene_score[order(gene_score$score, decreasing = TRUE), ]
gene_score$vote_count <- gene_votes$count[match(gene_score$gene, gene_votes$gene)]
print(head(gene_score, 20))

### 5. Final model construction with consensus features --------------------

## 5.1 Final feature set
final_features <- consensus_2              # selected by >= 2 methods
# 更保守可换用: final_features <- consensus_3
X_train_final <- X_train[, final_features, drop = FALSE]
X_test_final  <- X_test[, final_features, drop = FALSE]
y_train_num <- ifelse(y_train == "Disease", 1, 0)
y_test_num  <- ifelse(y_test  == "Disease", 1, 0)

## 5.2 Donor-grouped CV (reuse the donor-level folds from Section 3)
foldid_final <- foldid

ctrl <- trainControl(
  method = "cv",
  number = k_folds,
  index = split(seq_len(nrow(X_train_final)), foldid_final),
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

## 5.3 LASSO logistic regression ------------------------------------------
set.seed(42)
cv_lasso_final <- cv.glmnet(
  x = X_train_final,
  y = y_train_num,
  family = "binomial",
  alpha = 1,
  foldid = foldid_final,
  parallel = TRUE
)
lasso_final_pred <- predict(cv_lasso_final, newx = X_test_final,
                            s = "lambda.min", type = "response")[, 1]
lasso_final_coef <- as.matrix(coef(cv_lasso_final, s = "lambda.min"))
lasso_final_genes <- rownames(lasso_final_coef)[
  rownames(lasso_final_coef) != "(Intercept)" & lasso_final_coef[, 1] != 0]

## 5.4 Random forest --------------------------------------------------------
set.seed(42)
model_rf <- train(
  x = X_train_final,
  y = y_train,
  method = "rf",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5,
  importance = TRUE
)
rf_final_pred <- predict(model_rf, newdata = X_test_final,
                         type = "prob")[, "Disease"]

## 5.5 Linear SVM -----------------------------------------------------------
set.seed(42)
model_svm <- train(
  x = X_train_final,
  y = y_train,
  method = "svmLinear",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5,
  preProc = c("center", "scale")
)
svm_final_pred <- predict(model_svm, newdata = X_test_final,
                          type = "prob")[, "Disease"]

## 5.6 XGBoost --------------------------------------------------------------
dtrain_final <- xgb.DMatrix(data = X_train_final, label = y_train_num)
dtest_final  <- xgb.DMatrix(data = X_test_final,  label = y_test_num)

# class weights to handle the imbalance
n_pos <- sum(y_train_num == 1)
n_neg <- sum(y_train_num == 0)
class_weights <- ifelse(y_train_num == 1, n_neg / n_pos, 1)
setinfo(dtrain_final, "weight", class_weights)

xgb_params <- list(
  objective       = "binary:logistic",
  eval_metric     = "auc",
  max_depth       = 4,
  eta             = 0.1,
  subsample       = 0.8,
  colsample_bytree = 0.8
)
set.seed(42)
cv_xgb <- xgb.cv(
  params = xgb_params,
  data = dtrain_final,
  nrounds = 200,
  nfold = k_folds,
  folds = split(seq_len(nrow(X_train_final)), foldid_final),
  early_stopping_rounds = 20,
  verbose = 0
)
best_nrounds <- cv_xgb$early_stop$best_iteration

set.seed(42)
xgb_final <- xgb.train(
  params = xgb_params,
  data = dtrain_final,
  nrounds = best_nrounds,
  verbose = 0
)
xgb_final_pred <- predict(xgb_final, dtest_final)

## 5.7 Collect models and predictions ---------------------------------------
models_results <- list(
  LASSO   = cv_lasso_final,
  RF      = model_rf,
  SVM     = model_svm,
  XGBoost = xgb_final
)
pred_list <- list(
  LASSO   = lasso_final_pred,
  RF      = rf_final_pred,
  SVM     = svm_final_pred,
  XGBoost = xgb_final_pred
)

### 6. Evaluation on the held-out donor test set ----------------------------

compute_pr_auc <- function(prob, truth) {
  pred_obj <- prediction(prob, truth)
  perf <- performance(pred_obj, "prec", "rec")
  x <- perf@x.values[[1]]   # recall
  y <- perf@y.values[[1]]   # precision
  valid <- is.finite(x) & is.finite(y)
  x <- x[valid]; y <- y[valid]
  o <- order(x)
  x <- x[o]; y <- y[o]
  sum(diff(x) * (y[-1] + y[-length(y)]) / 2)   # trapezoidal integration
}

roc_list <- lapply(pred_list, function(p) roc(y_test_num, p, quiet = TRUE))

eval_results <- do.call(rbind, lapply(names(pred_list), function(nm) {
  p <- pred_list[[nm]]
  r <- roc_list[[nm]]
  thr <- coords(r, "best", ret = "threshold", best.method = "youden")
  thr <- as.numeric(thr$threshold[1])
  pred_class <- factor(ifelse(p >= thr, "Disease", "Normal"),
                       levels = levels(y_test))
  cm <- confusionMatrix(pred_class, y_test, positive = "Disease")
  data.frame(
    Model             = nm,
    AUC               = as.numeric(r$auc),
    PR_AUC            = compute_pr_auc(p, y_test_num),
    Threshold         = thr,
    Accuracy          = cm$overall["Accuracy"],
    Sensitivity       = cm$byClass["Sensitivity"],
    Specificity       = cm$byClass["Specificity"],
    Balanced_Accuracy = cm$byClass["Balanced Accuracy"],
    Precision         = cm$byClass["Precision"],
    F1                = cm$byClass["F1"]
  )
}))
rownames(eval_results) <- NULL
print(eval_results)

### 7. Final-model feature importance ---------------------------------------

lasso_imp_final <- data.frame(gene = rownames(lasso_final_coef),
                              coef = lasso_final_coef[, 1]) |>
  filter(gene != "(Intercept)", coef != 0) |>
  mutate(abs_coef = abs(coef))

rf_imp_vals <- data.frame(gene = rownames(varImp(model_rf)$importance),
                          importance = varImp(model_rf)$importance[, 1]) |>
  arrange(desc(importance))

svm_coef <- model_svm$finalModel@coef
svm_sv   <- model_svm$finalModel@xmatrix
if (is.list(svm_coef)) svm_coef <- svm_coef[[1]]
if (is.list(svm_sv))   svm_sv   <- svm_sv[[1]]

svm_w_final <- t(svm_coef) %*% svm_sv

svm_imp_vals <- data.frame(gene = colnames(svm_w_final),
                           weight = as.numeric(svm_w_final)) |>
  mutate(abs_weight = abs(weight)) |>
  arrange(desc(abs_weight))

xgb_imp_vals <- xgb.importance(model = xgb_final)

top_genes_list <- list(
  LASSO   = lasso_imp_final |> arrange(desc(abs_coef)) |> slice_head(n = 20),
  RF      = rf_imp_vals   |> slice_head(n = 20),
  SVM     = svm_imp_vals  |> slice_head(n = 20),
  XGBoost = xgb_imp_vals  |> arrange(desc(Gain)) |> slice_head(n = 20)
)

imp_plot_df <- bind_rows(
  data.frame(Model = "LASSO",   Gene = top_genes_list$LASSO$gene,      Imp = top_genes_list$LASSO$abs_coef),
  data.frame(Model = "RF",      Gene = top_genes_list$RF$gene,         Imp = top_genes_list$RF$importance),
  data.frame(Model = "SVM",     Gene = top_genes_list$SVM$gene,        Imp = top_genes_list$SVM$abs_weight),
  data.frame(Model = "XGBoost", Gene = top_genes_list$XGBoost$Feature, Imp = top_genes_list$XGBoost$Gain)
)
# per-panel ordering without extra packages
imp_plot_df$key <- reorder(paste(imp_plot_df$Model, imp_plot_df$Gene, sep = "___"),
                           imp_plot_df$Imp)
p_bar <- ggplot(imp_plot_df, aes(x = key, y = Imp, fill = Model)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  facet_wrap(~ Model, scales = "free") +
  scale_x_discrete(labels = function(x) sub("^.*___", "", x)) +
  coord_flip() +
  labs(title = "Top 20 Genes by Model Importance", x = NULL, y = "Importance") +
  theme_minimal(base_size = 12)
p_bar

### 9. Save results ----------------------------------------------------------
saveRDS(models_results, file.path(results_dir, "final_models.rds"))
write.csv(eval_results, file.path(results_dir, "test_metrics.csv"), row.names = FALSE)
write.csv(gene_score,   file.path(results_dir, "consensus_gene_scores.csv"), row.names = FALSE)

