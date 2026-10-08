library(dplyr)
library(tidyr)
library(ggplot2)
library(pROC)
library(parallel)
library(doParallel)
library(foreach)
library(kernelshap)
library(treeshap)
library(ggdist)
library(ggbeeswarm)
library(glmnet)

ML_results <- readRDS("results/final_models.rds")
X_train <- ML_results$RF$trainingData
X_train$.outcome <- NULL       
y_train <- factor(ML_results$RF$trainingData$.outcome)

###SHAP analysis
pred_fun <- function(object, newdata) {
  predict(object, newdata = as.matrix(newdata), type = "prob")[, "Disease"]
}
set.seed(42)
n_bg   <- 100
n_test <- 100
bg_X   <- X_train[sample(nrow(X_train), n_bg), , drop = FALSE]
X_shap <- X_train[1:n_test, , drop = FALSE]

pred_fun_lasso <- function(object, newdata) {
  predict(object, newx = as.matrix(newdata), s = "lambda.min",
          type = "response")[, 1]
}
set.seed(42)
shap_lasso <- kernelshap(ML_results$LASSO, X = X_shap,
                         bg_X = bg_X, pred_fun = pred_fun_lasso,
                         parallel = FALSE, verbose = FALSE)
rf_unified <- randomForest.unify(ML_results$RandomForest$finalModel, X_train)
rf_shap_raw <- treeshap(rf_unified, X_shap, verbose = FALSE)
shap_rf <- list(S = as.matrix(rf_shap_raw$shaps[, colnames(X_shap)]),
                X = X_shap, baseline = attr(rf_shap_raw$shaps, "baseline"))

xgb_shap_matrix <- predict(ML_results$XGBoost,
                           newdata = xgb.DMatrix(as.matrix(X_shap)),
                           predcontrib = TRUE)
colnames(xgb_shap_matrix)[1:ncol(X_shap)] <- colnames(X_shap)
shap_xgb <- list(S = xgb_shap_matrix[, 1:ncol(X_shap)],
                 X = X_shap, baseline = xgb_shap_matrix[1, ncol(X_shap) + 1])

shap_svm <- kernelshap(ML_results$SVM, X = X_shap,
                       bg_X = bg_X, pred_fun = pred_fun, parallel = FALSE, verbose = FALSE)

extract_imp <- function(obj, name) {
  obj$S |> as.data.frame() |>
    summarise(across(everything(), ~ mean(abs(.x)))) |>
    pivot_longer(everything(), names_to = "Gene", values_to = "SHAP_Importance") |>
    mutate(Model = name, .before = 1)
}

shap_importance <- bind_rows(
  extract_imp(shap_lasso, "LASSO"),
  extract_imp(shap_rf,    "RandomForest"),
  extract_imp(shap_xgb,   "XGBoost"),
  extract_imp(shap_svm,   "SVM")
)

shap_long <- bind_rows(
  shap_lasso$S |> as.data.frame() |> mutate(.id = row_number(), Model = "LASSO"),
  shap_rf$S    |> as.data.frame() |> mutate(.id = row_number(), Model = "RF"),
  shap_xgb$S   |> as.data.frame() |> mutate(.id = row_number(), Model = "XGBoost"),
  shap_svm$S   |> as.data.frame() |> mutate(.id = row_number(), Model = "SVM")
) |>
  pivot_longer(-c(.id, Model), names_to = "Gene", values_to = "SHAP")

top_genes_overall <- shap_importance |>
  group_by(Gene) |>
  summarise(MeanImp = mean(SHAP_Importance), .groups = "drop") |>
  slice_max(MeanImp, n = 15) |>
  pull(Gene)

shap_top <- shap_long |> filter(Gene %in% top_genes_overall)

p_shap_summary <- ggplot(shap_top, aes(x = SHAP, y = reorder(Gene, abs(SHAP), FUN = median))) +
  ggdist::stat_halfeye(aes(fill = Model), alpha = 0.6,
                       position = position_dodge(0.7), width = 0.6) +
  labs(x = "SHAP Value", y = NULL, title = "SHAP value distribution (top 15 genes)") +
  scale_fill_manual(values = c("LASSO" = "red", "RF" = "#000000",
                                "XGBoost" = "#069DFF", "SVM" = "#A4E048")) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), axis.line = element_line(color = "black"),
        legend.position = "bottom")
p_shap_summary

top20_shap <- shap_importance |>
  group_by(Gene) |>
  summarise(MeanSHAP = mean(SHAP_Importance), .groups = "drop") |>
  slice_max(MeanSHAP, n = 20)

p_shap_imp <- ggplot(top20_shap, aes(x = MeanSHAP, y = reorder(Gene, MeanSHAP))) +
  geom_col(width = 0.7, fill = "#069DFF", alpha = 0.85) +
  labs(x = "Mean |SHAP|", y = NULL, title = "Top 20 genes by SHAP importance") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), axis.line = element_line(color = "black"))

X_shap_long <- X_shap |>
  as.data.frame() |>
  mutate(.id = row_number()) |>
  pivot_longer(-.id, names_to = "Gene", values_to = "FeatureValue") |>
  group_by(Gene) |>
  mutate(FeatureNorm = (FeatureValue - min(FeatureValue)) /
           (max(FeatureValue) - min(FeatureValue))) |>
  ungroup() |>
  select(-FeatureValue)

shap_beeswarm <- shap_long |>
  left_join(X_shap_long, by = c(".id", "Gene"))

beeswarm_top_genes <- shap_importance |>
  group_by(Gene) |>
  summarise(MeanImp = mean(SHAP_Importance), .groups = "drop") |>
  slice_max(MeanImp, n = 15) |>
  pull(Gene)
shap_beeswarm <- shap_beeswarm |>
  select(-any_of("MeanImp")) |>
  left_join(gene_imp, by = "Gene")

## Figure 3E
p_shap_beeswarm <- shap_beeswarm |>
  filter(Gene %in% beeswarm_top_genes) |>
  ggplot(aes(x = SHAP, y = reorder(Gene, MeanImp), color = FeatureNorm)) +
  geom_quasirandom(alpha = 0.75, size = 1.6, width = 0.2) +
  scale_color_gradient2(low = "#0068FF", mid = "#FFFFFF", high = "#FF2E00",
                        midpoint = 0.5, name = "Feature value\n(normalized)") +
  facet_wrap(~ Model, nrow = 2) +
  labs(x = "SHAP value", y = NULL,
       title = "SHAP beeswarm plot") +
  theme_minimal(base_size = 12) +
  theme(
  panel.grid = element_blank(),
  axis.line = element_line(color = "black"),
  axis.text = element_text(color = "black", size = 14, face = "plain"),
  axis.title = element_text(color = "black", size = 14, face = "plain"),
  strip.background = element_blank(),
  legend.position = "right"
)
p_shap_beeswarm
