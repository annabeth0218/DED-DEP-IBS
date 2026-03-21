library(tidyverse)

g.gwas <- read.table("/Users/annabethlu/Projects/gwas_for_prs", header = T)
g.assoc <- read.table("/Users/annabethlu/Projects/matching_QC_assoc_adjust.assoc.adjusted", header = T)

g.gwas_select <- g.gwas[, c("SNP", "A1", "BETA", "P")]
g.assoc_select <- g.assoc[, c("SNP", "CHR")]
arvo <- inner_join(
  x = g.gwas_select,
  y = g.assoc_select,
  by = "SNP"
)

write.table(
  x = arvo,
  file = "26arvo_fuma.txt",
  sep = " ", 
  row.names = FALSE, 
  quote = FALSE 
)

cus.raw <- read_csv("~/Projects/cus_1203.csv")
cus <- cus.raw |> select("Release_No", "GLAUCOMA", "GLAUCOMA_L", "GLAUCOMA_R",
                         "SIT_1_SYSTOLIC_PRESSURE", "SIT_2_SYSTOLIC_PRESSURE", "SIT_3_SYSTOLIC_PRESSURE", "SIT_1_DIASTOLIC_PRESSURE", "SIT_2_DIASTOLIC_PRESSURE", "SIT_3_DIASTOLIC_PRESSURE",
                         "RCCA_plaque", "RCCA_MaxF", "RCCA_AvgF", "RCCA_MinF", "RCCA_AvgMaxF", "RCCA_Ca", "RCCA_r", "RCCA_diameter", "RCCA_abnormal", "RCCA.IMT", "RICA_plaque", "RICA_MaxF", "RICA_AvgF", "RICA_MinF", "RICA_AvgMaxF", "RICA_Ca", "RICA_r", "RICA_diameter", "RICA_abnormal", "LCCA_plaque", "LCCA_MaxF", "LCCA_AvgF", "LCCA_MinF", "LCCA_AvgMaxF", "LCCA_Ca", "LCCA_r", "LCCA_diameter", "LCCA_abnormal", "LCCA.IMT", "LICA_plaque", "LICA_MaxF", "LICA_AvgF", "LICA_MinF", "LICA_AvgMaxF", "LICA_Ca", "LICA_r", "LICA_diameter", "LICA_abnormal")

# cleaning, for R/L CCA/ICA ----
cus <- cus |> mutate(
  LICA_plaque_mm = case_when(
    is.na(LICA_plaque) ~
      NA_real_,
    str_detect(LICA_plaque, "^\\+$") ~
      0.1,
    str_detect(LICA_plaque, regex("cm", ignore_case = TRUE)) ~
      as.numeric(str_replace_all(LICA_plaque, regex("cm|\\s", ignore_case = TRUE), "")) * 10,
    str_detect(LICA_plaque, regex("mm", ignore_case = TRUE)) ~
      as.numeric(str_replace_all(LICA_plaque, regex("mm|\\s", ignore_case = TRUE), "")),
    str_detect(LICA_plaque, fixed("+")) ~ 
      0.1,
    TRUE ~
      as.numeric(LICA_plaque)
  )
)

# cleaning, bp ----
cus <- cus |> rowwise() |> mutate(
  flag = sum(!is.na(SIT_1_SYSTOLIC_PRESSURE), !is.na(SIT_2_SYSTOLIC_PRESSURE), !is.na(SIT_3_SYSTOLIC_PRESSURE)),
  sum_bp = sum(SIT_1_SYSTOLIC_PRESSURE, SIT_2_SYSTOLIC_PRESSURE, SIT_3_SYSTOLIC_PRESSURE, na.rm = TRUE),
  SBP = if_else(
    flag == 0,
    NA_real_,
    sum_bp / flag
  )
) |>
  select(-sum_bp, -flag) |> ungroup()

cus <- cus |> rowwise() |> mutate(
  flag = sum(!is.na(SIT_1_DIASTOLIC_PRESSURE), !is.na(SIT_2_DIASTOLIC_PRESSURE), !is.na(SIT_3_DIASTOLIC_PRESSURE)),
  sum_bp = sum(SIT_1_DIASTOLIC_PRESSURE, SIT_2_DIASTOLIC_PRESSURE, SIT_3_DIASTOLIC_PRESSURE, na.rm = TRUE),
  DBP = if_else(
    flag == 0,
    NA_real_,
    sum_bp / flag
  )
) |>
  select(-sum_bp, -flag) |> ungroup()



# model vars ----
cus.rm.plaque <- c("RCCA_MaxF", "RCCA_AvgF", "RCCA_MinF", "RCCA_AvgMaxF", "RCCA_Ca", "RCCA_r", "RCCA_diameter", "RCCA_abnormal", "RCCA.IMT", 
                      "RICA_MaxF", "RICA_AvgF", "RICA_MinF", "RICA_AvgMaxF", "RICA_Ca", "RICA_r", "RICA_diameter", "RICA_abnormal", 
                      "LCCA_MaxF", "LCCA_AvgF", "LCCA_MinF", "LCCA_AvgMaxF", "LCCA_Ca", "LCCA_r", "LCCA_diameter", "LCCA_abnormal", "LCCA.IMT", 
                      "LICA_MaxF", "LICA_AvgF", "LICA_MinF", "LICA_AvgMaxF", "LICA_Ca", "LICA_r", "LICA_diameter", "LICA_abnormal")
target_var <- "GLAUCOMA"
adjustment_vars <- c("Age", "Sex")

# install.packages(c("pROC", "lmtest", "pscl"))
library(pROC)
library(lmtest) 
library(pscl)

# clean
model_vars <- c(target_var, cus.rm.plaque)
df <- cus |> select(all_of(model_vars)) |> na.omit()
nrow(df[df$GLAUCOMA == 1, ]) # 471
# na.count <- colSums(!is.na(df_complete))
# cus <- cus[, -which(names(cus) == "BP")]

# eval ----
eval <- function(model_obj){
  cat("AIC:", model_obj$aic, "\n")
  cat("Residual Deviance:", model_obj$deviance, "\n")
  cat("BIC:", BIC(model_obj), "\n") # Bayesian Information Criterion
  mcfadden_r2 <- pR2(model_obj)['McFadden']
  cat("r2:", mcfadden_r2, "\n")
  probabilities <- fitted(model_obj)
  matched_outcome <- model_obj$model$GLAUCOMA
  roc_obj <- roc(matched_outcome, probabilities, quiet = T)
  auc_value <- auc(roc_obj)
  cat("AUC: ", auc_value, "\n")
}


# basic model ----
f_full <- as.formula(paste(target_var, " ~ ", paste(cus.rm.plaque, collapse = " + ")))
full <- glm(f_full, data = df, family = binomial(link = "logit"))
eval(full_model)

# AIC: 2388.296 
# Residual Deviance: 2318.296 
# BIC: 2590.032 
# r2: 0.01619088 
# AUC:  0.5889628 

back <- step(full, direction = "backward", trace = FALSE)
eval(back)
# sum.back_htn <- tidy(back_htn)
# sum.back_htn$OR <- exp(sum.back_htn$estimate)
sum <- exp(cbind(OR = coef(back), confint(back)))
sum <- round(sum, 3)
view(sum)

# AIC: 2346.063 
# Residual Deviance: 2330.063 
# BIC: 2392.174 
# r2: 0.01119722 
# AUC:  0.5742643 

# univar ----
library(broom) 
library(dplyr) 

results_list <- list()

for (feature in cus.rm.plaque) {
  f_uni <- as.formula(paste(target_var, " ~ ", feature))
  
  # Fit the logistic regression model using tryCatch for robust error handling
  uni_model <- tryCatch({
    glm(f_uni, data = df_complete, family = binomial(link = "logit"))
  }, error = function(e) {
    cat("Skipping feature due to fundamental glm error:", feature, "\n")
    return(NULL) # Return NULL if glm itself throws an error
  })
  
  # Check if model failed entirely
  if (is.null(uni_model)) {
    next 
  }
  
  # Check for Convergence Failure (often signals separation)
  if (!uni_model$converged) {
    cat("Warning: Model for", feature, "did not converge (possible separation).\n")
    # You can choose to skip or proceed cautiously
    next # Skip this feature to prevent the error
  }
  
  # --- A. Extract Tidy Results ---
  tidy_results <- tidy(uni_model) %>%
    # Filter for the predictor variable (not the intercept)
    filter(term == feature)
  
  # --- B. CRITICAL CHECK: Ensure the predictor coefficient was found ---
  # If the filter returned 0 rows, we skip it to prevent the "differing number of rows" error.
  if (nrow(tidy_results) == 0) {
    cat("Skipping feature:", feature, "- Predictor coefficient not found in model summary.\n")
    next # Skip to the next iteration
  }
  
  # --- C. Proceed with calculations only if the tidy result is valid (1 row) ---
  
  # Extract single values
  log_odds_coef <- tidy_results$estimate
  p_value_val <- tidy_results$p.value
  
  # Calculate Odds Ratio (OR)
  OR <- exp(log_odds_coef)
  
  # Calculate AUC
  probabilities <- fitted(uni_model)
  # The outcome (df_complete$GLAUCOMA) is already matched to probabilities
  roc_obj <- roc(df_complete$GLAUCOMA, probabilities, quiet = TRUE)
  AUC <- auc(roc_obj)
  AUC_numeric <- as.numeric(AUC)
  
  # --- D. Combine and Store ---
  results <- data.frame(
    Feature = feature,
    Log_Odds_Coefficient = log_odds_coef,
    Odds_Ratio = OR,
    AUC = AUC_numeric,
    P_value = p_value_val
  )
  results_list[[feature]] <- results
}

# Combine all results into a single data frame
g.univar <- bind_rows(results_list)





df <- cus |> select(SBP, DBP, GLAUCOMA) |> na.omit()
cat(nrow(df)) # 2780, :)
df$GLAUCOMA <- factor(df$GLAUCOMA, levels = c(0, 1))
levels(df$GLAUCOMA)
model <- glm(GLAUCOMA ~ SBP + DBP, data = df, family = "binomial")
summary(model)
print(exp(coef(model))) # OR
print(exp(confint(model))) # CI
# SBP 1.002 [0.9958-1.0096]
# DBP 0.994 [0.9815-1.0063]

# model <- glm(GLAUCOMA ~ SBP + I(SBP^2) + DBP + I(DBP^2), data = df, family = "binomial")
model_quad <- glm(GLAUCOMA ~ SBP + DBP + I(DBP^2), data = df, family = "binomial")
library(rms)
dd <- datadist(df); options(datadist = 'dd')
model_spline <- lrm(GLAUCOMA ~ rcs(SBP, 4) * rcs(DBP, 4), data = df)
plot_data_sbp <- Predict(model_spline, SBP, ref.zero = TRUE)
plot(plot_data_sbp)
plot_data_dbp <- Predict(model_spline, DBP, ref.zero = TRUE)
plot(plot_data_dbp)

# sbp > 130, dbp >= 80 ----
model_vars <- c(target_var, cus.rm.plaque)
df.htn <- cus |> 
  filter(SBP > 130 & DBP >= 80) |>
  select(all_of(model_vars)) |> na.omit()
nrow(df.htn[df.htn$GLAUCOMA == 1, ]) # 108

f_full <- as.formula(paste(target_var, " ~ ", paste(cus.rm.plaque, collapse = " + ")))
full_htn <- glm(f_full, data = df.htn, family = binomial(link = "logit"))
eval(full_htn)

back_htn <- step(full_htn, direction = "backward", trace = FALSE)
eval(back_htn)
sum.back.htn <- exp(cbind(OR = coef(back_htn), confint(back_htn)))
view(sum.back.htn)

# sbp > 130, dbp < 80----
model_vars <- c(target_var, cus.rm.plaque)
df.weird <- cus |> 
  filter(SBP > 130) |>
  filter(DBP < 80) |>
  select(all_of(model_vars)) |> na.omit()
nrow(df.weird[df.weird$GLAUCOMA == 1, ]) # 91

f_full <- as.formula(paste(target_var, " ~ ", paste(cus.rm.plaque, collapse = " + ")))
full_weird <- glm(f_full, data = df.weird, family = binomial(link = "logit"))
eval(full_weird)

back_weird <- step(full_weird, direction = "backward", trace = FALSE)
eval(back_weird)
sum.back.weird <- exp(cbind(OR = coef(back_weird), confint(back_weird)))
view(sum.back.weird)



# sbp <= 130 ----
model_vars <- c(target_var, cus.rm.plaque)
df.try1 <- cus |> 
  filter(SBP <= 130) |>
  filter(DBP < 70) |>
  select(all_of(model_vars)) |> na.omit()
nrow(df.try1[df.try1$GLAUCOMA == 1, ]) # 208

f_full <- as.formula(paste(target_var, " ~ ", paste(cus.rm.plaque, collapse = " + ")))
full_try1 <- glm(f_full, data = df.try1, family = binomial(link = "logit"))
eval(full_try1)

back_norm <- step(full_norm, direction = "backward", trace = FALSE)
eval(back_norm)
sum.back.norm <- exp(cbind(OR = coef(back_norm), confint(back_norm)))
view(sum.back.norm)

# w/plaque ----
df <- cus
df$RCCA_plaque_mm[is.na(df$RCCA_plaque_mm)] <- 0
df$LCCA_plaque_mm[is.na(df$LCCA_plaque_mm)] <- 0
df$RICA_plaque_mm[is.na(df$RICA_plaque_mm)] <- 0
df$LICA_plaque_mm[is.na(df$LICA_plaque_mm)] <- 0

cus.full <- c(cus.rm.plaque, "RCCA_plaque_mm", "LCCA_plaque_mm", "LICA_plaque_mm", "RICA_plaque_mm")
model_vars <- c(target_var, cus.full)
df.htn.plaque <- df |> 
  filter(SBP > 130 & DBP > 80) |>
  select(all_of(model_vars)) |> na.omit()
nrow(df.htn.plaque[df.htn.plaque$GLAUCOMA == 1, ]) # 104
f_full <- as.formula(paste(target_var, " ~ ", paste(cus.full, collapse = " + ")))
full_htn.plaque <- glm(f_full, data = df.htn.plaque, family = binomial(link = "logit"))
eval(full_htn.plaque)

# AIC: 572.4414 
# Residual Deviance: 496.4414 
# BIC: 736.8351 
# r2: 0.07575389 
# AUC:  0.687743 

# sbp < 120, dbp < 70----
model_vars <- c(target_var, cus.rm.plaque)
df.lbp <- cus |> 
  filter(SBP < 120) |>
  filter(DBP < 70) |>
  select(all_of(model_vars)) |> na.omit()
nrow(df.lbp[df.lbp$GLAUCOMA == 1, ]) # 114

f_full <- as.formula(paste(target_var, " ~ ", paste(cus.rm.plaque, collapse = " + ")))
full_lbp <- glm(f_full, data = df.lbp, family = binomial(link = "logit"))
eval(full_lbp)

back_lbp <- step(full_htn, direction = "backward", trace = FALSE)
eval(back_lbp)
sum.back.lbp <- exp(cbind(OR = coef(back_lbp), confint(back_lbp)))
view(sum.back.lbp)



# multivar
df <- cus |>
  select(all_of(model_vars)) |> na.omit()
vars <- cus.rm.plaque      # your vector of predictor names
target <- target_var

models <- list()

for (i in seq_along(vars)) {
  f <- as.formula(
    paste(target, "~", paste(vars[1:i], collapse = " + "))
  )
  models[[i]] <- glm(f, data = df, family = binomial())
}
sapply(models, AIC)
lr_tests <- lapply(2:length(models), function(i) {
  anova(models[[i-1]], models[[i]], test = "Chisq")
})
lr_tests

cum_coefs <- lapply(models, function(m) {
  tidy(m, exponentiate = TRUE, conf.int = TRUE)
})

library(dplyr)
library(broom)
library(ggplot2)

df_std <- df %>% mutate(across(all_of(cus.rm.plaque), scale))

m_std <- glm(
  as.formula(paste(target_var, "~", paste(cus.rm.plaque, collapse = " + "))),
  data = df_std, family = binomial()
)

coef_df <- tidy(m_std) |> 
  filter(term != "(Intercept)") |> 
  mutate(term = reorder(term, abs(estimate)))

ggplot(coef_df, aes(x = term, y = estimate)) +
  geom_col() +
  coord_flip() +
  labs(x = "", y = "Standardized Coefficient")

or_df <- tidy(full, exponentiate = TRUE, conf.int = TRUE) |>
  filter(term != "(Intercept)") |>
  filter(term != "LICA_r") |> filter(term != "LCCA_r") |> filter(term != "RCCA_r") |>
  mutate(term = reorder(term, estimate))

ggplot(or_df, aes(estimate, term)) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high)) +
  labs(x = "Odds Ratio", y = "")




