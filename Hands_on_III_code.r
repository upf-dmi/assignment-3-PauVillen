# Assignment 3
setwd("/Users/pauvillen14/Desktop/BIOINFO/DMI/assignment-3-PauVillen")

# ---------------------------------------------------------------
# --- EXERCISE 1 ------------------------------------------------
# ---------------------------------------------------------------
# Load all the required libraries
library(tidyverse)
library(readxl)
library(janitor) 
library(dplyr)
library(caret)
library(recipes)
library(pROC)
library(rsample)
library(randomForest)
library(mice)

# --- LOAD DATA ---
df3 <- read_excel("table_s3.xlsx", sheet = 2)
df3_clean <- clean_names(df3)

df1 <- read_excel("table_s1.xlsx", sheet = 2)
df1_clean <- clean_names(df1)

# --- CLEAN DATA ---
new_headers <- df3_clean[1, ] %>%
  select(patient_id:last_col()) %>% 
  unlist() %>%
  as.character()

colnames(df3_clean)[3:ncol(df3_clean)] <- new_headers

df3_clean_final <- df3_clean %>%
  slice(-1) %>%                         # Remove that first row now that we've used it for names
  pivot_longer(
    cols = all_of(new_headers),         # Pivot all the patient columns (XG1, XG2...)
    names_to = "patient_id",            # New column to hold the patient names
    values_to = "expression"            # New column to hold the numeric values
  ) %>%
  select(-gene_symbol) %>%
  pivot_wider(
    names_from = proteins_metabolites,          
    values_from = expression
  ) %>%
  # Convert the protein columns to numeric for the random forest
  mutate(across(-patient_id, as.numeric))

colnames(df1_clean)[1] <- "patient_id" #changed the name from patient_id_a to patient_id 
df1_clean <- df1_clean %>% select(patient_id, group_d)

# Merge clinical metadata
df_final <- left_join(df3_clean_final, df1_clean, by = "patient_id")
df_final <- df_final %>%
  filter(group_d %in% c(2, 3)) %>%
  mutate(group_d = factor(group_d, levels = c(2,3), labels = c("non_Severe", "Severe"))) %>%
  select(-patient_id) #remove ID so it's not used as a predictor
table(df_final$group_d) #to see how many non-severe and severe patients are there and if we have removed the other groups correctly
#Now df_final contains all proteins in columns, patients in rows and in the last column it appears group_d (and only contians non-severe and severe patients)

# --- MISSING DATA HANDLING ---
# Exploratory NA analysis
sum(is.na(df_final)) #Total NA values

# Remove proteins with > 20% missing values
na_threshold <- 0.2
cols_to_keep <- colMeans(is.na(df_final)) <= na_threshold
df_final_0.2_removed <- df_final[ , cols_to_keep]


# Median imputation
df_final_clean <- df_final_0.2_removed %>%
  mutate(across(where(is.numeric), ~ifelse(is.na(.), median(., na.rm = TRUE), .)))
df_final_clean <- df_final_clean %>%
  clean_names()

#  --- INITIAL RF MODEL ---
rf_model <- randomForest(as.factor(group_d) ~ .,
                         data = df_final_clean,
                         importance = TRUE, #essential for exercise 2, where need to identify the top 25 proteins. so it tells R to calculate which proteins were most helpful in splitting the severe from non-severe groups
                         ntree = 500)

print(rf_model)

# --- STRUCTURED ML WORKFLOW ---
df_final_preclean <- df_final %>% janitor::clean_names()

#--- COHORT SPLIT & INTERNAL DIVISION ------
set.seed(123)

#Create an 80/20 split stratified by the outcome (group_d)
data_split <- initial_split(df_final_preclean, prop = 0.8, strata = group_d)
#Development cohort (80%)
training_data <- training(data_split) #will have 24 observation, which is a 0.8% of the 31 ones initially
#Internal testing (20%)
internal_test_data <- testing(data_split)

# Cross-Validation setup (10-Fold CV)

#Define the recipe (blueprint) (THIS IS DONE IN THE TEACHER SLIDES!!!!)
rf_recipe <- recipe(group_d ~ ., data = training_data) %>%
  step_nzv(all_predictors()) %>%
  step_impute_median(all_numeric_predictors()) %>%
  step_normalize(all_numeric_predictors())

# # set up training controls
# This handles the 10-fold CV and collects data for the ROC curve
ctrl <- trainControl(
  method = "cv",
  number = 10,
  savePredictions = "final",       # Necessary to plot the ROC later
  classProbs = TRUE,               # Necessary for AUC calculation
)

# Hyperparameter tuning
# Define tuning grid
# Testing different numbers of predictors at each split
p <- ncol(df_final_preclean) - 1 
mtry_grid <- data.frame(mtry = c(1, 2, floor(sqrt(p)), floor(p/3), 10, 30))

# Model training with caret
rf_model_2 <- train(
  rf_recipe, 
  data = training_data, 
  method = "rf",
  trControl = ctrl,
  tuneGrid = mtry_grid,
  metric = "Accuracy",                 
  ntree = 500,
  importance = TRUE
)

# --- MODEL EVALUATION ---
# Performance during Cross-validation (the 80%)
print(rf_model_2)

# Performance on the Internal Hold-out (the 20%)
## This proves the model works on "unseen" data before going to the external cohort
internal_preds <- predict(rf_model_2, newdata = internal_test_data)
internal_cm <- confusionMatrix(internal_preds, internal_test_data$group_d)

print("----Internal 20% test results:-----")
print(internal_cm)

# Plot the Realistic ROC Curve
# Using the cross-validated predictions (prevents the "perfect square" look)
roc_obj <- roc(rf_model_2$pred$obs, rf_model_2$pred$Severe)

# Draw the plot (this creates 'plot.new' inside the file)
{
  plot(roc_obj, col = "blue", lwd = 3, main = "ROC Curve (10-Fold Cross-Validation)")
  
  #Add the lines (these will now work)
  abline(a = 0, b = 1, lty = 2, col = "grey")
  
  # Add the legend
  legend("bottomright", 
         legend = paste("AUC =", round(auc(roc_obj), 3)), 
         col = "blue", 
         lwd = 3, 
         bty = "n")
  
}


# ---------------------------------------------------------------
# --- EXERCISE 2 ------------------------------------------------
# ---------------------------------------------------------------
library(pheatmap)
library(tidyverse)
library(caret)
library(janitor)
library(dplyr) #for calibration plot
library(ggplot2) #for calibration plot
library(iml) #for SHAP analysis
library(fastshap) #for SHAP analysis

# --- TOP 25 PROTEINS AND IMPORTING DATA ---
# extract importance
importance_results <- varImp(rf_model_2, scale = FALSE)

# create top 25 list
top_25_list <- importance_results$importance %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_ID") %>%
  # Fix the case: Ensure Protein_ID is Uppercase to match the metadata
  mutate(Protein_ID = toupper(Protein_ID)) %>%
  arrange(desc(Severe)) %>% 
  slice_head(n = 25)

#prepare comparison data
# We load the metadata and ensure the mapping ID is also Uppercase
df5 <- read_excel("table_s5.xlsx", sheet = 2) %>% 
  clean_names() %>%
  slice(-1) %>% # Remove the row used for patient IDs
  select(proteins_metabolites, gene_symbol) %>%
  distinct() %>%
  mutate(proteins_metabolites = toupper(proteins_metabolites))

# merge
# Note: 'Protein_ID' from your list matches 'proteins_metabolites' from the metadata
top_25_named <- top_25_list %>% 
  left_join(df5, by = c("Protein_ID" = "proteins_metabolites"))

# --- RESULTS VISUALIZATION ---
# view results
print(top_25_named)

# plot
# If the gene symbol exists, use it. If not, fallback to the Protein_ID.
top_25_named <- top_25_named %>%
  mutate(Label = ifelse(!is.na(gene_symbol) & gene_symbol != "", gene_symbol, Protein_ID))

ggplot(top_25_named, aes(x = reorder(Label, Severe), y = Severe)) +
  
  # Draw the lines 
  geom_segment(aes(xend = Label, yend = 0), color = "grey50", linewidth = 1) +
  
  # Draw the points at the end 
  geom_point(color = "steelblue", size = 2) +
  
  coord_flip() +
  theme_minimal() +
  
  labs(
    title = "Exercise 2: Top 25 Protein Predictors",
    x = NULL,              
    y = "Importance"       
  ) +
  
  theme(
    panel.grid = element_blank(), 
    
    # Draw border
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black")
  )

# Volcano plot
library(dplyr)
library(ggplot2)
library(ggrepel)

#Pick the top 25 proteins from your RF importance table 
top_proteins <- tolower(top_25_named$Protein_ID)
# Choose the dataset that contains protein values + group labels 
df <- df_final_clean 
# Compute log2FC (Severe vs non_Severe) + p-values for those proteins -
volcano_df <- lapply(top_proteins, function(p) {
  
  if (!p %in% colnames(df)) return(NULL)
  
  x_sev <- df %>% filter(group_d == "Severe") %>% pull(!!sym(p))
  x_non <- df %>% filter(group_d == "non_Severe") %>% pull(!!sym(p))
  
  x_sev <- x_sev[is.finite(x_sev)]
  x_non <- x_non[is.finite(x_non)]
  
  if (length(x_sev) < 2 || length(x_non) < 2) return(NULL)
  
  # Robust effect size: median log2 fold change
  log2fc <- log2(median(x_sev) + 1e-9) - log2(median(x_non) + 1e-9)
  
  # Robust test for omics distributions
  pval <- suppressWarnings(wilcox.test(x_sev, x_non)$p.value)
  
  data.frame(Protein_ID = p, log2FC = log2fc, p_value = pval)
}) %>%
  bind_rows() %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    negLog10FDR = -log10(p_adj),
    Significant = p_adj < 0.05 & abs(log2FC) >= 1
  )

# Add gene symbols from mapping table 
volcano_df <- volcano_df %>%
  left_join(
    top_25_named %>%
      mutate(Protein_ID = tolower(Protein_ID)) %>%
      select(Protein_ID, gene_symbol),
    by = "Protein_ID"
  ) %>%
  mutate(Label = ifelse(!is.na(gene_symbol) & gene_symbol != "", gene_symbol, Protein_ID))

ggplot(volcano_df, aes(x = log2FC, y = negLog10FDR)) +
  geom_point(aes(shape = Significant), size = 3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  ggrepel::geom_text_repel(
    data = volcano_df %>% arrange(p_adj) %>% slice(1:10),
    aes(label = Label),
    max.overlaps = 30
  ) +
  labs(
    title = "Volcano plot (Top 25 RF proteins): Severe vs non-Severe",
    x = "log2 Fold Change (Severe / non-Severe)",
    y = "-log10(FDR)"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank()
  )

# --- HEATMAP OF SIGNIFICANT FEATURES ---
library(pheatmap)
library(dplyr)

# 1) Matrix (patients x proteins)
significant_proteins <- volcano_df %>%
  filter(p_adj < 0.05 & abs(log2FC) >= 1) %>%
  pull(Protein_ID)

heatmap_data <- training_data %>%
  select(all_of(significant_proteins), group_d)

expr_matrix <- heatmap_data %>%
  select(-group_d) %>%
  as.matrix()

# 2) Give the matrix valid rownames (must match annotation_row rownames)
rownames(expr_matrix) <- paste0("S", seq_len(nrow(expr_matrix)))

# 3) Scale by protein (column-wise) -> keep same rownames
expr_matrix_scaled <- scale(expr_matrix)

# 4) Annotation with matching rownames
matched_labels <- volcano_df$Label[match(colnames(expr_matrix_scaled), volcano_df$Protein_ID)]
colnames(expr_matrix_scaled) <- matched_labels
annotation_df <- data.frame(Group = heatmap_data$group_d)
rownames(annotation_df) <- rownames(expr_matrix_scaled)

# 5) Plot
pheatmap(
  expr_matrix_scaled,
  annotation_row = annotation_df,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  main = "Heatmap of Significant Proteins"
)

# --- BIOLOGICAL INTERPRETATION ---
# Make sure Protein_ID in volcano_df matches case
volcano_df$Protein_ID <- tolower(volcano_df$Protein_ID)
top_25_named$Protein_ID <- tolower(top_25_named$Protein_ID)

feature_mapping <- top_25_named %>%
  left_join(
    volcano_df %>%
      select(Protein_ID, log2FC, p_adj),
    by = "Protein_ID"
  ) %>%
  arrange(desc(Severe))   # Severe column = RF importance score

print(feature_mapping)


# ---------------------------------------------------------------
# --- EXERCISE 3 ------------------------------------------------
# ---------------------------------------------------------------
# Independent test cohort preparation
df_test_raw <- read_excel("table_s4.xlsx", sheet = 2) %>%
  clean_names()
#Comparison based on supplementary table 4
df1_clean_2 <- df1_clean %>%
  mutate(
    group = case_when(
      group_d %in% c(0, 1) ~ "Non-COVID-19",
      group_d == 2 ~ "Non-severe",
      group_d == 3 ~ "Severe",
      TRUE ~ NA_character_
    ),
    group = factor(group, levels = c("Non-severe", "Severe"))
  ) %>%
  select(patient_id, group) %>%
  filter(group %in% c("Non-severe", "Severe"))

# --- Data reshaping and formatting ---
new_headers <- df_test_raw[1, ] %>%
  select(patient_id:last_col()) %>% 
  unlist() %>%
  as.character()
colnames(df_test_raw)[3:ncol(df_test_raw)] <- new_headers

#pivot to tidy format
df_test_raw_final <- df_test_raw %>%
  slice(-1) %>%                     
  pivot_longer(
    cols = all_of(new_headers),         
    names_to = "patient_id",         
    values_to = "expression"           
  ) %>%
  select(-gene_symbol) %>%
  pivot_wider(
    names_from = proteins_metabolites,          
    values_from = expression
  ) %>%
  mutate(across(-patient_id, as.numeric))

df_test_with_labels <- df_test_raw_final %>%
  inner_join(df1_clean_2, by = "patient_id") %>%
  mutate(group = recode(group, "Non-severe" = "non_Severe")) %>%
  mutate(group = factor(group, levels = c("non_Severe", "Severe"))) %>%
  clean_names()

# --- Feature alignment with training model ---
# Extracting model predictors
model_predictors <- rf_model_2$finalModel$xNames

# Handling missing features
missing_features <- setdiff(model_predictors, colnames(df_test_with_labels))
# missing_features

df_test_with_labels[missing_features] <- NA #Add the missing columns to the test set, filled with NAs
#This is done because RF requires the exact same column names in the exact order as the training data.

x_test <- df_test_with_labels %>%
  dplyr::select(all_of(model_predictors))

# --- Model application to independent cohort ---
# Generating class predictions and predicted probabilities
# Predict disease severity
# use the rf_model_2 to see if it can correctly guess which proteins are severe
# predictions (make sure these are vectors, not data frames) --> Apply trained model to the test cohort
test_predictions <- predict(rf_model_2, newdata = x_test)
test_probabilities <- predict(rf_model_2, newdata = x_test, type = "prob")

# If predict() returned a data frame (tidymodels), extract the column:
if (is.data.frame(test_predictions)) test_predictions <- test_predictions$.pred

# Probabilities column name can differ: "Severe" (caret) vs ".pred_Severe" (tidymodels)
prob_severe <- if (is.data.frame(test_probabilities)) {
  if ("Severe" %in% colnames(test_probabilities)) test_probabilities$Severe else test_probabilities$.pred_Severe
} else {
  test_probabilities[, "Severe"]
}
# Results table
test_results_table <- tibble(
  patient_id = df_test_with_labels$patient_id,
  Actual = df_test_with_labels$group,
  Predicted_Class = test_predictions,
  Prob_Severe = prob_severe
) %>%
  mutate(Correct = ifelse(Actual == Predicted_Class, "Correct", "Wrong"))

print(test_results_table)

# Calibration plot and Brierer score
# Convert actual class to numeric (1 = Severe, 0 = non_Severe)
calibration_df <- test_results_table %>%
  mutate(
    Actual_numeric = ifelse(Actual == "Severe", 1, 0)
  )

# Fit logistic regression of observed outcome vs predicted probability
calibration_model <- glm(Actual_numeric ~ Prob_Severe,
                         data = calibration_df,
                         family = binomial)

# Create smooth calibration curve
prob_seq <- seq(0, 1, length.out = 100)

calibration_curve <- data.frame(
  Prob_Severe = prob_seq,
  Observed = predict(calibration_model,
                     newdata = data.frame(Prob_Severe = prob_seq),
                     type = "response")
)

# Plot
ggplot(calibration_df, aes(x = Prob_Severe, y = Actual_numeric)) +
  geom_point(size = 3) +
  geom_line(data = calibration_curve,
            aes(x = Prob_Severe, y = Observed),
            color = "blue", linewidth = 1) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed",
              color = "red") +
  labs(
    title = "Calibration Plot – External Cohort",
    x = "Predicted Probability (Severe)",
    y = "Observed Proportion (Severe)"
  ) +
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

brier_score <- mean((calibration_df$Prob_Severe - calibration_df$Actual_numeric)^2)
brier_score

# --- Performance evaluation on external cohort ---
# Confusion matrix (predicted vs truth)
final_metrics <- confusionMatrix(
  factor(test_results_table$Predicted_Class, levels = levels(test_results_table$Actual)),
  test_results_table$Actual,
  positive = "Severe"
)

print(final_metrics)
