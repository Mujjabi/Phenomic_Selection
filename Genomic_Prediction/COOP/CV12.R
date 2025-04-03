
library(dplyr)
library(sommer)
library(splitTools)


# SCA relationship matrix (Hadamard product)
inbred_names <- New_crosses$Parent1
tester_names <- New_crosses$Tester

# Ensure inbred and tester names exist in their respective matrices
valid_hybrids <- inbred_names %in% rownames(GRM_New_inv) & tester_names %in% rownames(Ginv.blend_tester)
hybrid_list <- paste(New_crosses$Parent1, New_crosses$Tester, sep = ":")
hybrid_list <- unique(hybrid_list)
# Filter only valid crosses
inbred_names <- inbred_names[valid_hybrids]
tester_names <- tester_names[valid_hybrids]
hybrid_list <- hybrid_list[valid_hybrids]  # Keep only valid hybrids

# Construct G_hybrid only for specified crosses
G_hybrid <- matrix(0, nrow = length(hybrid_list), ncol = length(hybrid_list))
rownames(G_hybrid) <- hybrid_list
colnames(G_hybrid) <- hybrid_list

for (i in seq_along(hybrid_list)) {
  for (j in seq_along(hybrid_list)) {
    G_hybrid[i, j] <- 0.5 * (Ginv.blend[inbred_names[i], inbred_names[j]] + 
                               Ginv.blend_tester[tester_names[i], tester_names[j]])
  }
}

# Ensure correct alignment with dataset

# Get the order of hybrids based on Parent1_Tester
hybrid_order <- match(unique(Pheno_org$Parent1_Tester), rownames(G_hybrid))

# Remove any NA values (in case some hybrids are missing from G_hybrid)
hybrid_order <- hybrid_order[!is.na(hybrid_order)]

# Reorder both rows and columns of G_hybrid
G_hybrid <- G_hybrid[hybrid_order, hybrid_order]



#Set folds
num_folds <- 5 
n_iterations <- 1


## CV1O, Training on conventional to predict Organic 

create_cv_folds_scheme1O <- function(data, num_folds) {
  #set.seed(6343)  # For reproducibility
  
  unique_hybrids <- unique(data$Parent1_Tester)
  
  # Assign hybrids to folds
  fold_assignments <- tibble(Parent1_Tester = unique_hybrids) %>%
    mutate(Fold = sample(rep(1:num_folds, length.out = n()), replace = FALSE))
  
  # Merge fold assignments into dataset
  data <- data %>%
    left_join(fold_assignments, by = "Parent1_Tester")
  
  train_data_list <- list()
  test_data_list <- list()
  
  for (i in 1:num_folds) {
    # Get hybrids assigned to the current fold
    test_hybrids <- fold_assignments %>% filter(Fold == i) %>% pull(Parent1_Tester)
    
    # Split data for the current fold
    # Training set: Include all data from one system for the test hybrids
    # Validation set: Include data from the other system for the test hybrids
    train_data <- data %>%
      filter((Parent1_Tester %in% test_hybrids & System == "Conventional") |  # Include Conventional for test hybrids
               !(Parent1_Tester %in% test_hybrids))  # Include all data for non-test hybrids
    
    test_data <- data %>%
      filter(Parent1_Tester %in% test_hybrids & System == "Organic")  # Include Organic for test hybrids
    
    # Store training and validation data for the current fold
    train_data_list[[i]] <- train_data
    test_data_list[[i]] <- test_data
  }
  
  return(list(train_data_list = train_data_list, test_data_list = test_data_list))
}


# Perform K-Fold Cross-Validation with CV1
cv_scheme_1O <- create_cv_folds_scheme1O(Pheno_org, num_folds = 5)

## CV Scheme 1
accuracy_results_scheme1O <- data.frame(
  Iteration = integer(),
  Fold = integer(),
  PSM1 = numeric(),
  PSM2 = numeric(),
  PSM3 = numeric(),
  SPM1 = numeric(),
  SPM2 = numeric(),
  SPM3 = numeric()
)

# Apply CV Schemes to the model
for (iter in 1:n_iterations) {
  for (i in 1:num_folds) {
    # Extract training and test data for the current fold
    training_data   <- cv_scheme_1O$train_data_list[[i]]
    validation_data <- cv_scheme_1O$test_data_list[[i]]
    
    # Fit the sommer model on the training data
    model <- mmer(
      BLUE ~ 1 + System,  # Fixed effects
      random = ~ vsr(Parent1, Gu = Ginv.blend ) +  # GCA for Parent1
        vsr(Tester, Gu = Ginv.blend_tester) +   # GCA for Tester
        vsr(Parent1_Tester, Gu = G_hybrid_inv), #SCA effect
      rcov = ~ diag(System),  # Allows different residual variances for each system
      weights = Weights,
      data = training_data)
    
    # Extract genetic effects
    GCA_female <- model$U$`u:Parent1`$BLUE
    GCA_male   <- model$U$`u:Tester`$BLUE
    SCA        <- model$U$`u:Parent1_Tester`$BLUE

    
    # Split validation data into Organic and Conventional
    Predicted_Org  <- validation_data %>% filter(System == "Organic")
 
    # Predictions for Organic system
    org_intercept         <- model$Beta$Estimate[1] + model$Beta$Estimate[2]  #conventional is the reference system
    Predicted_Org$Model1  <- org_intercept + GCA_female[Predicted_Org$Parent1]
    Predicted_Org$Model2  <- org_intercept + GCA_female[Predicted_Org$Parent1] + GCA_male[Predicted_Org$Tester]
    Predicted_Org$Model3  <- org_intercept + GCA_female[Predicted_Org$Parent1] + GCA_male[Predicted_Org$Tester] + SCA[Predicted_Org$Parent1_Tester]
    
    
    # Combine organic and conventional predictions
    #Predicted_Org <- org
    
    # Compute Pearson and Spearman correlations for model accuracy
    Pearson_Model1  <- cor(Predicted_Org$BLUP, Predicted_Org$Model1, use = "complete.obs")
    Pearson_Model2  <- cor(Predicted_Org$BLUP, Predicted_Org$Model2, use = "complete.obs")
    Pearson_Model3  <- cor(Predicted_Org$BLUP, Predicted_Org$Model3, use = "complete.obs")
  
    
    spearman_Model1 <- cor(Predicted_Org$BLUP, Predicted_Org$Model1, use = "complete.obs", method = "spearman")
    spearman_Model2 <- cor(Predicted_Org$BLUP, Predicted_Org$Model2, use = "complete.obs", method = "spearman")
    spearman_Model3 <- cor(Predicted_Org$BLUP, Predicted_Org$Model3, use = "complete.obs", method = "spearman")
    
      # Store accuracy metrics for this fold and iteration
    accuracy_results_scheme1O  <- rbind(accuracy_results_scheme1O, data.frame(
      Iteration = iter,
      Fold = i,
      PSM1  = Pearson_Model1, 
      PSM2  = Pearson_Model2,
      PSM3  = Pearson_Model3,
      SPM1  = spearman_Model1, 
      SPM2  = spearman_Model2, 
      SPM3  = spearman_Model3
    ))
  }
}


# Step 5: Plot Accuracies
library(tidyr)
library(ggplot2)

# Reshape accuracy_results for plotting
PA_results_long <-    accuracy_results_scheme1O %>%  
  pivot_longer(cols = c(PSM1,PSM2,PSM3,
                        SPM1,SPM2,SPM3,), 
               names_to = "Method", 
               values_to = "Accuracy")


# Plot boxplots of accuracy metrics
ggplot(PA_results_long, aes(x = Method, y = Accuracy, fill = Method)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot(outlier.colour = "Black", outlier.shape = 1, outlier.alpha = 0, # Remove outliers
               outlier.size = 2, notch = FALSE, varwidth = FALSE) +
  labs(x = "Prediction Models", y = "Accuracy", title = "Predictive Ability for Grain Yield") +
  theme_minimal() +
  theme(
    legend.position = "right",  # Place legend at the bottom
    strip.text = element_text(size = 14, face = "bold"),  # Increase facet label font size
    axis.title.x = element_text(size = 14, face = "bold"),  # Increase x-axis title font size
    axis.title.y = element_text(size = 14, face = "bold"),  # Increase y-axis title font size
    axis.text.x = element_text(size = 12),  # Increase x-axis tick label font size
    axis.text.y = element_text(size = 12),
    panel.grid.minor = element_blank()) + 
  theme(legend.position = "none") +  # Rem
  scale_y_continuous(limits = c(0.1, 1), breaks = seq(0.1, 1, by = 0.2)) +
  #scale_fill_manual(values = c("Pearson" = "#1f78b4", "Spearman" = "#e31a1c")) +
  scale_color_manual(values = c("#a6cee3","#33a02c","#1f78b4","#b2df8a","#fb9a99",
                                "#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a", "yellow"))

## CV1C, Training on Organic to predict conventional

create_cv_folds_scheme1C <- function(data, num_folds) {
  #set.seed(6343)  # For reproducibility
  
  unique_hybrids <- unique(data$Parent1_Tester)
  
  # Assign hybrids to folds
  fold_assignments <- tibble(Parent1_Tester = unique_hybrids) %>%
    mutate(Fold = sample(rep(1:num_folds, length.out = n()), replace = FALSE))
  
  # Merge fold assignments into dataset
  data <- data %>%
    left_join(fold_assignments, by = "Parent1_Tester")
  
  train_data_list <- list()
  test_data_list <- list()
  
  for (i in 1:num_folds) {
    # Get hybrids assigned to the current fold
    test_hybrids <- fold_assignments %>% filter(Fold == i) %>% pull(Parent1_Tester)
    
    # Split data for the current fold
    # Training set: Include all data from the Conventional system for the test hybrids
    # Validation set: Include data from the Organic system for the test hybrids
    train_data <- data %>%
      filter((Parent1_Tester %in% test_hybrids & System == "Organic") |  # Include Conventional for test hybrids
               !(Parent1_Tester %in% test_hybrids))  # Include all data for non-test hybrids
    
    test_data <- data %>%
      filter(Parent1_Tester %in% test_hybrids & System == "Conventional")  # Include Conventional for test hybrids
    
    # Store training and validation data for the current fold
    train_data_list[[i]] <- train_data
    test_data_list[[i]] <- test_data
  }
  
  return(list(train_data_list = train_data_list, test_data_list = test_data_list))
}


# Perform K-Fold Cross-Validation with Scheme 1C (Train on Conventional, Validate on Organic)
cv_scheme_1C <- create_cv_folds_scheme1C(Pheno_org, num_folds = 5)

## CV Scheme 1C
accuracy_results_scheme1C <- data.frame(
  Iteration = integer(),
  Fold = integer(),
  PSM1 = numeric(),
  PSM2 = numeric(),
  PSM3 = numeric(),
  SPM1 = numeric(),
  SPM2 = numeric(),
  SPM3 = numeric()
)

## Apply CV1C to model
for (iter in 1:n_iterations) {
  for (i in 1:num_folds) {
    # Extract training and test data for the current fold
    training_data <- cv_scheme_1C$train_data_list[[i]]
    validation_data <- cv_scheme_1C$test_data_list[[i]]
    
    # Fit the sommer model on the training data
    model <- mmer(
      BLUE ~ 1 + System,  # Fixed effects
      random = ~ vsr(Parent1, Gu = Ginv.blend ) +  # GCA for Parent1
        vsr(Tester, Gu = Ginv.blend_tester) +   # GCA for Tester
        vsr(Parent1_Tester, Gu = G_hybrid_inv), #Ginv.blend_Ihybrid), #SCA effect
      rcov = ~ diag(System),  # Allows different residual variances for each system
      weights = Weights,
      data = training_data)
    
    # Extract genetic effects
    GCA_female <- model$U$`u:Parent1`$BLUE
    GCA_male   <- model$U$`u:Tester`$BLUE
    SCA        <- model$U$`u:Parent1_Tester`$BLUE
    
    # Split validation data into Organic and Conventional
    Predicted_conv  <- validation_data %>% filter(System == "Conventional")
    
    # Predictions for Conventional system
    Predicted_conv$Model1  <-   model$Beta$Estimate[1]+ GCA_female[Predicted_conv$Parent1]
    Predicted_conv$Model2  <-   model$Beta$Estimate[1]+ GCA_female[Predicted_conv$Parent1] + GCA_male[Predicted_conv$Tester]
    Predicted_conv$Model3  <-   model$Beta$Estimate[1]+ GCA_female[Predicted_conv$Parent1] + GCA_male[Predicted_conv$Tester] + SCA[Predicted_conv$Parent1_Tester]
    
    # Compute Pearson and Spearman correlations for model accuracy
    Pearson_Model1  <- cor(Predicted_conv$BLUP, Predicted_conv$Model1, use = "complete.obs")
    Pearson_Model2  <- cor(Predicted_conv$BLUP, Predicted_conv$Model2, use = "complete.obs")
    Pearson_Model3  <- cor(Predicted_conv$BLUP, Predicted_conv$Model3, use = "complete.obs")
    spearman_Model1 <- cor(Predicted_conv$BLUP, Predicted_conv$Model1, use = "complete.obs", method = "spearman")
    spearman_Model2 <- cor(Predicted_conv$BLUP, Predicted_conv$Model2, use = "complete.obs", method = "spearman")
    spearman_Model3 <- cor(Predicted_conv$BLUP, Predicted_conv$Model3, use = "complete.obs", method = "spearman")
    
    # Store accuracy metrics for this fold and iteration
    accuracy_results_scheme1C  <- rbind(accuracy_results_scheme1C, data.frame(
      Iteration = iter,
      Fold = i,
      PSM1  = Pearson_Model1, 
      PSM2  = Pearson_Model2,
      PSM3  = Pearson_Model3,
      SPM1  = spearman_Model1, 
      SPM2  = spearman_Model2, 
      SPM3  = spearman_Model3
    ))
  }
}

# Reshape accuracy_results for plotting accuracies

PA_results_long2 <-    accuracy_results_scheme1C  %>%  
  pivot_longer(cols = c(PSM1,PSM2,PSM3,
                        SPM1,SPM2,SPM3,), 
               names_to = "Method", 
               values_to = "Accuracy")


# Plot boxplots of accuracy metrics
ggplot(PA_results_long2, aes(x = Method, y = Accuracy, fill = Method)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot(outlier.colour = "Black", outlier.shape = 1, outlier.alpha = 0, # Remove outliers
               outlier.size = 2, notch = FALSE, varwidth = FALSE) +
  labs(x = "Prediction Models", y = "Accuracy", title = "Predictive Ability for Grain Yield") +
  theme_minimal() +
  theme(
    legend.position = "right",  # Place legend at the bottom
    strip.text = element_text(size = 14, face = "bold"),  # Increase facet label font size
    axis.title.x = element_text(size = 14, face = "bold"),  # Increase x-axis title font size
    axis.title.y = element_text(size = 14, face = "bold"),  # Increase y-axis title font size
    axis.text.x = element_text(size = 12),  # Increase x-axis tick label font size
    axis.text.y = element_text(size = 12),
    panel.grid.minor = element_blank()) + 
  theme(legend.position = "none") +  # Rem
  scale_y_continuous(limits = c(0.1, 1), breaks = seq(0.1, 1, by = 0.2)) +
  #scale_fill_manual(values = c("Pearson" = "#1f78b4", "Spearman" = "#e31a1c")) +
  scale_color_manual(values = c("#a6cee3","#33a02c","#1f78b4","#b2df8a","#fb9a99",
                                "#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a", "yellow"))


##Combine the organic and conv to compare predictive ability
Predicted_combine <- rbind(Predicted_Org, Predicted_conv)

# Compute correlation for each system
cor_results <- Predicted_combine %>%
  group_by(System) %>%
  summarise(R_value = cor(BLUP, Model3, use = "complete.obs"),
            RMSE = sqrt(mean((BLUP - Model3)^2, na.rm = TRUE)))

print(cor_results)  # View correlation results

# Create scatter plot with regression lines and display R-values
ggplot(Predicted_combine , aes(x = Model3, y = BLUP, color = System)) +
  geom_point(size = 2, alpha = 0.7) +  # Scatter points
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Regression line
  facet_wrap(~ System) +  # Separate plots for each system
  theme_minimal() +
  labs(title = "Predicted vs. Actual BLUPs Across Systems",
       x = "Predicted Value",
       y = "Actual BLUP") +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 14, face = "bold"), # Style facet labels
        axis.title.x = element_text(size = 14),               # X-axis label size
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 14),         # Legend title size
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size = 12),  # X-axis values size
        axis.text.y = element_text(size = 12)) +  # Y-axis values size
  geom_text(data = cor_results, 
            aes(x = min(Predicted_combine$Model3), 
                y = max(Predicted_combine$BLUP), 
                    label = paste("R =", round(R_value, 2), "\nRMSE =", round(RMSE, 2))),
                    hjust = 0, vjust = 1, size = 4, color = "black")



## CV2,Hybrids in validation set are seen before in the training set for both systems
create_cv_folds_scheme2 <- function(data, num_folds) {
  #set.seed(343)  # For reproducibility
  
  unique_hybrids <- unique(data$Parent1_Tester)
  
  # Assign hybrids to folds
  fold_assignments <- tibble(Parent1_Tester = unique_hybrids) %>%
    mutate(Fold = sample(rep(1:num_folds, length.out = n()), replace = FALSE))
  
  # Merge fold assignments into dataset
  data <- data %>%
    left_join(fold_assignments, by = "Parent1_Tester")
  
  train_data_list <- list()
  test_data_list <- list()
  
  for (i in 1:num_folds) {
    # Get hybrids assigned to the current fold
    test_hybrids <- fold_assignments %>% filter(Fold == i) %>% pull(Parent1_Tester)
    
    # Training set: Include data from both systems for hybrids not in the current fold
    train_data <- data %>%
      filter(!(Parent1_Tester %in% test_hybrids))  # Exclude hybrids in the current fold
    
    # Validation set: Include data from one system for hybrids in the current fold
    # Here, we use the Organic system for validation (can be swapped to Conventional if needed)
    test_data <- data %>%
      filter(Parent1_Tester %in% test_hybrids & System == c("Organic","Conventional"))  # Use both
    
    # Store training and validation data for the current fold
    train_data_list[[i]] <- train_data
    test_data_list[[i]] <- test_data
  }
  
  return(list(train_data_list = train_data_list, test_data_list = test_data_list))
}


# Perform K-Fold Cross-Validation with Scheme 2 (CV2)
cv_scheme_2 <- create_cv_folds_scheme2(Pheno_org, num_folds = 5)

## CV Scheme 2
accuracy_results_scheme2 <- data.frame(
  Iteration = integer(),
  Fold = integer(),
  PSM1 = numeric(),
  PSM2 = numeric(),
  PSM3 = numeric(),
  SPM1 = numeric(),
  SPM2 = numeric(),
  SPM3 = numeric()
)


# Apply CV Schemes to the model
for (iter in 1:n_iterations) {
  for (i in 1:num_folds) {
    # Extract training and test data for the current fold
    training_data <- cv_scheme_2$train_data_list[[i]]
    validation_data <- cv_scheme_2$test_data_list[[i]]
    
    # Fit the sommer model on the training data
    model <- mmer(
      BLUE ~ 1 + System,  # Fixed effects
      random = ~ vsr(Parent1, Gu =Ginv.blend) +  # GCA for Parent1
        vsr(Tester, Gu = Ginv.blend_tester) +   # GCA for Tester
        vsr(Parent1_Tester, Gu = G_hybrid_inv), #SCA effect
      rcov = ~ diag(System),  # Allows different residual variances for each system
      weights = Weights,
      data = Pheno_org)
    
    # Extract genetic effects
    GCA_female <- model$U$`u:Parent1`$BLUE
    GCA_male   <- model$U$`u:Tester`$BLUE
    SCA        <- model$U$`u:Parent1_Tester`$BLUE

    
    # Split validation data into Organic and Conventional
    org   <- validation_data %>% filter(System == "Organic")
    conv  <- validation_data %>% filter(System == "Conventional")
    
    # Predictions for Conventional system
    conv$Model1  <-   model$Beta$Estimate[1]+ GCA_female[conv$Parent1]
    conv$Model2  <-   model$Beta$Estimate[1]+ GCA_female[conv$Parent1] + GCA_male[conv$Tester]
    conv$Model3  <-   model$Beta$Estimate[1]+ GCA_female[conv$Parent1] + GCA_male[conv$Tester] + SCA[conv$Parent1_Tester]

    
    # Predictions for Organic system
    org_intercept   <- model$Beta$Estimate[1] + model$Beta$Estimate[2]  #conventional is the reference system
    org$Model1  <- org_intercept + GCA_female[org$Parent1]
    org$Model2  <- org_intercept + GCA_female[org$Parent1] + GCA_male[org$Tester]
    org$Model3  <- org_intercept + GCA_female[org$Parent1] + GCA_male[org$Tester] + SCA[org$Parent1_Tester]
    
# Combine organic and conventional predictions
    Predicted_cv2 <- rbind(org, conv)
    
# Compute Pearson and Spearman correlations for model accuracy
    Pearson_Model1  <- cor(Predicted_cv2$BLUP, Predicted_cv2$Model1, use = "complete.obs")
    Pearson_Model2  <- cor(Predicted_cv2$BLUP, Predicted_cv2$Model2, use = "complete.obs")
    Pearson_Model3  <- cor(Predicted_cv2$BLUP, Predicted_cv2$Model3, use = "complete.obs")
    
    spearman_Model1 <- cor(Predicted_cv2$BLUP, Predicted_cv2$Model1, use = "complete.obs", method = "spearman")
    spearman_Model2 <- cor(Predicted_cv2$BLUP, Predicted_cv2$Model2, use = "complete.obs", method = "spearman")
    spearman_Model3 <- cor(Predicted_cv2$BLUP, Predicted_cv2$Model3, use = "complete.obs", method = "spearman")
    
# Store accuracy metrics for this fold and iteration
    accuracy_results_scheme2  <- rbind(accuracy_results_scheme2, data.frame(
      Iteration = iter,
      Fold = i,
      PSM1  = Pearson_Model1, 
      PSM2  = Pearson_Model2,
      PSM3  = Pearson_Model3,
      SPM1  = spearman_Model1, 
      SPM2  = spearman_Model2, 
      SPM3  = spearman_Model3
    ))
  }
}

##Reshape to long format

PA_results_long_cv2 <-   accuracy_results_scheme2 %>%  
                         pivot_longer(cols = c(PSM1,PSM2,PSM3,
                        SPM1,SPM2,SPM3,), 
                        names_to = "Method", 
                        values_to = "Accuracy")


# Plot boxplots of accuracy metrics
ggplot(PA_results_long_cv2, aes(x = Method, y = Accuracy, fill = Method)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot(outlier.colour = "Black", outlier.shape = 1, outlier.alpha = 0, # Remove outliers
               outlier.size = 2, notch = FALSE, varwidth = FALSE) +
  labs(x = "Prediction Models", y = "Accuracy", title = "Predictive Ability for Grain Yield") +
  theme_minimal() +
  theme(
    legend.position = "right",  # Place legend at the bottom
    strip.text = element_text(size = 14, face = "bold"),  # Increase facet label font size
    axis.title.x = element_text(size = 14, face = "bold"),  # Increase x-axis title font size
    axis.title.y = element_text(size = 14, face = "bold"),  # Increase y-axis title font size
    axis.text.x = element_text(size = 12),  # Increase x-axis tick label font size
    axis.text.y = element_text(size = 12),
    panel.grid.minor = element_blank()) + 
  theme(legend.position = "none") +  # Rem
  scale_y_continuous(limits = c(0.1, 1), breaks = seq(0.1, 1, by = 0.2)) +
  #scale_fill_manual(values = c("Pearson" = "#1f78b4", "Spearman" = "#e31a1c")) +
  scale_color_manual(values = c("#a6cee3","#33a02c","#1f78b4","#b2df8a","#fb9a99",
                                "#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a", "yellow"))

# Compute correlation for each system
cor_results <-   Predicted_cv2 %>%
                 group_by(System) %>%
                 summarise(R_value = cor(BLUE, Model3, use = "complete.obs"),
                 RMSE = sqrt(mean((BLUE - Model3)^2, na.rm = TRUE)))

print(cor_results)  # View correlation results

# Create scatter plot with regression lines and display R-values
ggplot(Predicted_cv2 , aes(x = Model3, y = BLUE, color = System)) +
  geom_point(size = 2, alpha = 0.7) +  # Scatter points
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Regression line
  facet_wrap(~ System) +  # Separate plots for each system
  theme_minimal() +
  labs(title = "Predicted vs. Actual BLUPs Across Systems",
       x = "Predicted Value",
       y = "Actual BLUP") +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 14, face = "bold"), # Style facet labels
        axis.title.x = element_text(size = 14),               # X-axis label size
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 14),         # Legend title size
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size = 12),  # X-axis values size
        axis.text.y = element_text(size = 12)) +  # Y-axis values size
  geom_text(data = cor_results, 
            aes(x = min(Predicted_cv2$Model3, na.rm = TRUE), 
                y = max(Predicted_cv2$BLUE, na.rm = TRUE), 
                label = paste("R =", round(R_value, 2), "\nRMSE =", round(RMSE, 2))),
            hjust = 0, vjust = 1, size = 4, color = "black")



## CV0,Hybrids in validation set are never seen before in the training set for both systems

create_cv_folds_scheme0 <- function(data, num_folds) {
set.seed(343)  # For reproducibility
  
unique_hybrids <- unique(data$Parent1_Tester)
  
# Assign hybrids to folds
  fold_assignments <- tibble(Parent1_Tester = unique_hybrids) %>%
    mutate(Fold = sample(rep(1:num_folds, length.out = n()), replace = FALSE))
  
  # Merge fold assignments into dataset
  data <- data %>%
    left_join(fold_assignments, by = "Parent1_Tester")
  
  train_data_list <- list()
  test_data_list <- list()
  
  for (i in 1:num_folds) {
    # Get hybrids assigned to the current fold
    test_hybrids <- fold_assignments %>% filter(Fold == i) %>% pull(Parent1_Tester)
    
    # Training set: Include data from both systems for hybrids not in the current fold
    train_data <- data %>%
      filter(!(Parent1_Tester %in% test_hybrids))  # Exclude hybrids in the current fold
    
    # Validation set: Include data from one system for hybrids in the current fold
    # Here, we use the Organic system for validation (can be swapped to Conventional if needed)
    test_data <- data %>%
      filter(Parent1_Tester %in% test_hybrids & System == c("Organic","Conventional"))  # Use both
    
    # Store training and validation data for the current fold
    train_data_list[[i]] <- train_data
    test_data_list[[i]] <- test_data
  }
  
  return(list(train_data_list = train_data_list, test_data_list = test_data_list))
}


# Perform K-Fold Cross-Validation with Scheme 0
cv_scheme_0 <- create_cv_folds_scheme0(Pheno_org, num_folds = 5)

## CV Scheme 0
accuracy_results_scheme0 <- data.frame(
  Iteration = integer(),
  Fold = integer(),
  PSM1 = numeric(),
  PSM2 = numeric(),
  PSM3 = numeric(),
  SPM1 = numeric(),
  SPM2 = numeric(),
  SPM3 = numeric()
)

# Apply CV Schemes to the model
for (iter in 1:n_iterations) {
  for (i in 1:num_folds) {
    # Extract training and test data for the current fold
    training_data0   <- cv_scheme_0$train_data_list[[i]]
    validation_data0 <- cv_scheme_0$test_data_list[[i]]
    
    # Fit the sommer model on the training data
    model0 <- mmer(
      BLUE ~ 1 + System,  # Fixed effects
      random = ~ vsr(Parent1, Gu =Ginv.blend) +  # GCA for Parent1
        vsr(Tester, Gu = Ginv.blend_tester) +   # GCA for Tester
        vsr(Parent1_Tester, Gu = Ginv.blend_Ihybrid), #SCA effect
      rcov = ~ diag(System),  # Allows different residual variances for each system
      weights = Weights,
      data = training_data0)
    
    # Extract genetic effects
    GCA_female <- model0$U$`u:Parent1`$BLUE
    GCA_male   <- model0$U$`u:Tester`$BLUE
    SCA        <- model0$U$`u:Parent1_Tester`$BLUE
    
    
    # Split validation data into Organic and Conventional
    org   <- validation_data0 %>% filter(System == "Organic")
    conv  <- validation_data0 %>% filter(System == "Conventional")
    
    # Predictions for Conventional system
    conv$Model1  <-   model0$Beta$Estimate[1]+ GCA_female[conv$Parent1]
    conv$Model2  <-   model0$Beta$Estimate[1]+ GCA_female[conv$Parent1] + GCA_male[conv$Tester]
    conv$Model3  <-   model0$Beta$Estimate[1]+ GCA_female[conv$Parent1] + GCA_male[conv$Tester] + SCA[conv$Parent1_Tester]
    
    
    # Predictions for Organic system
    org_intercept   <- model0$Beta$Estimate[1] + model0$Beta$Estimate[2]  #conventional is the reference system
    org$Model1  <- org_intercept + GCA_female[org$Parent1]
    org$Model2  <- org_intercept + GCA_female[org$Parent1] + GCA_male[org$Tester]
    org$Model3  <- org_intercept + GCA_female[org$Parent1] + GCA_male[org$Tester] + SCA[org$Parent1_Tester]
    
    # Combine organic and conventional predictions
    Predicted_cv0 <- rbind(org, conv)
    
    # Compute Pearson and Spearman correlations for model accuracy
    Pearson_Model1  <- cor(Predicted_cv0$BLUE, Predicted_cv0$Model1, use = "complete.obs")
    Pearson_Model2  <- cor(Predicted_cv0$BLUE, Predicted_cv0$Model2, use = "complete.obs")
    Pearson_Model3  <- cor(Predicted_cv0$BLUE, Predicted_cv0$Model3, use = "complete.obs")
    spearman_Model1 <- cor(Predicted_cv0$BLUE, Predicted_cv0$Model1, use = "complete.obs", method = "spearman")
    spearman_Model2 <- cor(Predicted_cv0$BLUE, Predicted_cv0$Model2, use = "complete.obs", method = "spearman")
    spearman_Model3 <- cor(Predicted_cv0$BLUE, Predicted_cv0$Model3, use = "complete.obs", method = "spearman")
    
    # Store accuracy metrics for this fold and iteration
    accuracy_results_scheme0  <- rbind(accuracy_results_scheme0, data.frame(
      Iteration = iter,
      Fold = i,
      PSM1  = Pearson_Model1, 
      PSM2  = Pearson_Model2,
      PSM3  = Pearson_Model3,
      SPM1  = spearman_Model1, 
      SPM2  = spearman_Model2, 
      SPM3  = spearman_Model3
    ))
  }
}

##Reshape to long format

PA_results_long_cv0 <-   accuracy_results_scheme0 %>%  
  pivot_longer(cols = c(PSM1,PSM2,PSM3,
                        SPM1,SPM2,SPM3,), 
               names_to = "Method", 
               values_to = "Accuracy")


# Plot boxplots of accuracy metrics
ggplot(PA_results_long_cv0, aes(x = Method, y = Accuracy, fill = Method)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot(outlier.colour = "Black", outlier.shape = 1, outlier.alpha = 0, # Remove outliers
               outlier.size = 2, notch = FALSE, varwidth = FALSE) +
  labs(x = "Prediction Models", y = "Accuracy", title = "Predictive Ability for Grain Yield") +
  theme_minimal() +
  theme(
    legend.position = "right",  # Place legend at the bottom
    strip.text = element_text(size = 14, face = "bold"),  # Increase facet label font size
    axis.title.x = element_text(size = 14, face = "bold"),  # Increase x-axis title font size
    axis.title.y = element_text(size = 14, face = "bold"),  # Increase y-axis title font size
    axis.text.x = element_text(size = 12),  # Increase x-axis tick label font size
    axis.text.y = element_text(size = 12),
    panel.grid.minor = element_blank()) + 
  theme(legend.position = "none") +  # Rem
  scale_y_continuous(limits = c(0.1, 1), breaks = seq(0.1, 1, by = 0.2)) +
  #scale_fill_manual(values = c("Pearson" = "#1f78b4", "Spearman" = "#e31a1c")) +
  scale_color_manual(values = c("#a6cee3","#33a02c","#1f78b4","#b2df8a","#fb9a99",
                                "#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a", "yellow"))

# Compute correlation for each system
cor_results0 <-   Predicted_cv0 %>%
  group_by(System) %>%
  summarise(R_value = cor(BLUP, Model3, use = "complete.obs"),
            RMSE = sqrt(mean((BLUP - Model3)^2, na.rm = TRUE)))

print(cor_results0)  # View correlation results

# Create scatter plot with regression lines and display R-values
ggplot(Predicted_cv0 , aes(x = Model3, y = BLUP, color = System)) +
  geom_point(size = 2, alpha = 0.7) +  # Scatter points
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Regression line
  facet_wrap(~ System) +  # Separate plots for each system
  theme_minimal() +
  labs(title = "Predicted vs. Actual BLUPs Across Systems",
       x = "Predicted Value",
       y = "Actual BLUP") +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 14, face = "bold"), # Style facet labels
        axis.title.x = element_text(size = 14),               # X-axis label size
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 14),         # Legend title size
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size = 12),  # X-axis values size
        axis.text.y = element_text(size = 12)) +  # Y-axis values size
  geom_text(data = cor_results0, 
            aes(x = min(Predicted_cv0$Model3, na.rm = TRUE), 
                y = max(Predicted_cv0$BLUP, na.rm = TRUE), 
                label = paste("R =", round(R_value, 2), "\nRMSE =", round(RMSE, 2))),
            hjust = 0, vjust = 1, size = 4, color = "black")



