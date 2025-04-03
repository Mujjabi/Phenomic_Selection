
                       ###########################
                       #  Asreml workflow        #
                       #     Chris Mujjabi       #
                       #          2025           #
                       ###########################
##Upload data
COOP <- read.csv("COOP_Phenotypic.csv")

## Set design variables as factors and data as numeric

COOP$HYB      <- as.factor(COOP$HYB)
COOP$Location <- as.factor(COOP$Location)
COOP$System <- as.factor(COOP$System)
COOP$Year <- as.factor(COOP$Year)
COOP$Rep <-  as.factor(COOP$Rep)
COOP$MST <- as.numeric(COOP$MST)
COOP$GWT <- as.numeric(COOP$GWT)
COOP$TWT <- as.numeric(COOP$TWT)
COOP$YLD <- as.numeric(COOP$YLD)


# Estimate location means

## Estimate BLUES by location and management system
library(emmeans)
# Fit the model
model <- lm(YLD ~ HYB + Year + Location + System, data = COOP)
anova(model)
# Estimate means for HYB within each Location
emm_options(rg.limit = 200000)  # 
means <- emmeans(model, ~ HYB | Location)
Location_means <- as.data.frame(means)
write.csv(means, "All_Location_Blues.csv")

## Estimate System Blups 
Location_means$Weights <- 1/(Location_means$SE^2) 


## Filter Pheno to select only hybrids included in modeling set
Pheno <- read.csv("System_means.csv")   #
Modeling_HYB <- read.csv("Modeling_Set.csv")  #hybrids with complete dataset
GRM_All <- read.csv("Kinship_Matrix.csv")
rownames(GRM_All) <- GRM_All$Taxa  #make taxa rowname
GRM_All <- GRM_All[, -which(names(GRM_All) == "Taxa")]  # Remove 'Taxa' column


Filtered_Pheno <- Pheno%>%
  filter(HYB %in% Modeling_HYB$HYB)

Filtered_Pheno <- Filtered_Pheno %>%
  filter(complete.cases(.))

#add parent ids
COOP_unique <- COOP %>%   #Deduplicate hybrids in original file
  distinct(HYB, .keep_all = TRUE) 
Filtered_Pheno <- merge(Filtered_Pheno, COOP_unique[, c("HYB", "Parent1", "Tester")], by = "HYB", all.x = TRUE)


#we remove KTC029 because parent wasnt genotyped
Filtered_Pheno  <- Filtered_Pheno  %>% filter(HYB != "KTC029")

Pheno_org <- Filtered_Pheno
length(unique(Pheno_org$Parent1))



#Select the GRM only for inbreds in filtered
unique_parents <- unique(Pheno_org$Parent1)
GRM_Tested_O <- GRM_All[unique_parents, unique_parents, drop = FALSE]
GRM_Tested_O[1:5, 1:5]
GRM_Tested_O  <- as.matrix(GRM_Tested_O)


## tester relationship matrix
Testers_o      <- unique(Pheno_org$Tester)
GRM_Tester_O   <- GRM_All[Testers_o, Testers_o, drop = FALSE]
GRM_Tester_O[1:5, 1:5]
GRM_Tester_O  <- as.matrix(GRM_Tester_O)


### Use the Sommer code
# Ensure Parent1_Tester interaction term is a factor with correct levels
Pheno_org$Parent1_Tester <- factor(interaction(Pheno_org$Parent1, Pheno_org$Tester, sep = ":"))

# Create an identity matrix for the interaction term with proper levels
n_hybrids <- length(unique(Pheno_org$Parent1_Tester))
I_hybrids <- diag(n_hybrids)
rownames(I_hybrids) <- levels(Pheno_org$Parent1_Tester)
colnames(I_hybrids) <- levels(Pheno_org$Parent1_Tester)


#Make columns as factors
Pheno_org <-Pheno_org %>% rename(BLUE = emmean)
Pheno_org <-Pheno_org %>% rename(System = Treatments)
Pheno_org$System <- as.factor(Pheno_org$System)
#Pheno_org$Location <- as.factor(Pheno_org$Location)
Pheno_org$Parent1 <- as.factor(Pheno_org$Parent1)
Pheno_org$Tester <- as.factor(Pheno_org$Tester)
Pheno_org$Parent1_Tester <- as.factor(Pheno_org$Parent1_Tester)

## Add weights toaccount for heterogeous erros

Pheno_org$Weights <- 1/(Pheno_org$SE^2)

library(asreml)
library(caret)  # For cross-validation folds

# Step 1: Create 5 folds for cross-validation
set.seed(123)
Pheno_org$Fold <- createFolds(Pheno_org$BLUE, k = 5, list = FALSE)

# Step 2: Initialize vectors for storing results
observed_all <- c()
predicted_all <- c()

# Step 3: Perform 5-fold cross-validation
for (fold in 1:5) {
  
  # Split data into training and validation
  train_data <- Pheno_org[Pheno_org$Fold != fold, ]
  test_data  <- Pheno_org[Pheno_org$Fold == fold, ]
  
  # Fit ASReml model on training data
  model_cv <- asreml(
    fixed = BLUE ~ Year + Location + System,  # Management system as a fixed effect
    random = ~ vm(Parent1, GRM_Tested_O) +    ##GCA for Parent1 (genomic effect)
               vm(Tester, GRM_Tester_O) +   ##GCA for Tester (genomic effect)
               vm(Parent1_Tester, I_hybrids) + #+ #SCA for hybrids (identity matrix)
               Location:Parent1_Tester, # Explicit interaction term,   
    weights = weight, #accounts for heterogeius environment
    data = train_data,
    na.action = na.method(x = "include", y = "include"),
    workspace = 6e8,  # Increase workspace if needed
    maxit = 5     # Maximum number of iterations
  )
predicted_values <- predict(model_cv, classify = "System:Parent1_Tester") #System:Parent1_Tester
predicted_pvals <- predicted_values$pvals
  
# Merge predicted values with test_data
test_data <- merge(test_data, predicted_pvals, by = c("System", "Parent1_Tester"), all.x = TRUE) #"System", "Parent1_Tester"
test_data$Predicted <- test_data$predicted.value
# Assign predicted values
test_data$Predicted <- test_data$predicted.value
  
# Store observed and predicted values
observed_all <- c(observed_all, test_data$BLUE)
predicted_all <- c(predicted_all, test_data$Predicted)
}

# Step 4: Evaluate Model Performance
cv_correlation <- cor(observed_all, predicted_all, use = "complete.obs")
cv_mse <- mean((observed_all - predicted_all)^2, na.rm = TRUE)


# Step 4: Evaluate Model Performance by System
##Conventional
Pred_conv  <- test_data %>% filter(System == "Conventional")
PA_conv <- cor(Pred_conv$BLUP,  Pred_conv$Predicted, use = "complete.obs")
mse_conv <- mean((Pred_conv$emmean -  Pred_conv$Predicted)^2, na.rm = TRUE)

##Organic
Pred_org  <- test_data %>% filter(System == "Organic")
PA_org <- cor(Pred_org$BLUP,  Pred_org$Predicted, use = "complete.obs")
mse_org <- mean((Pred_org$BLUP -  Pred_org$Predicted)^2, na.rm = TRUE)

# Calculate Spearman's rank correlation
spearman_corr <- cor(test_data$BLUE, test_data$Predicted, use = "complete.obs", method = "spearman")
# Calculate Kendall's rank correlation
kendall_corr <- cor(test_data$BLUE, test_data$Predicted, use = "complete.obs", method = "kendall")










# Step 5: Fit final ASReml model on all tested hybrids
final_model <- asreml(
  fixed = emmean ~ System,  # Management system as a fixed effect
  random = ~ vm(Parent1, GRM_Tested_O) +    ##GCA for Parent1 (genomic effect)
    vm(Tester, GRM_Tester_O) +   ##GCA for Tester (genomic effect)
    vm(Parent1_Tester, I_hybrids) + #SCA for hybrids (identity matrix)
    System:Parent1_Tester, # Explicit interaction term, 
                            
    data = Pheno_org,
  weights = Weights,  # Specify the column containing weights
  na.action = na.method(x = "include", y = "include"),
  workspace = 6e8,  # Increase workspace if needed
  maxit = 100     # Maximum number of iterations
)



# Step 6: Extract GCA estimates from the model
# Extract BLUPs (Best Linear Unbiased Predictions) for random effects
blups <- summary(final_model, coef = TRUE)$coef.random


# Extract GCA and SCA effects
GCA_female <- blups[grep("Parent1", rownames(blups)), "solution"]  # GCA for Parent1
GCA_male   <- blups[grep("Tester", rownames(blups)), "solution"]  # GCA for Tester
SCA        <- blups[grep("Parent1_Tester", rownames(blups)), "solution"]  # SCA for Parent1:Tester

# Predict performance for the validation set

fixed_effects <- summary(final_model, coef = TRUE)$coef.fixed

# Extract the intercept
intercept <- fixed_effects["(Intercept)", "solution"]

predicted_performance <- intercept  +
  GCA_female[Pheno_org$Parent1] +
  GCA_male[Pheno_org$Tester] +
  SCA[Pheno_org$Parent1_Tester]

Pheno_org$Predicted_asreml <- predicted_performance
predictive_ability <- cor(Pheno_org$emmean, Pheno_org$Predicted_asreml, use = "complete.obs")

# Calculate Spearman's rank correlation
spearman_corr <- cor(Pheno_org$emmean, Pheno_org$Predicted_asreml, use = "complete.obs", method = "spearman")
# Calculate Kendall's rank correlation
kendall_corr <- cor(Pheno_org$emmean, Pheno_org$Predicted_asreml, use = "complete.obs", method = "kendall")





# Step 7: Predict GCA for untested inbreds using GRM
untested_inbreds <- setdiff(rownames(GRM_All_Mat), unique(Pheno_org$Parent1))
GCA_female_pred <- GRM_All_Mat[untested_inbreds, unique(Pheno_org$Parent1)] %*% GCA_female  

# Step 8: Assign predicted GCA values to the untested dataframe
untested_crosses$GCA_Female <- GCA_female_pred[untested_crosses$Female]
untested_crosses$GCA_Male <- GCA_male[untested_crosses$Male]  

# Step 9: Estimate SCA (Assumption: Use mean SCA from tested hybrids)
mean_SCA <- mean(final_model$U$`u:Parent1_Tester`$emmean, na.rm = TRUE)
untested_crosses$SCA <- mean_SCA  # If specific SCA isn't available

# Step 10: Predict hybrid yield for untested crosses
mu <- final_model$Beta$Estimate[1]  # Intercept (overall mean)
untested_crosses$Predicted_Yield <- mu + untested_crosses$GCA_Female + untested_crosses$GCA_Male + untested_crosses$SCA

# Step 11: Print top predicted hybrids
untested_crosses <- untested_crosses[order(-untested_crosses$Predicted_Yield), ]
print(head(untested_crosses, 10))


# Step 4: Evaluate Model Performance by System
##Conventional
Pred_conv  <- test_data %>% filter(System == "Conventional")
PA_conv <- cor(Pred_conv$emmean,  Pred_conv$Predicted, use = "complete.obs")
mse_conv <- mean((Pred_conv$emmean -  Pred_conv$Predicted)^2, na.rm = TRUE)

##Organic
Pred_org  <- test_data %>% filter(System == "Organic")
PA_org <- cor(Pred_org$emmean,  Pred_org$Predicted, use = "complete.obs")
mse_org <- mean((Pred_org$emmean -  Pred_org$Predicted)^2, na.rm = TRUE)



# Initialize vectors for storing results
observed_all <- c()
predicted_all <- c()
fold_accuracies <- c()  # Vector to store accuracies for each fold

# Perform 5-fold cross-validation
for (fold in 1:5) {
  # Split data into training and validation sets
  train_data <- Pheno_org[Pheno_org$Fold != fold, ]
  test_data  <- Pheno_org[Pheno_org$Fold == fold, ]
  
  # Fit ASReml model on training data
  model_cv <- asreml(
    fixed = emmean ~ System,  # Management system as a fixed effect
    random = ~ vm(Parent1, GRM_Tested_O) +    # GCA for Parent1 (genomic effect)
      vm(Tester, GRM_Tester_O) +     # GCA for Tester (genomic effect)
      vm(Parent1_Tester, I_hybrids) + # SCA for hybrids (identity matrix)
      System:Parent1_Tester,  # Explicit interaction term
    data = train_data,
    weights = "Weights",  # Specify the column containing weights
    na.action = na.method(x = "include", y = "include"),
    workspace = 6e8,  # Increase workspace if needed
    maxit = 100       # Maximum number of iterations
  )
  
  # Predict on validation set
  predicted_values <- predict(model_cv, classify = "System:Parent1_Tester")
  predicted_pvals <- predicted_values$pvals
  
  # Merge predictions with test_data
  test_data <- merge(test_data, predicted_pvals, by = c("System", "Parent1_Tester"), all.x = TRUE)
  test_data$Predicted <- test_data$predicted.value
  
  # Store observed and predicted values
  observed_all <- c(observed_all, test_data$emmean)
  predicted_all <- c(predicted_all, test_data$Predicted)
  
  # Calculate accuracy (correlation) for the current fold
  fold_accuracy <- cor(test_data$emmean, test_data$Predicted, use = "complete.obs")
  
  # Store the accuracy for the current fold
  fold_accuracies <- c(fold_accuracies, fold_accuracy)
  
  # Print accuracy for the current fold
  cat("Fold", fold, "Accuracy (r):", round(fold_accuracy, 3), "\n")
}

# Step 5: Evaluate Model Performance
cv_correlation <- cor(observed_all, predicted_all, use = "complete.obs")
cv_mse <- mean((observed_all - predicted_all)^2, na.rm = TRUE)


# Print mean accuracy across all folds
mean_accuracy <- mean(fold_accuracies, na.rm = TRUE)
cat("Mean Accuracy (r) across 5 folds: ", round(mean_accuracy, 3), "\n")

boxplot(fold_accuracies)

fold_accuracies



