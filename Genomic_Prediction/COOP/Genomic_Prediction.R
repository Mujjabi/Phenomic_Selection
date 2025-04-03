## Genomic Prediction of Single Crosses

# Load necessary libraries
library(lme4)
library(rrBLUP)
library(Matrix)
library(dplyr)
library(tidyr)
library(emmeans)
library(ggplot2)
library(car)
library(agricolae)


## Load all our phenotypic Data

COOP <- read.csv("COOP_Phenotypic.csv")

## Set design variables as factors and data as numeric

COOP$HYB      <- as.factor(COOP$HYB)
COOP$Location <- as.factor(COOP$Location)
COOP$Treatments <- as.factor(COOP$Treatments)
COOP$Year <- as.factor(COOP$Year)
COOP$Rep <-  as.factor(COOP$Rep)
COOP$MST <- as.numeric(COOP$MST)
COOP$GWT <- as.numeric(COOP$GWT)
COOP$TWT <- as.numeric(COOP$TWT)
COOP$YLD <- as.numeric(COOP$YLD)

# Load your datasets
# Phenotypic data: single-cross performance

# Calculate GCA and SCA effects
# Fit a mixed model to estimate GCA and SCA effects using lme4

Pheno <- read.csv("System_BlUEs.csv")
Modeling_HYB <- read.csv("Modeling_Set.csv")
GRM_All <- read.csv("Kinship_Matrix.csv")


## Filter Pheno to select only hybrids included in modeling set

Filtered_Pheno <- Pheno %>%
  filter(HYB %in% Modeling_HYB$HYB)

Filtered_Pheno <- Filtered_Pheno %>%
  filter(complete.cases(.))

## Add parental IDs to Pheno

COOP_unique <- COOP %>%   #Deduplicate hybrids in original file
  distinct(HYB, .keep_all = TRUE) 

Filtered_Pheno <- Filtered_Pheno %>%
  left_join(COOP_unique %>% select(HYB, Parent1, Tester), by = "HYB")


#Select the GRM only for inbreds in filtered
tested_inbreds <- intersect(unique(Filtered_Pheno$Parent1), GRM_All$Taxa)
tested_inbreds <- tested_inbreds[tested_inbreds %in% GRM_All$Taxa]
# Step 2: Subset the GRM to include only tested inbreds

GRM_Tested <- GRM_All[match(tested_inbreds, GRM_All$Taxa), 
                      c(TRUE, match(tested_inbreds, GRM_All$Taxa) + 1)]

GRM_Tested[1:5, 1:5]

## tester relationship matrix
Testers <- intersect(unique(Filtered_Pheno$Tester), GRM_All$Taxa)
Testers <- Testers[Testers %in% GRM_All$Taxa]

# Step 2: Subset the GRM to include only tested inbreds

GRM_Tester <- GRM_All[match(Testers, GRM_All$Taxa), 
                      c(TRUE, match(Testers, GRM_All$Taxa) + 1)]
GRM_Tester[1:5, 1:5]

#this repeats inbreds in the same order as they appear in the hybrid list
#GRM_Tested <- GRM_All[match(tested_inbreds, GRM_All$Taxa), match(tested_inbreds, GRM_All$Taxa)]

#set parents as factors
Filtered_Pheno$Parent1 <- as.factor(Filtered_Pheno$Parent1)
Filtered_Pheno$Tester   <- as.factor(Filtered_Pheno$Tester)




# Prepare the data
# Assuming pheno_data has columns: CrossID, Female, Male, TraitValue
# geno_data has columns: ParentID, SNP1, SNP2, ..., SNPN
# GRM is a matrix with ParentIDs as row and column names

# Step 1: Calculate GCA and SCA effects
# Fit a mixed model to estimate GCA and SCA effects using lme4
pheno_data$Female <- as.factor(pheno_data$Female)
pheno_data$Male <- as.factor(pheno_data$Male)


# Fit the model using lme4
model <- lmer(emmean ~ 1 + (1|System) + (1|Parent1) + (1|Tester) + (1|Parent1:Tester), weights = 1/Filtered_Pheno$SE^2, data = Filtered_Pheno)

# Extract GCA and SCA effects
summary_model <- summary(model)
GCA_female <- summary_model$varcor$Parent1[1]
GCA_male <- summary_model$varcor$Tester[1]
SCA <- summary_model$varcor$`Parent1:Tester`[1]
GCA_female2 <-  coef(model)$Parent1[1]

random_effects <- ranef(model)
GCA_female2  <- random_effects$Parent1
GCA_male2  <- random_effects$Tester
SCA2  <- random_effects$`Parent1:Tester`


# Step 2: Genomic Prediction
# Use the rrBLUP package to predict GCA and SCA effects based on genomic data
# Assuming geno_data is in the format: ParentID, SNP1, SNP2, ..., SNPN
# Convert geno_data to a matrix
geno_matrix <- as.matrix(geno_data[, -1])  # Remove the ParentID column
rownames(geno_matrix) <- geno_data$ParentID

# Calculate the genomic relationship matrix (if not already provided)
# GRM <- A.mat(geno_matrix)  # Uncomment if you need to calculate GRM

# Predict GCA effects using rrBLUP
GCA_female_pred <- mixed.solve(y = GCA_female, Z = GRM)
GCA_male_pred <- mixed.solve(y = GCA_male, Z = GRM)

# Predict SCA effects using rrBLUP
SCA_pred <- mixed.solve(y = SCA, Z = GRM)

# Step 3: Predict Single-Cross Performance
# Combine GCA and SCA predictions to predict single-cross performance
predict_single_cross <- function(female, male, GCA_female_pred, GCA_male_pred, SCA_pred) {
  GCA_female <- GCA_female_pred$u[female]
  GCA_male <- GCA_male_pred$u[male]
  SCA <- SCA_pred$u[paste(female, male, sep = ":")]
  return(GCA_female + GCA_male + SCA)
}

# Example: Predict performance for a new single cross
new_female <- "Parent1"  # Replace with actual female parent ID
new_male <- "Parent2"    # Replace with actual male parent ID
predicted_performance <- predict_single_cross(new_female, new_male, GCA_female_pred, GCA_male_pred, SCA_pred)
print(predicted_performance)

# Step 4: Cross-Validation (Optional)
# Implement cross-validation to estimate prediction accuracy
# This is a simplified example; you may need to adapt it to your specific dataset
n <- nrow(pheno_data)
accuracy <- numeric(n)

for (i in 1:n) {
  # Leave-one-out cross-validation
  train_data <- pheno_data[-i, ]
  test_data <- pheno_data[i, ]
  
  # Fit the model on the training data
  model_cv <- lmer(TraitValue ~ 1 + (1|Female) + (1|Male) + (1|Female:Male), data = train_data)
  
  # Predict the performance of the test cross
  predicted_value <- predict_single_cross(test_data$Female, test_data$Male, GCA_female_pred, GCA_male_pred, SCA_pred)
  
  # Calculate accuracy (correlation between observed and predicted values)
  accuracy[i] <- cor(test_data$TraitValue, predicted_value)
}

mean_accuracy <- mean(accuracy, na.rm = TRUE)
print(paste("Mean prediction accuracy:", mean_accuracy))