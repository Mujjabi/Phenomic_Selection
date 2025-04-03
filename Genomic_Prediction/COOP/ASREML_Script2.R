## Load Packages
library(dplyr)
library(tidyr)
library(lme4)
library(emmeans)
library(ggplot2)
library(car)
library(agricolae)


#Upload data
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
Pheno <- read.csv("Combined_System_Blups.csv")   #
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

# Create a system-specific identity matrix
I_hybrids_system <- diag(nrow(Pheno_org))
rownames(I_hybrids_system) <- paste(Pheno_org$System, Pheno_org$Parent1_Tester, sep = ":")
colnames(I_hybrids_system) <- paste(Pheno_org$System, Pheno_org$Parent1_Tester, sep = ":")


#Make columns as factors
Pheno_org$System <- as.factor(Pheno_org$System)
Pheno_org$Parent1 <- as.factor(Pheno_org$Parent1)
Pheno_org$Tester <- as.factor(Pheno_org$Tester)
Pheno_org$Parent1_Tester <- as.factor(Pheno_org$Parent1_Tester)

## Add weights toaccount for heterogeous erros

Pheno_org$Weights <- 1/(Pheno_org$SE^2)


## ASREML Model

# Load required libraries
library(dplyr)
library(asreml)
library(caret)

# Step 1: Prepare the data
# Assuming Pheno_org is your dataset with columns: Parent1_Tester, System, Location, BLUE, Weights, etc.
# Add a unique identifier for hybrid-system combinations
Pheno_org <- Pheno_org %>%
  mutate(hybrid_system = paste(Parent1_Tester, System, sep = "_"))


# Step 2: Set up CV1 cross-validation
set.seed(1343)  # For reproducibility
k <- 5  # Number of folds
n_iterations <- 1  # Number of iterations

# Initialize a data frame to store accuracy metrics for each fold in each iteration
accuracy_results <- data.frame(
  Iteration = integer(),
  Fold = integer(),
  Pearson = numeric(),
  Spearman = numeric()
)

# Step 3: Perform 100 iterations of k-fold CV
for (iter in 1:n_iterations) {
  # Create folds for this iteration
  hybrid_systems <- unique(Pheno_org$hybrid_system)
  library(splitTools)
  #folds <- createFolds(hybrid_systems, k = k, list = TRUE, returnTrain = FALSE)
  folds <- split(sample(nrow(Pheno_org)), rep(1:k, length.out = nrow(Pheno_org)))
  #folds <- create_folds(Pheno_org$hybrid_system, k = k)
  
  # Perform k-fold CV for this iteration
  for (i in 1:length(folds)) {
    # Get hybrid-system combinations for the validation set
    validation_hybrid_systems <- hybrid_systems[folds[[i]]]
    
    # Split data into training and validation sets
    validation_data <- Pheno_org %>% filter(hybrid_system %in% validation_hybrid_systems)
    training_data   <- Pheno_org %>% filter(!hybrid_system %in% validation_hybrid_systems)
    
    # Fit the model
    model_cv <- asreml(
      fixed = BLUE ~ 1,  # Fixed effects
      random = ~ vm(Parent1, GRM_Tested_O) +  # GCA for Parent1
        vm(Tester, GRM_Tester_O) +   # GCA for Tester
        vm(Parent1_Tester, I_hybrids) +  # SCA for hybrids
        vm(System:Parent1_Tester, Gu = I_hybrids_system) +
        System:Parent1_Tester,  # Interaction between System and Hybrid
      weights = Weights,  # Weights for heterogeneous variances
      data = training_data,
      na.action = na.method(x = "include", y = "include"),
      workspace = 6e8,  # Increase workspace if needed
      maxit = 100  # Increase the number of iterations
    )
 
    # Predict on the validation set
    predicted_values <- predict(model_cv, newdata = validation_data, classify = "System:Parent1_Tester")
    predicted_pvals  <- predicted_values$pvals
    
    # Merge predicted values with validation_data
    validation_data <- merge(validation_data, predicted_pvals, by = c("Parent1_Tester", "System"), all.x = FALSE)
    
    # Assign predicted values to a new column
    validation_data$Predicted <- validation_data$predicted.value
    
    # Calculate accuracy metrics
    pearson_corr <- cor(validation_data$BLUE, validation_data$Predicted, use = "complete.obs")
    spearman_corr <- cor(validation_data$BLUE, validation_data$Predicted, use = "complete.obs", method = "spearman")
    
    # Store accuracy metrics for this fold and iteration
    accuracy_results <- rbind(accuracy_results, data.frame(
      Iteration = iter,
      Fold = i,
      Pearson = pearson_corr,
      Spearman = spearman_corr
    ))
  }
  
}


# Step 5: Plot Accuracies
library(tidyr)
library(ggplot2)

# Reshape accuracy_results for plotting
accuracy_results_long <- accuracy_results %>%
  pivot_longer(cols = c(Pearson, Spearman), 
               names_to = "Method", 
               values_to = "Accuracy")

# Plot boxplots of accuracy metrics
ggplot(accuracy_results_long, aes(x = Method, y = Accuracy, fill = Method)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot(outlier.colour = "Black", outlier.shape = 1, outlier.alpha = 0, # Remove outliers
               outlier.size = 2, notch = FALSE, varwidth = FALSE) +
  labs(x = "Correlation Method", y = "Accuracy", title = "Predictive Ability for Grain Yield") +
  theme_minimal() +
  scale_fill_manual(values = c("Pearson" = "#1f78b4", "Spearman" = "#e31a1c"))
