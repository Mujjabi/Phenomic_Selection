## Genomic Prediction of Single Crosses

## First estimate the blues of each study independently. Get Blues and errors associated to each study
 
# Load the data
COOP <- read.csv("ALL_COOP.csv")


COOP$HYB      <- as.factor(COOP$HYB)
COOP$Location <- as.factor(COOP$Location)
COOP$Treatments <- as.factor(COOP$Treatments)
COOP$Year <- as.factor(COOP$Year)
COOP$REP <-  as.factor(COOP$Rep)
COOP$MST <- as.numeric(COOP$MST)
COOP$GWT <- as.numeric(COOP$GWT)
COOP$TWT <- as.numeric(COOP$TWT)
COOP$YLD <- as.numeric(COOP$YLD)


COOP2 <- COOP %>%    #exclude the urbana location
  filter(!(Year == 2024 & Location == "Urbana24"))


# Load necessary libraries
# Load necessary libraries
library(lme4)
library(emmeans)
library(dplyr)

# Example data structure
# Assuming your data is in a data frame called 'data'
# Columns: Year, Location, Rep, Hybrid, Yield

# Fit a mixed model for each year separately
# Hybrid is fixed, Location and Rep are random
years <- unique(COOP2$Year)  # Get unique years
blue_results_list <- list() # Store BLUEs for each year

for (year in years) {
  # Subset data for the current year
  year_data <- subset(COOP2[, c(2, 6, 9,12, 26)], Year == year)
  
  # Fit the mixed model
  model <- lmer(YLD ~ HYB + (1 | Location/Rep), data = year_data)
  
  # Estimate BLUEs for hybrids using emmeans
  blue <- emmeans(model, ~ HYB)
  
  # Convert emmeans object to a dataframe
  blue_df <- as.data.frame(blue)
  blue_df$Year <- year  # Add year column
  
  # Store the results in the list
  blue_results_list[[as.character(year)]] <- blue_df
}

# Combine all BLUEs into a single dataframe
blue_results <- bind_rows(blue_results_list)

write.csv(blue_results, "Hybrid_Year_BLUES.csv")


# Filter data for 2024
BLUES_2022 <- blue_results %>% filter(Year == 2022)
BLUES_2023 <- blue_results %>% filter(Year == 2023)
BLUES_2024 <- blue_results %>% filter(Year == 2024)


# Extract specific hybrids for 2024
highlight_2022 <- BLUES_2022 %>%
  filter(HYB %in% c("UIUC502", "UIUC524", "VIKING51", "BR59H44"))

highlight_2023 <- BLUES_2023 %>%
  filter(HYB %in% c("UIUC502", "UIUC524", "VIKING51", "BR59H44"))

highlight_2024 <- BLUES_2024 %>%
  filter(HYB %in% c("UIUC502", "UIUC524", "ISU174","ISU176"))


# Plot for 2022
ggplot(BLUES_2022, aes(x = emmean)) +
  geom_histogram(alpha = 0.7, binwidth = 12, fill = "darkgoldenrod") +
  geom_vline(data = highlight_2022, aes(xintercept = emmean, color = HYB),
             linetype = "dashed", size = 1) +
  labs(x = "Grain Yield BLUES [bu/acre]", y = "Count", title = "Yield Distribution for 2022") +
  theme_minimal() +
  theme(
    legend.position = "right",  # Place legend at the bottom
    strip.text = element_text(size = 14, face = "bold"),  # Increase facet label font size
    axis.title.x = element_text(size = 14, face = "bold"),  # Increase x-axis title font size
    axis.title.y = element_text(size = 14, face = "bold"),  # Increase y-axis title font size
    axis.text.x = element_text(size = 12),  # Increase x-axis tick label font size
    axis.text.y = element_text(size = 12),
    panel.grid.minor = element_blank()) + 
  theme(legend.position = "right") +  # Remove legend
  scale_x_continuous(limits = c(80, 230), breaks = seq(80, 230, by = 20)) 



# Plot for 2023
ggplot(BLUES_2023, aes(x = emmean)) +
  geom_histogram(alpha = 0.5, binwidth = 12, fill = "#33A02C") +
  geom_vline(data = highlight_2023, aes(xintercept = emmean, color = HYB),
             linetype = "dashed", size = 1) +
  labs(x = "Grain Yield BLUES [bu/acre]", y = "Count", title = "Yield Distribution for 2022") +
  theme_minimal() +
  theme(
    legend.position = "bottom",  # Place legend at the bottom
    strip.text = element_text(size = 14, face = "bold"),  # Increase facet label font size
    axis.title.x = element_text(size = 14, face = "bold"),  # Increase x-axis title font size
    axis.title.y = element_text(size = 14, face = "bold"),  # Increase y-axis title font size
    axis.text.x = element_text(size = 12),  # Increase x-axis tick label font size
    axis.text.y = element_text(size = 12),
    panel.grid.minor = element_blank()) + 
  theme(legend.position = "right") +  # Remove legend
  scale_x_continuous(limits = c(60, 220), breaks = seq(60, 220, by = 20)) 

# Plot for 2024

ggplot(BLUES_2024, aes(x = emmean)) +
  geom_vline(data = highlight_2024, aes(xintercept = emmean, color = HYB),
             linetype = "dashed", size = 1, alpha = 0.7) +  
  geom_histogram(alpha = 1, bins = 15, fill = "darkgray", color = "black") +  
  
  scale_color_manual(values = c("#1b9e77","#e41a1c","#377eb8","#e7298a")) +
  labs(x = "Grain Yield BLUES [bu/acre]", y = "Count", title = "Yield Distribution Across 2024 Locations") +
  theme_minimal() +
  theme(
    legend.position = "right",  # Remove legend
    strip.text = element_text(size = 14, face = "bold"),  # Increase facet label font size
    axis.title.x = element_text(size = 14, face = "bold"),  # Increase x-axis title font size
    axis.title.y = element_text(size = 14, face = "bold"),  # Increase y-axis title font size
    axis.text.x = element_text(size = 12),  # Increase x-axis tick label font size
    axis.text.y = element_text(size = 12),
    panel.grid.minor = element_blank()
  ) +  
  scale_x_continuous(limits = c(50, 220), breaks = seq(50, 220, by = 50))



# Define the checks. These will be removed from predictions since we dont have parental information
Checks <- c("PH6878",
            "BR59H44", 
            "BR70A47", 
            "BR65G55",  
            "VIKING51", 
            "VIKING72")


## Genomic Selection 


# Load necessary libraries
library(lme4)
library(rrBLUP)
library(Matrix)

# Load your datasets
# Phenotypic data: single-cross performance

pheno_data <- read.csv("pheno_data.csv")  # These would be the BLUES 

# Genotypic data: SNP markers for parents
geno_data <- read.csv("geno_data.csv")   #This is the snp data

# Genomic Relationship Matrix (GRM) 
GRM <- read.csv("GRM.csv")  

# Prepare the data
# Assuming pheno_data has columns: CrossID, Female, Male, TraitValue
# geno_data has columns: ParentID, SNP1, SNP2, ..., SNPN
# GRM is a matrix with ParentIDs as row and column names

# Step 1: Calculate GCA and SCA effects
# Fit a mixed model to estimate GCA and SCA effects using lme4
pheno_data$Female <- as.factor(pheno_data$Female)
pheno_data$Male <- as.factor(pheno_data$Male)

# Fit the model using lme4

#Weights: The reciprical of the error variance.
#Values that have a lot of error will be contributing less and values will less error will contribute more

model <- lmer(TraitValue ~ 1 + (1|Female) + (1|Male) + (1|Female:Male), data = pheno_data)
model <- lmer(TraitValue ~ 1 + (1|Female) + (1|Male) + (1|Female:Male), weights = 1/BLUES$SE^2, data = pheno_data)  #this includes weights

# Extract GCA and SCA effects
summary_model <- summary(model)
GCA_female <- summary_model$varcor$Female[1]
GCA_male <- summary_model$varcor$Male[1]
SCA <- summary_model$varcor$`Female:Male`[1]

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