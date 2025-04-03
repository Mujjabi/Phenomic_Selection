## Load Packages
library(dplyr)
library(tidyr)
library(lme4)
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

### Exploratory Analysis

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

##Estimate BLUES by System
model2 <- lm(YLD ~ HYB + Year + Treatments, data = COOP)
anova(model2)
emm_options(rg.limit = 200000)  # 
System_means <- emmeans(model2, ~ HYB | Treatments)
System_means <- as.data.frame(System_means)
write.csv(System_means, "System_means.csv")

#store Blues separately
Data_Organic <- System_means %>% filter(Treatments %in% c("Organic"))
Data_Conv    <- System_means %>% filter(Treatments %in% c("Conventional"))
                                 
                                                
## Estimate BLUPs for Management system and add a colomn to Blues
library(dplyr)
library(lme4)
Org_Loc_Blues  <- means %>% filter(Location %in% c("MRSE22", "MRSF23", "Urbana23", "MRSC24"))  # remove "Urbana24"
Conv_Loc_Blues <- means %>% filter(Location %in% c("BRKAb22", "CRL0622", "CRW1422", "KYS0222", "BRKC23", "BRKB24"))

## Fit BLUP Model to organic
mixed_Org <- lmer(emmean  ~ 1 + (1|HYB), weights = 1/Org_Loc_Blues$SE^2, data = Org_Loc_Blues)
blup_org = as.data.frame(coef(mixed_Org)$HYB)
blup_org <- data.frame(HYB = rownames(blup_org), BLUP = blup_org[,1], stringsAsFactors = FALSE)
Data_Organic <- merge(Data_Organic, blup_org, by = "HYB", all.x = TRUE)

## Fit BLUP Model to conventional
mixed_conv <- lmer(emmean  ~ 1 + (1|HYB), weights = 1/Conv_Loc_Blues$SE^2, data = Conv_Loc_Blues)
blup_conv = as.data.frame(coef(mixed_conv)$HYB)
blup_conv <- data.frame(HYB = rownames(blup_conv), BLUP = blup_conv[,1], stringsAsFactors = FALSE)
Data_Conv <- merge(Data_Conv, blup_conv, by = "HYB", all.x = TRUE)

#Combine both data for both systems
Combined_System_Blups <- rbind(Data_Organic, Data_Conv)
write.csv(Combined_System_Blups, "Combined_System_Blups.csv")

## Visualization of System BLUES
library(ggplot2)
ggplot(Combined_System_Blups, aes(x = emmean, color = Treatments, fill = Treatments)) +
  geom_histogram(alpha = 0.7, binwidth = 12, position = "identity") +
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

## Visualization of System BLUPS
library(ggplot2)
ggplot(Combined_System_Blups, aes(x = BLUP, color = Treatments, fill = Treatments)) +
  geom_histogram(alpha = 0.7, binwidth = 12, position = "identity") +
  labs(x = "Grain Yield BLUPs [bu/acre]", y = "Count", title = "Yield Distribution for 2022") +
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


## Visualization of Location BLUEs
library(ggplot2)
ggplot(Location_means, aes(x = Treatments, y = emmean, fill = Location)) +
  stat_boxplot(geom="errorbar") +
  geom_boxplot(outlier.colour="Black", outlier.shape=1,outlier.alpha = 0, #this removes outliers. 
               outlier.size=2, notch = FALSE,varwidth = FALSE) +
  labs(x = "Management Systems", y = "Grain Yield [Bu/Acre]", title = "Location BLUES") +
  theme_minimal() +
  #scale_fill_manual(values = c("BLUE" = "#1f78b4", "BLUP" = "#33a02c", "Predicted" = "#e31a1c"))

scale_color_manual(values = c("#a6cee3","#33a02c","#1f78b4","#b2df8a","#fb9a99",
                             "#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a", "yellow"))


library(ggplot2)
ggplot(Location_means, aes(x = Year, y = emmean, fill = Year)) +
  stat_boxplot(geom="errorbar") +
  geom_boxplot(outlier.colour="Black", outlier.shape=1,outlier.alpha = 0, #this removes outliers. 
               outlier.size=2, notch = FALSE,varwidth = FALSE) +
  labs(x = "Year", y = "Grain Yield [Bu/Acre]", title = "BLUEs by Year") +
  theme_minimal() +
  #scale_fill_manual(values = c("BLUE" = "#1f78b4", "BLUP" = "#33a02c", "Predicted" = "#e31a1c"))
  
  scale_color_manual(values = c("#a6cee3","#33a02c","#1f78b4","#b2df8a","#fb9a99",
                                "#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a", "yellow"))




## Load your SNP data 

### Convert genotypic data to numeric centered
### TASSEL converted SNPS into 0,0.5, 1 but I want 0, 1,1 

data <- read.csv("Markers_Numeric.csv", header = TRUE)


### Function to convert SNPs from 0, 0.5, 1 to -1, 0, 1. This helps for centering and scaling

convert_snps <- function(x) {
  # Mark the original 1s (e.g., by setting them to a temporary value like 999)
  original_zeros <- x == 0
  x[x == 0.5] <- 0  # Convert 0.5s to 1s
  x[original_zeros] <- -1  # Convert the original 1s to 2s
  return(x)
}

### Apply the conversion function to the data 
Geno_data <- apply(data[,-1], 2, convert_snps)


### Save the converted dataset
write.csv(Geno_data, "Geno_centered.csv", row.names = FALSE)

### Print a preview of the converted data
dim(Geno_data)


# Estimate kinship or Genomic Relationship Matrix (GRM) 
Geno_data  <- read.csv("Geno_centered.csv")
library(rrBLUP)
GRM  <- as.matrix(Geno_data)
GRM  <- apply(GRM, 2, as.numeric) 
GRM  <- A.mat(GRM )
rownames(GRM) <- colnames(GRM) <- Geno_data$Taxa
GRM[1:5, 1:5]
write.csv(GRM, "GRM.csv", row.names = TRUE, col.names = TRUE)

#Or use Gmatrix function
library(AGHmatrix)
Geno_matrix <- as.matrix(Geno_data)
Geno_matrix <- apply(Geno_matrix, 2, as.numeric)  # Convert all columns to numeric
Geno_matrix <- Geno_matrix + 1  #change it 0,1,2
GMatrix <- Gmatrix(Geno_matrix, method = "VanRaden")
GMatrix  <- as.data.frame(GMatrix )
rownames(GMatrix) <- colnames(GMatrix) <- Geno_data$Taxa
GMatrix [1:5, 1:5]

library(pheatmap)
# Assuming your kinship matrix is called 'kinship_matrix'

library(tibble)
Kinship_Tassel <- read.csv("Kinship_Matrix.csv", header = TRUE, row.names = 1)  # Load with first column as row names
Kinship_Tassel <- as.matrix(Kinship_Tassel)  # Convert to matrix


dim(Kinship_Tassel)
str(Kinship_Tassel)
sum(is.na(Kinship_Tassel))

# Define breaks to range from 0 to 2
breaks_seq <- seq(0,1, length.out = 100)

pheatmap(Kinship_Tassel, 
         cluster_rows = FALSE, #can enable clustering, but this rearrages genotypes
         cluster_cols = FALSE, 
         show_rownames = TRUE,  
         show_colnames = TRUE,  
         color = colorRampPalette(c("yellow", "white", "red"))(100),
         breaks = breaks_seq,  # Adjust scale range
         main = "Kinship Matrix Heatmap",
         fontsize_row = 1,  
         fontsize_col = 1,
         legend = TRUE,  # Show color legend
         legend_breaks = seq(0, 1, 0.5),  # Define key intervals
         legend_labels = seq(0, 1, 0.5),  
         display_numbers = FALSE)  # Remove numbers if needed for clarity



## Using GGplot2 for the same plot

library(ggplot2)
library(reshape2)
library(tibble)
Kinship_Tassel <- read.csv("Kinship_Matrix.csv", header = TRUE, row.names = 1)  # Load with first column as row names
Kinship_Tassel <- as.matrix(Kinship_Tassel)  # Convert to matrix

# Convert matrix to long format
kinship_long <- melt(Kinship_Tassel)
# Plot heatmap
ggplot(kinship_long, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "orange", midpoint = 0) +
  theme_minimal() +
  labs(title = "Kinship Matrix Heatmap", x = "Samples", y = "Samples") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 1),  # Adjust font size
    axis.text.y = element_text(size = 1)  )



# Calculate GCA and SCA effects
# Fit a mixed model to estimate GCA and SCA effects using lme4
COOP <- read.csv("COOP_Phenotypic.csv")
Location_means <- read.csv("All_Location_Blues.csv")
Modeling_HYB <- read.csv("Modeling_Set.csv")
GRM_All <- read.csv("Kinship_Matrix.csv")
rownames(GRM_All) <- GRM_All$Taxa  #make taxa rowname
GRM_All <- GRM_All[, -which(names(GRM_All) == "Taxa")]  # Remove 'Taxa' column


## Filter Pheno to select only hybrids included in modeling set
Filtered_Pheno <- Location_means %>%
  filter(HYB %in% Modeling_HYB$HYB)

Filtered_Pheno <- Filtered_Pheno %>%
  filter(complete.cases(.))

## Add parental IDs to Pheno

COOP_unique <- COOP %>%   #Deduplicate hybrids in original file
  distinct(HYB, .keep_all = TRUE) 

library(dplyr)
#Filtered_Pheno <- Filtered_Pheno %>%
                   #left_join(COOP_unique %>% select(HYB, Parent1, Tester), by = "HYB")  #it doesnt work sometimes

Filtered_Pheno <- merge(Filtered_Pheno, COOP_unique[, c("HYB", "Parent1", "Tester")], by = "HYB", all.x = TRUE)


### Divide the dataset into conventional and Organic to predict separately

library(dplyr)
Pheno_org   <- Filtered_Pheno %>% filter(Treatments == "Organic")
Pheno_conv  <- Filtered_Pheno %>% filter(Treatments == "Conventional")

## Organic
#we remove KTC029 because parent wasnt genotyped
Pheno_org <- Filtered_Pheno 
Pheno_org <- Pheno_org %>% filter(HYB != "KTC029")
write.csv(Pheno_org, "BLUEs_Organic.csv")


#Select the GRM only for inbreds in filtered
unique_parents <- unique(Pheno_org$Parent1)
length(unique_parents )
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



Pheno_org$Weights <- 1/(Pheno_org$SE^2)
# Fit the model using sommer
library(sommer)
model <- mmer(
  BLUE ~ 1 + System,  # Fixed effect (intercept only)
  random = ~ vsr(Parent1, Gu = GRM_Tested_O) + 
             vsr(Tester, Gu = GRM_Tester_O) + 
             vsr(Parent1_Tester, Gu = I_hybrids),
  weights = Weights,
  data = Pheno_org)

# Extract GCA and SCA effects
GCA_female <- model$U$`u:Parent1`$BLUE   # GCA for Parent1 (inbreds)
GCA_male   <- model$U$`u:Tester`$BLUE  # GCA for Tester (testers)
SCA        <- model$U$`u:Parent1_Tester`$BLUE  # SCA for Parent1:Tester (interaction)

# Match the levels in the prediction
predicted_performance <- model$Beta$Estimate[1] + 
                         GCA_female[Pheno_org$Parent1]+
                           GCA_male[Pheno_org$Tester] +
                                SCA[Pheno_org$Parent1_Tester]

# Add predicted performance to the dataset
Pheno_org$Predicted <- predicted_performance
PA_blue       <- cor(Pheno_org$BLUE, Pheno_org$Predicted, use = "complete.obs")
spearman_corr <- cor(Pheno_org$BLUE, Pheno_org$Predicted, use = "complete.obs", method = "spearman")


#Plot predicted, blue and blup
# Combine into one data frame
Pheno_org  <-  Pheno_org  %>%
  rename(BLUE = emmean)
library(tidyr)
plot_data  <- Pheno_org %>%
              pivot_longer(cols = c(BLUE, Predicted), 
              names_to = "Type", 
              values_to = "Value")

# Plot using ggplot
##Treant
# Create an index for x-axis
plot_data$Index <- seq_along(plot_data$Value)
library(ggplot2)
ggplot(plot_data, aes(x = Index, y = Value, color = Type)) +
  geom_point(size = 3) +                    # Adds points for each value
  geom_line(aes(group = Type), size = 1) +   # Adds lines connecting points for each type
  labs(x = "Index", y = "Value", title = "Comparison of Blues, BLUPs, and Predicted Values for Organic Locations") +
  theme_minimal() +
  scale_color_manual(values = c("BLUE" = "#1f78b4", "BLUP" = "#33a02c", "Predicted" = "#e31a1c"))

##Boxplot
library(ggplot2)
ggplot(plot_data, aes(x = System, y = Value, fill = Type)) +
  stat_boxplot(geom="errorbar") +
  geom_boxplot(outlier.colour="Black", outlier.shape=1,outlier.alpha = 0, #this removes outliers. 
               outlier.size=2, notch = FALSE,varwidth = FALSE) +
  labs(x = "Treatments", y = "Value", title = "Comparison of Blues, BLUPs, and Predicted Values") +
  theme_minimal() +
  scale_fill_manual(values = c("BLUE" = "#1f78b4", "BLUP" = "#33a02c", "Predicted" = "#e31a1c"))


### Predict GcA of new inbreds
# Fit the model using sommer

# Identify untested inbreds
GRM_All_Mat  <- as.matrix(GRM_All)
untested_inbreds <- setdiff(rownames(GRM_All_Mat), unique_parents)

# Predict GCA for untested inbreds using the relationship matrix
GCA_female_pred <- GRM_All_Mat[untested_inbreds, unique_parents] %*% model$U$`u:Parent1`$BLUE

# Store predictions in a dataframe
predicted_GCA_untested <- data.frame(Inbred = untested_inbreds, Predicted_GCA = GCA_female_pred)





### CV with 5 fold and 100 ITERACTIONS
# Step 1: Set up parameters
k <- 5  # Number of folds
n_iterations <- 100  # Number of iterations
set.seed(123)  # For reproducibility

# Initialize a data frame to store predictive ability values across all iterations
predictive_ability_all <- data.frame(
  Iteration = integer(),
  Fold = integer(),
  PA_Blue = numeric(),
  PA_Blup = numeric()
)

# Step 2: Perform 100 iterations of k-fold cross-validation
for (iter in 1:n_iterations) {
  # Assign each observation to a fold (randomly for each iteration)
  folds <- sample(1:k, nrow(Pheno_org), replace = TRUE)
  
  # Initialize a vector to store predictive ability values for this iteration
   PA_Blue <- numeric(k)
   PA_Blup <- numeric(k)
  
  # Step 3: Perform k-fold cross-validation for this iteration
  for (i in 1:k) {
    # Split the data into training and validation sets
    train_data <- Pheno_org[folds != i, ]
    valid_data <- Pheno_org[folds == i, ]
    
    # Fit the model on the training set
    model <- mmer(
      BLUE ~ 1,  # Fixed effect (intercept only)
      random = ~ vsr(Parent1, Gu = GRM_Tested_O) + 
                 vsr(Tester, Gu = GRM_Tester_O) + 
                 vsr(Parent1_Tester, Gu = I_hybrids),
      data = train_data
    )
    
    # Extract GCA and SCA effects
    GCA_female <- model$U$`u:Parent1`$BLUE  # GCA for Parent1 (inbreds)
    GCA_male   <- model$U$`u:Tester`$BLUE  # GCA for Tester (testers)
    SCA        <- model$U$`u:Parent1_Tester`$BLUE  # SCA for Parent1:Tester (interaction)
    
    # Predict performance for the validation set
    predicted <- model$Beta$Estimate + 
      GCA_female[valid_data$Parent1] +
      GCA_male[valid_data$Tester] +
      SCA[valid_data$Parent1_Tester]
    
    # Calculate predictive ability (correlation between observed and predicted values)
    PA_Blue[i] <- cor(valid_data$BLUE, predicted, use = "complete.obs")
    PA_Blup[i] <- cor(valid_data$BLUP, predicted, use = "complete.obs")
  }
  
  # Step 4: Store predictive ability values for this iteration
  predictive_ability_all <- rbind(
    predictive_ability_all,
    data.frame(
      Iteration = iter,
      Fold = 1:k,
      PA_Blue = PA_Blue,
      PA_Blup =  PA_Blup
      
    )
  )
}

# Step 5: Summarize and plot the results
# Calculate mean predictive ability across all iterations
PA_blup_mean <- mean(predictive_ability_all$PA_Blup)
PA_blue_mean <- mean(predictive_ability_all$PA_Blue)

# Plot the distribution of predictive ability values
PA_Plot   <- predictive_ability_all %>%
             pivot_longer(cols = c(PA_Blup, PA_Blue), 
               names_to = "Method", 
               values_to = "PA")

##Boxplot
library(ggplot2)
ggplot(PA_Plot, aes(x = PA, y = Method, fill = Method)) +
  stat_boxplot(geom="errorbar") +
  geom_boxplot(outlier.colour="Black", outlier.shape=1,outlier.alpha = 0, #this removes outliers. 
               outlier.size=2, notch = FALSE,varwidth = FALSE) +
  labs(y = "Method", x = "Predictive Ability", title = "Predictive Ability For Yield in Organic") +
  theme_minimal() +
  scale_fill_manual(values = c("PA_Blue" = "#1f78b4", "PA_Blup" = "#33a02c")) +
  theme(axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5))  # Remove x-axis labels




## Predictive ability for conventional fields 
Pheno_conv  <- Filtered_Pheno %>% filter(Treatments == "Conventional")
Pheno_org  <- Pheno_conv

## Organic
#we remove KTC029 because parent wasnt genotyped
Pheno_org <- Pheno_org %>% filter(HYB != "KTC029")

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

# Fit the model using sommer
library(sommer)
model <- mmer(
  emmean ~ 1,  # Fixed effect (intercept only)
  random = ~ vs(Parent1, Gu = GRM_Tested_O) + 
    vs(Tester, Gu = GRM_Tester_O) + 
    vs(Parent1_Tester, Gu = I_hybrids),
  data = Pheno_org)

#check if we have to remove rownames

# Extract GCA and SCA effects
GCA_female <- model$U$`u:Parent1`$emmean  # GCA for Parent1 (inbreds)
GCA_male   <- model$U$`u:Tester`$emmean # GCA for Tester (testers)
SCA        <- model$U$`u:Parent1_Tester`$emmean  # SCA for Parent1:Tester (interaction)

# Match the levels in the prediction
predicted_performance <- model$Beta$Estimate + 
  GCA_female[Pheno_org$Parent1]+
  GCA_male[Pheno_org$Tester] +
  SCA[Pheno_org$Parent1_Tester]

# Add predicted performance to the dataset
Pheno_org$Predicted <- predicted_performance
PA_blue <- cor(Pheno_org$emmean, Pheno_org$Predicted, use = "complete.obs")
PA_blup  <- cor(Pheno_org$BLUP, Pheno_org$Predicted, use = "complete.obs")


#Plot predicted, blue and blup
# Combine into one data frame
Pheno_org  <-  Pheno_org  %>%
  rename(BLUE = emmean)
plot_data  <- Pheno_org %>%
  pivot_longer(cols = c(BLUE, BLUP, Predicted), 
               names_to = "Type", 
               values_to = "Value")

# Plot using ggplot
##Treant
# Create an index for x-axis
plot_data$Index <- seq_along(plot_data$Value)
ggplot(plot_data, aes(x = Index, y = Value, color = Type)) +
  geom_point(size = 3) +                    # Adds points for each value
  geom_line(aes(group = Type), size = 1) +   # Adds lines connecting points for each type
  labs(x = "Index", y = "Value", title = "Comparison of Blues, BLUPs, and Predicted Values for Organic Locations") +
  theme_minimal() +
  scale_color_manual(values = c("BLUE" = "#1f78b4", "BLUP" = "#33a02c", "Predicted" = "#e31a1c"))

##Boxplot
library(ggplot2)
ggplot(plot_data, aes(x = Treatments, y = Value, fill = Type)) +
  stat_boxplot(geom="errorbar") +
  geom_boxplot(outlier.colour="Black", outlier.shape=1,outlier.alpha = 0, #this removes outliers. 
               outlier.size=2, notch = FALSE,varwidth = FALSE) +
  labs(x = "Treatments", y = "Value", title = "Comparison of Blues, BLUPs, and Predicted Values") +
  theme_minimal() +
  scale_fill_manual(values = c("BLUE" = "#1f78b4", "BLUP" = "#33a02c", "Predicted" = "#e31a1c"))


### Predict GcA of new inbreds
# Fit the model using sommer

# Identify untested inbreds
GRM_All_Mat  <- as.matrix(GRM_All)
untested_inbreds <- setdiff(rownames(GRM_All_Mat), unique_parents)

# Predict GCA for untested inbreds using the relationship matrix
GCA_female_pred <- GRM_All_Mat[untested_inbreds, unique_parents] %*% model$U$`u:Parent1`$emmean

# Store predictions in a dataframe
predicted_GCA_untested <- data.frame(Inbred = untested_inbreds, Predicted_GCA = GCA_female_pred)





### CV with 5 fold and 100 ITERACTIONS
# Step 1: Set up parameters
k <- 5  # Number of folds
n_iterations <- 100  # Number of iterations
set.seed(123)  # For reproducibility

# Initialize a data frame to store predictive ability values across all iterations
predictive_ability_all <- data.frame(
  Iteration = integer(),
  Fold = integer(),
  PA_Blue = numeric(),
  PA_Blup = numeric()
)

# Step 2: Perform 100 iterations of k-fold cross-validation
for (iter in 1:n_iterations) {
  # Assign each observation to a fold (randomly for each iteration)
  folds <- sample(1:k, nrow(Pheno_org), replace = TRUE)
  
  # Initialize a vector to store predictive ability values for this iteration
  PA_Blue <- numeric(k)
  PA_Blup <- numeric(k)
  
  # Step 3: Perform k-fold cross-validation for this iteration
  for (i in 1:k) {
    # Split the data into training and validation sets
    train_data <- Pheno_org[folds != i, ]
    valid_data <- Pheno_org[folds == i, ]
    
    # Fit the model on the training set
    model <- mmer(
      BLUE ~ 1,  # Fixed effect (intercept only)
      random = ~ vsr(Parent1, Gu = GRM_Tested_O) + 
        vsr(Tester, Gu = GRM_Tester_O) + 
        vsr(Parent1_Tester, Gu = I_hybrids),
      data = train_data
    )
    
    # Extract GCA and SCA effects
    GCA_female <- model$U$`u:Parent1`$BLUE  # GCA for Parent1 (inbreds)
    GCA_male   <- model$U$`u:Tester`$BLUE  # GCA for Tester (testers)
    SCA        <- model$U$`u:Parent1_Tester`$BLUE  # SCA for Parent1:Tester (interaction)
    
    # Predict performance for the validation set
    predicted <- model$Beta$Estimate + 
      GCA_female[valid_data$Parent1] +
      GCA_male[valid_data$Tester] +
      SCA[valid_data$Parent1_Tester]
    
    # Calculate predictive ability (correlation between observed and predicted values)
    PA_Blue[i] <- cor(valid_data$BLUE, predicted, use = "complete.obs")
    PA_Blup[i] <- cor(valid_data$BLUP, predicted, use = "complete.obs")
  }
  
  # Step 4: Store predictive ability values for this iteration
  predictive_ability_all <- rbind(
    predictive_ability_all,
    data.frame(
      Iteration = iter,
      Fold = 1:k,
      PA_Blue = PA_Blue,
      PA_Blup =  PA_Blup
      
    )
  )
}

# Step 5: Summarize and plot the results
# Calculate mean predictive ability across all iterations
PA_blup_mean <- mean(predictive_ability_all$PA_Blup)
PA_blue_mean <- mean(predictive_ability_all$PA_Blue)

# Plot the distribution of predictive ability values
PA_Plot   <- predictive_ability_all %>%
  pivot_longer(cols = c(PA_Blup, PA_Blue), 
               names_to = "Method", 
               values_to = "PA")

##Boxplot
library(ggplot2)
ggplot(PA_Plot, aes(x = PA, y = Method, fill = Method)) +
  stat_boxplot(geom="errorbar") +
  geom_boxplot(outlier.colour="Black", outlier.shape=1,outlier.alpha = 0, #this removes outliers. 
               outlier.size=2, notch = FALSE,varwidth = FALSE) +
  labs(y = "Method", x = "Predictive Ability", title = "Predictive Ability For Yield in Organic") +
  theme_minimal() +
  scale_fill_manual(values = c("PA_Blue" = "#1f78b4", "PA_Blup" = "#33a02c")) +
  theme(axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5))  # Remove x-axis labels

