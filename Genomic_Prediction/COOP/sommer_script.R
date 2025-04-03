

###########################
#  Asreml workflow        #
#     Chris Mujjabi       #
#          2025           #
###########################
# Load Packages
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
COOP$System<- as.factor(COOP$System)
COOP$Year <- as.factor(COOP$Year)
COOP$Rep <-  as.factor(COOP$Rep)
COOP$MST <- as.numeric(COOP$MST)
COOP$GWT <- as.numeric(COOP$GWT)
COOP$TWT <- as.numeric(COOP$TWT)
COOP$YLD <- as.numeric(COOP$YLD)

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


## Visualization of Location BLUEs
Location_blues <- read.csv("All_Location_Blues.csv")
library(ggplot2)
ggplot(Location_blues, aes(x = System, y = BLUE, fill = Location)) +
  stat_boxplot(geom="errorbar") +
  geom_boxplot(outlier.colour="Black", outlier.shape=1,outlier.alpha = 0, #this removes outliers. 
               outlier.size=2, notch = FALSE,varwidth = FALSE) +
  labs(x = "Management Systems", y = "Grain Yield [Bu/Acre]", title = "Adjusted Location Means") +
  theme_minimal() +
  theme(
    legend.position = "right",  # Place legend at the bottom
    strip.text = element_text(size = 14, face = "bold"),  # Increase facet label font size
    axis.title.x = element_text(size = 14, face = "bold"),  # Increase x-axis title font size
    axis.title.y = element_text(size = 14, face = "bold"),  # Increase y-axis title font size
    axis.text.x = element_text(size = 12),  # Increase x-axis tick label font size
    axis.text.y = element_text(size = 12),
    panel.grid.minor = element_blank()) + 
  theme(legend.position = "right") +  # Rem
  scale_y_continuous(limits = c(0, 250), breaks = seq(0, 250, by = 50)) +
  
  scale_color_manual(values = c("#a6cee3","#33a02c","#1f78b4","#b2df8a","#fb9a99",
                                "#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a", "yellow"))

## BLUES
library(ggplot2)
Location_blues$Year <- as.factor(Location_blues$Year)
ggplot(Location_blues, aes(x = Year, y = BLUE, fill = Year)) +
  stat_boxplot(geom="errorbar") +
  geom_boxplot(outlier.colour="Black", outlier.shape=1,outlier.alpha = 0, #this removes outliers. 
               outlier.size=2, notch = FALSE,varwidth = FALSE) +
  labs(x = "Year", y = "Grain Yield [Bu/Acre]", title = "Adjusted Means Across Years") +
  theme_minimal() +
  theme(
    legend.position = "right",  # Place legend at the bottom
    strip.text = element_text(size = 14, face = "bold"),  # Increase facet label font size
    axis.title.x = element_text(size = 14, face = "bold"),  # Increase x-axis title font size
    axis.title.y = element_text(size = 14, face = "bold"),  # Increase y-axis title font size
    axis.text.x = element_text(size = 12),  # Increase x-axis tick label font size
    axis.text.y = element_text(size = 12),
    panel.grid.minor = element_blank()) + 
  theme(legend.position = "right") +  # Rem
  scale_y_continuous(limits = c(0, 250), breaks = seq(0, 250, by = 50)) +
  scale_color_manual(values = c("#a6cee3","#33a02c","#1f78b4","#b2df8a","#fb9a99",
                                "#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a", "yellow"))


## Visualization of System BLUPS
Combined_System_Blups <- read.csv("Combined_System_Blups.csv")
library(ggplot2)
ggplot(Combined_System_Blups, aes(x = BLUP, color = System, fill = System)) +
  geom_histogram(alpha = 0.7, binwidth = 12, position = "identity") +
  labs(x = "Grain Yield BLUPs [bu/acre]", y = "Count", title = "Grain Yield BLUP Distribution") +
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

## System Boxplot
library(ggplot2)
ggplot(Combined_System_Blups, aes(x = System, y = BLUP, fill = System)) +
  stat_boxplot(geom="errorbar") +
  geom_boxplot(outlier.colour="Black", outlier.shape=1,outlier.alpha = 0, #this removes outliers. 
               outlier.size=2, notch = FALSE,varwidth = FALSE) +
  labs(x = "Management Systems", y = "Grain Yield [Bu/Acre]", title = "Grain Yield BLUPs By System") +
  theme_minimal() +
  theme(
    legend.position = "right",  # Place legend at the bottom
    strip.text = element_text(size = 14, face = "bold"),  # Increase facet label font size
    axis.title.x = element_text(size = 14, face = "bold"),  # Increase x-axis title font size
    axis.title.y = element_text(size = 14, face = "bold"),  # Increase y-axis title font size
    axis.text.x = element_text(size = 12),  # Increase x-axis tick label font size
    axis.text.y = element_text(size = 12),
    panel.grid.minor = element_blank()) + 
  theme(legend.position = "right") +  # Rem
  scale_y_continuous(limits = c(80, 220), breaks = seq(80, 220, by = 40)) +
  scale_color_manual(values = c("#a6cee3","#33a02c","#1f78b4","#b2df8a","#fb9a99",
                                "#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a", "yellow"))





                                  ####################################
                                  ##       Marker Data Wrangling    ##
                                  ##                                ##
                                  ####################################
## Load your SNP data 
data <- read.csv("Markers_Numeric.csv", header = TRUE)

### Convert genotypic data to numeric centered, TASSEL converted SNPS into 0,0.5, 1 but I want to center it to -1, 0, 1
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


# Estimate kinship or Genomic Relationship Matrix (GRM) 
Geno_data  <- read.csv("Geno_centered.csv")
library(rrBLUP)
GRM  <- A.mat(Geno_data)
GRM[1:5, 1:5]
rownames(GRM) <- colnames(GRM) <- rownames(Geno_data)
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

install.packages("ASRgenomics")
library(ASRgenomics)
GMatrix <- as.matrix(GMatrix)
Ghat.blend <- G.tuneup(G=GMatrix, blend=TRUE, pblend=0.02)$Gb
Ghat.blend[1:8,1:8]



# Assuming your kinship matrix is called 'kinship_matrix'
library(pheatmap)
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



               ####################################
               ##       Genomic Prediction       ##
               ##                                ##
               ####################################


# Fit a mixed model to estimate GCA and SCA effects using lme4
library(dplyr)
COOP <- read.csv("COOP_Phenotypic.csv")
System_Blues <- read.csv("Combined_System_Blups.csv")  #has both Blues and Blups for each system
Modeling_HYB <- read.csv("Modeling_Set.csv")
New_crosses <- read.csv("Validation_Hybrids.csv")

GRM_All <- read.csv("Kinship_Matrix.csv")
rownames(GRM_All) <- GRM_All$Taxa  #make taxa rowname
GRM_All <- GRM_All[, -which(names(GRM_All) == "Taxa")]  # Remove 'Taxa' column


## Filter Pheno to select only hybrids included in modeling set
Filtered_Pheno <- System_Blues %>%
  filter(HYB %in% Modeling_HYB$HYB)

Filtered_Pheno <- Filtered_Pheno %>%
  filter(complete.cases(.))

## Add parental IDs to Pheno
COOP_unique <- COOP %>%   #Deduplicate hybrids in original file
  distinct(HYB, .keep_all = TRUE) 

#merge columns
Filtered_Pheno <- merge(Filtered_Pheno, COOP_unique[, c("HYB", "Parent1", "Tester")], by = "HYB", all.x = TRUE)


## Organic
#we remove KTC029 because parent wasnt genotyped
Pheno_org <- Filtered_Pheno 
Pheno_org <- Pheno_org %>% filter(HYB != "KTC029") #no parent genotic data found


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
Pheno_org$Weights <- 1/(Pheno_org$SE^2)  #Get weights to account for heterogeonous environments


# Compute the inverse of the relationship matrices
#Blend
install.packages("ASRgenomics")
library(ASRgenomics)

Ghat.blend <- G.tuneup(G=GRM_Tested_O, blend=TRUE, pblend=0.02)$Gb
Ghat.blend[1:8,1:8]

Ghat.blend_tester <- G.tuneup(G=GRM_Tester_O, blend=TRUE, pblend=0.02)$Gb
Ghat.blend_tester[1:8,1:8]

Ghat.blend_Ihybrid <- G.tuneup(G=I_hybrids, blend=TRUE, pblend=0.02)$Gb
Ghat.blend_Ihybrid[1:8,1:8]


#get inverse
Ginv.blend <- G.inverse(G=Ghat.blend)$Ginv
Ginv.blend[1:8,1:8]

Ginv.blend_tester <- G.inverse(G=Ghat.blend_tester)$Ginv

Ginv.blend_Ihybrid <- G.inverse(G=Ghat.blend_Ihybrid)$Ginv


# 5-fold cross-validation:
# Load required libraries
library(dplyr)
library(sommer)
library(splitTools)


# Set seed for reproducibility
set.seed(1343)  
#Number of folds and iterations
num_folds <- 5 
n_iterations <- 2  

# Step 1: Identify hybrids with both systems
hybrids_with_both_systems <- Pheno_org %>%
  group_by(Parent1_Tester) %>%
  summarise(Systems = n_distinct(System)) %>%
  filter(Systems == 2) %>%
  pull(Parent1_Tester)

# Step 2: Assign Hybrids to Folds (Ensuring No Hybrid Appears in Both Systems in the Same Fold)
fold_assignments <- data.frame(Parent1_Tester = unique(Pheno_org$Parent1_Tester))
fold_assignments$Fold <- sample(rep(1:num_folds, length.out = nrow(fold_assignments)))

# Assign folds to the dataset
Pheno_org <- Pheno_org %>%
  left_join(fold_assignments, by = "Parent1_Tester")

# Step 3: Ensure hybrids with both systems do not appear in both systems within the same fold
test_data_list <- list()
train_data_list <- list()

for (i in 1:num_folds) {
  # Select test hybrids for the current fold
  test_hybrids <- fold_assignments %>%
    filter(Fold == i) %>%
    pull(Parent1_Tester)
  
  # Select test data ensuring only one system per hybrid
  test_data_list[[i]] <- Pheno_org %>%
    filter(Parent1_Tester %in% test_hybrids) %>%
    group_by(Parent1_Tester) %>%
    filter(System == sample(System, 1))  # Randomly select one system per hybrid
  
  # Remove test_data from Pheno_org to create train_data for this fold
  train_data_list[[i]] <- Pheno_org %>%
    anti_join(test_data_list[[i]], by = c("Parent1_Tester", "System"))
}


# Initialize accuracy results storage
accuracy_results <- data.frame(
  Iteration = integer(),
  Fold = integer(),
  PSM1 = numeric(),
  PSM2 = numeric(),
  PSM3 = numeric(),
  SPM1 = numeric(),
  SPM2 = numeric(),
  SPM3 = numeric()
)

# Step 4: Perform K-Fold Cross-Validation with `sommer`
for (iter in 1:n_iterations) {
  for (i in 1:num_folds) {
    # Extract training and test data for the current fold
    training_data <- train_data_list[[i]]
    validation_data <- test_data_list[[i]]
    
    # Fit the sommer model on the training data
    
    model <- mmer(
      BLUE ~ 1 + System,  # Fixed effects
      random = ~ vsr(Parent1, Gu = Ginv.blend ) +  # GCA for Parent1
                 vsr(Tester, Gu = Ginv.blend_tester) +   # GCA for Tester
                 vsr(Parent1_Tester, Gu = Ginv.blend_Ihybrid), #SCA effect
      #rcov = ~ diag(System),  # Allows different residual variances for each system
      weights = Weights,
      data = training_data)
    
  
    # Extract genetic effects
    GCA_female <- model$U$`u:Parent1`$BLUE
    GCA_male   <- model$U$`u:Tester`$BLUE
    SCA        <- model$U$`u:Parent1_Tester`$BLUE
    System_SCA <- model$U$`u:System:Parent1_Tester`$BLUE
    
    # Split validation data into Organic and Conventional
    org   <- validation_data %>% filter(System == "Organic")
    conv  <- validation_data %>% filter(System == "Conventional")
    
    
    # Predictions for Conventional system
    conv$Model1  <-   model$Beta$Estimate[1]+ GCA_female[conv$Parent1]
    conv$Model2  <-   model$Beta$Estimate[1]+ GCA_female[conv$Parent1] + GCA_male[conv$Tester]
    conv$Model3  <-   model$Beta$Estimate[1]+ GCA_female[conv$Parent1] + GCA_male[conv$Tester] + SCA[conv$Parent1_Tester]
  
    # Predictions for Organic system
    org_intercept   <- model$Beta$Estimate[1] + model$Beta$Estimate[2]  #conventional is the reference system
    org$Model1  <- organic_intercept + GCA_female[org$Parent1]
    org$Model2  <- organic_intercept + GCA_female[org$Parent1] + GCA_male[org$Tester]
    org$Model3  <- organic_intercept + GCA_female[org$Parent1] + GCA_male[org$Tester] + SCA[org$Parent1_Tester]
  
    
     
    # Combine organic and conventional predictions
    Predicted <- rbind(org, conv)
    
    # Compute Pearson and Spearman correlations for model accuracy
    Pearson_Model1  <- cor(Predicted$BLUP, Predicted$Model1, use = "complete.obs")
    Pearson_Model2  <- cor(Predicted$BLUP, Predicted$Model2, use = "complete.obs")
    Pearson_Model3  <- cor(Predicted$BLUP, Predicted$Model3, use = "complete.obs")
 
    
    spearman_Model1 <- cor(Predicted$BLUP, Predicted$Model1, use = "complete.obs", method = "spearman")
    spearman_Model2 <- cor(Predicted$BLUP, Predicted$Model2, use = "complete.obs", method = "spearman")
    spearman_Model3 <- cor(Predicted$BLUP, Predicted$Model3, use = "complete.obs", method = "spearman")
   
    
    # Store accuracy metrics for this fold and iteration
    accuracy_results <- rbind(accuracy_results, data.frame(
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

#estimated variance components
summary(model)

#plot residuals
plot( model )
anova(model, object2=NULL, type=1)

# Step 5: Plot Accuracies
library(tidyr)
library(ggplot2)

# Reshape accuracy_results for plotting
PA_results_long <- accuracy_results %>%
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

#Plot Predicted vs Actual values

library(tidyr)
plot_data  <- Predicted %>%
  pivot_longer(cols = c(BLUE, BLUP, Model1, Model2, Model3, Model4), 
               names_to = "Type", 
               values_to = "Value")

# Plot using ggplot
##Treant
# Create an index for x-axis
plot_data$Index <- seq_along(plot_data$Value)
library(ggplot2)
ggplot(plot_data, aes(x = Index, y = Value, color = System)) +
  geom_point(size = 3) +                    # Adds points for each value
  geom_line(aes(group = Type), size = 1) +   # Adds lines connecting points for each type
  labs(x = "Index", y = "Value", title = "Comparison of BLUEs and GEBVs") +
  theme_minimal() +
  scale_color_manual(values = c("#a6cee3","#33a02c","#1f78b4","#b2df8a","#fb9a99",
                                "#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a", "yellow"))


##Boxplot
library(ggplot2)
ggplot(plot_data, aes(x = System, y = Value, fill = Type)) +
  stat_boxplot(geom="errorbar") +
  geom_boxplot(outlier.colour="Black", outlier.shape=1,outlier.alpha = 0, #this removes outliers. 
               outlier.size=2, notch = FALSE,varwidth = FALSE) +
  labs(x = "Treatments", y = "Value", title = "Comparison of BLUES and GEBV") +
  theme_minimal() +
  #scale_fill_manual(values = c("BLUE" = "#1f78b4", "BLUP" = "green", "GEBV" = "#e31a1c")) +
scale_color_manual(values = c("#a6cee3","#33a02c","#1f78b4","#b2df8a","#fb9a99",
                              "#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a", "yellow"))


# Compute correlation for each system
cor_results <- Predicted %>%
  group_by(System) %>%
  summarise(R_value = cor(BLUP, Model3, use = "complete.obs"))

print(cor_results)  # View correlation results

# Create scatter plot with regression lines and display R-values
ggplot(Predicted, aes(x = Model3, y = BLUP, color = System)) +
  geom_point(size = 2, alpha = 0.7) +  # Scatter points
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +  # Regression line
  facet_wrap(~ System) +  # Separate plots for each system
  theme_minimal() +
  labs(title = "Predicted vs. Actual BLUPs Across Systems",
       x = "Predicted Value",
       y = "Actual BLUP") +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 12, face = "bold")) +  # Style facet labels
  geom_text(data = cor_results, aes(x = min(Predicted$Model3), 
                                    y = max(Predicted$BLUP), 
                                    label = paste("R =", round(R_value, 3))),
            hjust = 0, vjust = 1, size = 5, color = "black")


                                  ####################################
                                  ## Predict Untested Testcrosses   ##
                                  ##                                ##
                                  ####################################

### Predict GCA of new inbreds
# Step 1: Load new testcross dataset
New_crosses <- read.csv("Validation_Hybrids.csv")  # Modify the path as needed

#check if all inbreds have genotypic data
New_parents <- unique(New_crosses$Parent1)
length(New_parents)
New_parents <- intersect(New_crosses$Parent1, rownames(GRM_All))
length(New_parents) #3 were removed because they were not genotyped

New_crosses <- New_crosses %>%
  filter(Parent1 %in% New_parents)

## Create interaction term in the new
New_crosses$Parent1_Tester <- factor(interaction(New_crosses$Parent1, New_crosses$Tester, sep = ":"))

# Merge tested and new to get common parents
Tested_parents <- unique(Pheno_org$Parent1)
Combined_parents <- union(New_parents, Tested_parents)
length(Combined_parents)

# Step 2: Subset the GRM for the combined list of parents
GRM_New <- GRM_All[Combined_parents, Combined_parents, drop = FALSE]
GRM_New[1:5, 1:5]
GRM_New <- as.matrix(GRM_New)

# Step 3: Predict GCA for untested inbreds using the new GRM
GCA_female <- GRM_New[New_crosses$Parent1, Tested_parents ] %*% model$U$`u:Parent1`$BLUE
New_crosses$GCA_female <- GCA_female

# Step 3: Retrieve GCA effects for the selected testers
GCA_male <- model$U$`u:Tester`$BLUE[New_crosses$Tester]
New_crosses$GCA_male <- GCA_male

New_crosses$Inter_Conv <- model$Beta$Estimate[1] 
New_crosses$Inter_Org <- model$Beta$Estimate[1] + model$Beta$Estimate[2]

New_crosses$Pred_Yield_Conv = New_crosses$Inter_Conv + New_crosses$GCA_female + New_crosses$GCA_male  
New_crosses$Pred_Yield_Org  = New_crosses$Inter_Org  + New_crosses$GCA_female + New_crosses$GCA_male  

# View the results
head(New_crosses)
write.csv(New_crosses, "Yield_New_crosses.csv")
