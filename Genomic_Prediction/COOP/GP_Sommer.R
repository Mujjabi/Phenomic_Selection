# Load necessary libraries
library(sommer)

# Example data
pheno_data <- data.frame(
  Hybrid = c("H1", "H2", "H3", "H4", "H5", "H6"),
  Parent1 = c("Inbred1", "Inbred2", "Inbred3", "Inbred1", "Inbred2", "Inbred3"),
  Tester = c("Tester1", "Tester1", "Tester1", "Tester2", "Tester2", "Tester2"),
  y = c(10.5, 12.3, 11.8, 9.7, 10.2, 11.5)  # Phenotypic values (e.g., yield)
)

# Simulate marker data
set.seed(123)
n_inbreds <- 3
n_testers <- 2
n_markers <- 100

# Marker data for inbreds and testers separately
Z_inbreds <- matrix(sample(0:2, n_inbreds * n_markers, replace = TRUE), nrow = n_inbreds)
Z_testers <- matrix(sample(0:2, n_testers * n_markers, replace = TRUE), nrow = n_testers)

rownames(Z_inbreds) <- c("Inbred1", "Inbred2", "Inbred3")
rownames(Z_testers) <- c("Tester1", "Tester2")

# Compute separate GRMs for inbreds and testers
K_inbreds <- A.mat(Z_inbreds)
K_testers <- A.mat(Z_testers)

# Ensure Parent1_Tester interaction term is a factor with correct levels
pheno_data$Parent1_Tester <- factor(interaction(pheno_data$Parent1, pheno_data$Tester, sep = ":"))

# Create an identity matrix for the interaction term with proper levels
n_hybrids <- length(unique(pheno_data$Parent1_Tester))
I_hybrids <- diag(n_hybrids)
rownames(I_hybrids) <- levels(pheno_data$Parent1_Tester)
colnames(I_hybrids) <- levels(pheno_data$Parent1_Tester)

# Fit the model using sommer
model <- mmer(
  y ~ 1,  # Fixed effect (intercept only)
  random = ~ vs(Parent1, Gu = K_inbreds) + vs(Tester, Gu = K_testers) + vs(Parent1_Tester, Gu = I_hybrids),
  data = pheno_data
)

# Extract GCA and SCA effects
GCA_female <- model$U$`u:Parent1`$y  # GCA for Parent1 (inbreds)
GCA_male   <- model$U$`u:Tester`$y # GCA for Tester (testers)
SCA        <- model$U$`u:Parent1_Tester`$y  # SCA for Parent1:Tester (interaction)

# Match the levels in the prediction
predicted_performance <- model$Beta$Estimate + 
                         GCA_female[pheno_data$Parent1] +
                         GCA_male[pheno_data$Tester] +
                         SCA[pheno_data$Parent1_Tester]

# Add predicted performance to the dataset
pheno_data$PredictedPerformance <- predicted_performance

predictive_ability <- cor(pheno_data$y, pheno_data$PredictedPerformance, use = "complete.obs")



### Predict GCA for new inbreds without phenotypic data

# Simulate genotypic data for new inbreds
Z_new <- matrix(sample(0:2, 2 * n_markers, replace = TRUE), nrow = 2)
rownames(Z_new) <- c("NewInbred1", "NewInbred2")

# Extend marker matrix and recompute GRM
Z_extended <- rbind(Z, Z_new)
K_extended <- A.mat(Z_extended)

# Fit model with extended genomic relationship matrix
model_extended <- mmer(
  y ~ 1,  
  random = ~ vs(Parent1, Gu = K_extended) + vs(Tester, Gu = K) + vs(Parent1_Tester, Gu = I_hybrids),
  data = pheno_data
)

# Extract GCA for all inbreds, including new ones
GCA_female_extended <- model_extended$u$Parent1

# Print GCA estimates for new inbreds
print("Predicted GCA for new inbreds:")
print(GCA_female_extended[c("NewInbred1", "NewInbred2"), ])

# Predict hybrid performance for new inbreds crossed with existing testers
new_hybrid_performance <- model_extended$beta$`(Intercept)` + 
  GCA_female_extended["NewInbred1", ] + 
  GCA_male["Tester1"]  # Example: NewInbred1 x Tester1

print(paste("Predicted performance for NewInbred1 x Tester1:", round(new_hybrid_performance, 3)))