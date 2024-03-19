# This code converts the spc file obtained from the DA 7200 NIR machine after scanning the samples
# into an excel format that is readable. 

install.packages("prospectr")
#library(prospectr)
library(hyperSpec)  # this is actually the package that supports read.spc function, not prospectr 
# Read the SPC file into a data frame
F1Samples  <- read.spc("Phenomic Selection.spc")

  
#Convert file to data.frame
F1Samples_df <- data.frame(F1Samples)
F1Samples_df <- F1Samples_df[-c(1,2,144,145)]

# Write the data frame to a CSV file
write.csv(F1Samples_df, file = "F1 Spectral data.csv", row.names = TRUE)

asreml.license.status(quiet = FALSE, task = "checkout", json = "")



