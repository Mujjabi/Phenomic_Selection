

## Converting excel SNP data to hapmap format for TASSEL

# 1. Install and Load Necessary Libraries
# Install any necessary libraries for reading, transforming, and writing data:

install.packages("dplyr")

#2. Read Your Data

library(dplyr)

# Load the dataset
Markers <- read.csv("ALL_COOP_SGS_Markers.csv", header = TRUE, stringsAsFactors = FALSE)
map <- read.csv("Maize_pMap.csv", header = TRUE, stringsAsFactors = FALSE)

# Check the structure of the data
head(Markers)


# Step 2: Merge the files based on SNP names
hapmap_data <- map %>%
  inner_join(Markers, by = "Name")  # Assuming both files have an "SNP" column


# Create HapMap format columns

library(dplyr)
hapmap_data2 <- hapmap_data %>%
  mutate(
    rs = Name,                # SNP IDs
    alleles = "N/N",               # Replace with actual alleles if known
    chrom = CHR,                     # Assign mock chromosome numbers
    pos = POS,      # Assign mock positions
    strand = "+",                  # Default strand
    assembly = "NA", center = "NA", protLSID = "NA",
    assayLSID = "NA", panelLSID = "NA", QCcode = "NA") %>%
  select(rs, alleles, chrom, pos, strand, assembly, center, protLSID, assayLSID, panelLSID, QCcode)


# Add genotype data (convert symbols to HapMap format)
genotype_data <- Markers[, -1] %>%
  mutate_all(~ case_when(
    . == "A" ~ "AA",
    . == "T" ~ "TT",
    . == "C" ~ "CC",
    . == "G" ~ "GG",
    . == "Y" ~ "CT",
    . == "R" ~ "AG",
    . == "K" ~ "GT",
    . == "M" ~ "AC",
    . == "failed" ~ "NN",
    TRUE ~ NA_character_
  ))

#populate the alleles column in your HapMap data by examining the genetic 
#data and extracting the two unique alleles for each SNP.

hapmap_data3 <- hapmap_data2 %>%
  rowwise() %>%
  mutate(
    alleles = paste(
      unique(na.omit(grep("^[AGCT]$", Markers[which(Markers$Name == rs), -1], value = TRUE))),
      collapse = "/")) %>%
  ungroup()



# Combine HapMap metadata with genotype data
hapmap_data4 <- cbind(hapmap_data3, genotype_data)



#Steps to Convert Chromosome Values
#strip the "CHR" prefix and leading zeros from the chromosome column in R. Here's an example:

hapmap_final <- hapmap_data4 %>%
  filter(grepl("^CHR", chrom)) %>%  # Keep only rows where chrom starts with 'CHR'.
  mutate(chrom = as.integer(gsub("CHR", "", chrom)))  # Remove 'CHR' and convert to integer
# 24 SNPs were removed since they were labelled as scaffold or NAs

# View the first few rows
head(hapmap_final)   

#Note
# Make sure there are no empty cells anywhere, otherwise, the file wont be recognized by TASSEL


#Export final file

write.table(hapmap_final, file = "hapmap.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(hapmap_final, "COOP.hapmap", sep = "\t", quote = FALSE, row.names = FALSE)

write.csv(hapmap_final, "COOP_hapmap.csv", row.names = FALSE)



