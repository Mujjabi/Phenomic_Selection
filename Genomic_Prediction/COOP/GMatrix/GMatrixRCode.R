##############################################
##  Apple Genomic Selection                 ##
##  Obtaining the G matrix and its inverse  ##
##############################################

# Install and load the R library ASRgenomics
install.packages("ASRgenomics")
library(ASRgenomics)

# Set working directory
setwd("C:/Users/JaneCohen/Documents/MyProjects")

# Loading Genotypic data (from R library)
markerFULL <- geno.apple

# View dimensions and the first 5 rows and 5 columns
dim(markerFULL)
markerFULL[1:5,1:5]

# Check for missing marker values
summary(colSums(is.na(markerFULL)))

# Filtering molecular matrix for downstream analyses
M_filter <- qc.filtering(M=markerFULL, maf=0.05, marker.callrate=0.2,
                         ind.callrate=0.20, impute=FALSE, plots=TRUE)

# View dimensions and the first records of filtered matrix
dim(M_filter$M.clean) 
M_filter$M.clean[1:5,1:5]

# Calculate the G matrix
Ghat <- G.matrix(M=M_filter$M.clean, method='VanRaden', na.string=NA)$G

# View the first 8 rows and 8 columns of Ghat
Ghat[1:8,1:8]

# Calculate the G-inverse matrix
Ginv<-solve(Ghat)

# View the first 8 rows and 8 columns of Ginv
Ginv[1:8,1:8]

# Determinant of Ghat
det(Ghat)

# Blending Ghat with the Identity, using p=0.02
Ghat.blend <- G.tuneup(G=Ghat, blend=TRUE, pblend=0.02)$Gb
Ghat.blend[1:8,1:8]

# Determinant and inverse of Ghat.blend
det(Ghat.blend)
Ginv.blend <- G.inverse(G=Ghat.blend)$Ginv
Ginv.blend[1:8,1:8]

# Bending Ghat
Ghat.bend <- G.tuneup(G=Ghat, bend=TRUE, eig.tol=1e-03)$Gb
Ghat.bend[1:8,1:8]

# Determinant and inverse of Ghat.bend
det(Ghat.bend)
Ginv.bend <- G.inverse(G=Ghat.bend)$Ginv
Ginv.bend[1:8,1:8]

# The G matrix variants
Ghat[1:8,1:8]
Ghat.blend[1:8,1:8]
Ghat.bend[1:8,1:8]

# Generate the three-column sparse form of the matrix Ginv.blend 
Ginv.sparse <- G.inverse(G=Ghat.blend, sparse=TRUE)$Ginv

# Compare Ginv.sparse and Ginv.blend
head(Ginv.sparse)
Ginv.blend[1:3,1:3]

# Save Ghat.blend & Ginv.sparse as RData file
save(Ghat.blend, Ginv.sparse, file = "Gmat.RData")
