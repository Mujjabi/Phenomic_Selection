
# Load the data
PSData <- read.csv("2022YieldTrial_PaulScott.csv")
PSData <- read.csv("Ours.csv")

PSData$PED <- as.factor(PSData$PED)
PSData$HYB <- as.factor(PSData$HYB)
PSData$LOC <- as.factor(PSData$LOC)
PSData$REP <- as.factor(PSData$REP)
PSData$MST <- as.numeric(PSData$MST)
PSData$GWT <- as.numeric(PSData$GWT)
PSData$TWT <- as.numeric(PSData$TWT)
PSData$YLD <- as.numeric(PSData$YLD)

#Cleaning out Outliers

boxplot(PSData[8:11], plot=TRUE)$out

outliers <- function(x) {
  
  Q1 <- quantile(x, probs=.25) # calculate first quantile
  Q3 <- quantile(x, probs=.75) # calculate third quantile
  iqr = Q3-Q1           # calculate inter quartile range
  
  upper_limit = Q3 + (iqr*1.5)  #Calculate upper limit
  lower_limit = Q1 - (iqr*1.5)  #calculate lower limit
  
  x > upper_limit | x < lower_limit   # return true or false
}

remove_outliers <- function(PSData, cols = names()) {  # for loop to traverse in columns vector
  for (col in cols) {
    PSData <- PSData[!outliers(PSData[[col]]),] # remove observation if it satisfies outlier function
  }
  PSData
}

#Remove outliers 
PSData <- remove_outliers(PSData, c(8:11))
boxplot(PSData[8:11], plot=TRUE)$out


#Yield Distribution
library(ggplot2)
ggplot(PSData, aes(x=YLD))  + 
  geom_histogram(fill="cadetblue", alpha=0.7, binwidth=10) + 
  labs(x="Grain Yield [bu/acre]", y="Count", title="Yield Distribution")

ggplot(PSData, aes(x=MST))  + 
  geom_histogram(fill="coral4", alpha=0.7, binwidth=0.5) + 
  labs(x="Moisture at Harvest [%]", y="Count", title="Moisture Distribution")

ggplot(PSData, aes(x=TWT))  + 
  geom_histogram(fill="blue", alpha=0.7, binwidth=0.5) + 
  labs(x="Test Weight [lbs/bu]", y="Count", title="Test Weight Distribution")


ggplot(PSData, aes(x=MST))  + 
  geom_histogram(fill="cadetblue", alpha=0.5, binwidth=0.5) + 
  labs(x="Moisture At Harvest", y="Count", title="Distribution of Yield")



traits <- colnames(PSData)[8:ncol(PSData)]



# Load the lme4 package
library(lme4)

# Fit the linear mixed effects model
model <- lmer(YLD ~  (1|HYB) + (1|LOC/REP) + (1|LOC:HYB), data= PSData )

# Check the model summary
summary(model)


# Estimate BLUPs for each hybrid
library(dplyr)
Means <- PSData %>%
  group_by(HYB) %>%
  summarize(Mean = mean(YLD, na.rm=TRUE))#Means of each HYB across locations


blups <- ranef(model) #These are the blups. 
predicted_yield <-  blups$HYB + mean(PSData$YLD) #this is the Blup-adjusted mean

#Another Way
varComp <-as.data.frame(VarCorr(model,comp="vcov")) #function calculates estimated variances between random-effects terms in a mixed-effects model  blup = coef(model)$Entry
YieldBlups = coef(model)$HYB   #this also extracts BLup-adjusted means

BLUPs <-as.data.frame(cbind(Means,predicted_yield,YieldBlups)) #compare the three
colnames(BLUPs)[3] <- "Predicted"
colnames(BLUPs)[4] <- "YieldBlup"
write.csv (BLUPs, "BLUPs.csv")
#Rank Hybrids from top to bottom
Ranking <- BLUPs[order(BLUPs$YieldBlup, decreasing = TRUE), ]

# Extract the top 10 hybrids
Top_ten  <- head(Ranking, n = 10)
print(Top_ten)



#BLUPs By location
varComp <-as.data.frame(VarCorr(model,comp="vcov")) #function calculates estimated variances between random-effects terms in a mixed-effects model  blup = coef(model)$Entry
BlUP2 = coef(model)  
Hyb_LOC <- BlUP2$`LOC:HYB`
write.csv (Hyb_LOC, "Hyb_LOC.csv")

Loc_Blups <- read.csv("Hyb_LOC.csv")

library(ggplot2)
ggplot(Loc_Blups, aes(x=LOC, y= BLUP, fill= LOC)) + 
  geom_boxplot()

ggplot (Loc_Blups, aes(x = HYB, y = BLUPs)) +
 geom_boxplot() +
 labs(title = "Boxplot of estimated BLUPs by hybrid") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))





## FOR OUR HYBRIDS AND THE CHECK
PSData <- read.csv("Ours.csv")

PSData$PED <- as.factor(PSData$PED)
PSData$HYB <- as.factor(PSData$HYB)
PSData$LOC <- as.factor(PSData$LOC)
PSData$REP <- as.factor(PSData$REP)
PSData$MST <- as.numeric(PSData$MST)
PSData$GWT <- as.numeric(PSData$GWT)
PSData$TWT <- as.numeric(PSData$TWT)
PSData$YLD <- as.numeric(PSData$YLD)

#Cleaning out Outliers

boxplot(PSData[8:11], plot=TRUE)$out

outliers <- function(x) {
  
  Q1 <- quantile(x, probs=.25) # calculate first quantile
  Q3 <- quantile(x, probs=.75) # calculate third quantile
  iqr = Q3-Q1           # calculate inter quartile range
  
  upper_limit = Q3 + (iqr*1.5)  #Calculate upper limit
  lower_limit = Q1 - (iqr*1.5)  #calculate lower limit
  
  x > upper_limit | x < lower_limit   # return true or false
}

remove_outliers <- function(PSData, cols = names()) {  # for loop to traverse in columns vector
  for (col in cols) {
    PSData <- PSData[!outliers(PSData[[col]]),] # remove observation if it satisfies outlier function
  }
  PSData
}

#Remove outliers 
PSData <- remove_outliers(PSData, c(8:11))
boxplot(PSData[8:11], plot=TRUE)$out


# Calculate the mean and standard deviation of the selected individuals
library(Rmisc)
library(ggplot2)
summary_data <- summarySE(PSData, measurevar = "YLD", groupvars = NULL)

ggplot(PSData, aes(x = factor(LOC), y = YLD, fill = LOC)) + 
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "blue") +
  stat_boxplot(geom="errorbar") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1 ))+
  scale_fill_manual(values = c("red", "green", "blue", "orange", "gray")) +
  labs(title = "Boxplot of Selected Hybrids",
     x = "Tested Location",
     y = "Grain Yield [Bu/Acre]")


ggplot(PSData, aes(x = HYB, y = YLD, fill = HYB)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "blue") +
  stat_boxplot(geom="errorbar") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(title = "Boxplot of Selected Hybrids",
       x = "UIUC Hybrids",
       y = "Grain Yield [Bu/Acre]")

ggplot(PSData, aes(x = HYB, y = YLD, fill = LOC)) +
  geom_boxplot() +
  stat_boxplot(geom="errorbar") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(title = "Boxplot of Selected Hybrids",
       x = "UIUC Hybrids",
       y = "Grain Yield [Bu/Acre]")

# Load the lme4 package
library(lme4)
library(tidyr)

# Fit the linear mixed effects model
model <- lmer(YLD ~  (1|HYB) + (1|LOC/REP) + (1|LOC:HYB), data= PSData )

# Check the model summary
summary(model)

#Extract Blups
varComp <-as.data.frame(VarCorr(model,comp="vcov")) #function calculates estimated variances between random-effects terms in a mixed-effects model  blup = coef(model)$Entry
BLUPs = coef(model)   #this also extracts BLup-adjusted means
Hyb_LOC <- BLUPs$`LOC:HYB`

HYB <- rownames(Hyb_LOC) # get the hybrid names
# match the intercepts with the hybrids
YLD <- Hyb_LOC[, 1]  # the first column contains the intercepts
BLUPs <- data.frame(HYB, YLD)

BLUPs <- separate(BLUPs, col = HYB, into = c("LOC", "HYB"), sep = ":")

#Rank Hybrids from top to bottom
Blup_mean <- aggregate(YLD ~ HYB, data = BLUPs, mean)
Ranking <- Blup_mean[order(Blup_mean$YLD, decreasing = TRUE), ]
print(Ranking)

