
# Load the data
Data <- read.csv("Combined_Dataset.csv")

Data$YR <- as.factor(Data$YR)
Data$HYB <- as.factor(Data$HYB)
Data$GROUP <- as.factor(Data$GROUP)
Data$LOC <- as.factor(Data$LOC)
Data$REP <- as.factor(Data$REP)

Data$Range <- as.numeric(Data$Range)
Data$Row <- as.numeric(Data$Row)
Data$DTA <- as.numeric(Data$DTA)
Data$DTS <- as.numeric(Data$DTS)
Data$PHT <- as.numeric(Data$PHT)
Data$EHT <- as.numeric(Data$EHT)
Data$RLD <- as.numeric(Data$RLD)
Data$SLD <- as.numeric(Data$SLD)
Data$MAH <- as.numeric(Data$MAH)
Data$GWT <- as.numeric(Data$GWT)
Data$TWT <- as.numeric(Data$TWT)
Data$YLD <- as.numeric(Data$YLD)

Data$MST <- as.numeric(Data$MST)
Data$PROT <- as.numeric(Data$PROT)
Data$FIB <- as.numeric(Data$FIB)
Data$STA <- as.numeric(Data$STA)
Data$ASH as.numeric(Data$ASH
Data$S <- as.numeric(Data$STA)


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
  labs(x="trait", y="Count", title="Distribution of Yield")



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
ggplot(Loc_Blups, aes(x=LOC, y= BLUPs, fill= LOC)) + 
  geom_boxplot()

ggplot (Loc_Blups, aes(x = HYB, y = BLUPs)) +
  geom_boxplot() +
  labs(title = "Boxplot of estimated BLUPs by hybrid") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

