fit <- glm(R ~., data = Reliability2, family = binomial)
summary(fit)
plot(fit, which = 1:2)


r <- randomForest(R~.,data=chiapas_train, nodesize=30)
plot(chiapas_test$R,predict(r,chiapas_train$R),col="green")
points(x,y)


rf <- randomForest(R~., data=chiapas_train, proximity=TRUE, nodesize=15) 
print(rf)

p1 <- predict(rf, chiapas_train)
confusionMatrix(p1, chiapas_train$R)

p2 <- predict(rf, Validation_Phase)
confusionMatrix(p2,Validation_Phase$R)




## Improving model perfomance 

#Stepwise regression  
step(fit,direction = "backward") 
step(fit,direction = "forward")


#models based on various information criteria.
#k=2 yields the AIC criterion, and k=log (n) refers to BIC.

step(fit,k=log(nrow(Reliability2)))  
step(fit,k=2)  

fit2 = step(fit,k=2,direction = "backward")
summary(fit2)

model_stepwise <- glm(formula = R ~ HYB + YEAR + FS + CC + Density + Rotation + CEC + pH + OM + EstNRel + Su + P + Ca + Mg + K + Na + SOC + TN + SoilCN + InorgN
                      + Bray1P + POMC + POMN + POMCN + PMN + POXC + Prcp + Tmax + Tmin + TGDD, data = Reliability2)
summary(model_stepwise)
plot(model_stepwise, which = 1:2)

#Data Split into training and testing 75/25% 
set.seed(1234)
train_index <- sample(seq_len(nrow(Reliability2)), size = 0.75*nrow(Reliability2))
chiapas_train<-Reliability2[train_index, ]
chiapas_test<-Reliability2[-train_index, ]

library(rpart)
chiapas.rpart<- rpart(fit, data=chiapas_train)
# install.packages("rpart.plot")
library(rpart.plot)
rpart.plot(chiapas.rpart, digits=3)

#fallen leaves. 
library(rpart.plot)
rpart.plot(chiapas.rpart, digits = 4, fallen.leaves = T, type=3, extra=101) 

#order and rules of split
library(rattle)
fancyRpartPlot(chiapas.rpart, cex = 0.8)

#Make a prediction
Pred1  <-predict(chiapas.rpart,chiapas_test )
summary(Pred1)

out.sol <-cbind(chiapas_test,Pred1)
out.sol2 <- out.sol[c(1,2,19)]
out.sol2 <- data.frame(out.sol2)

#Evaluating Model Performance 
# 1. Accuracy. This measures the correlation btn predicted and actual values. 
cor(Pred1, chiapas_test$R)  

#The metric mean absolute error measures 
#To measure the distance between predicted value and the true value, we can use a measurement called mean absolute error (MAE). 

MAE<-function(obs, pred){mean(abs(obs-pred))}
MAE(chiapas_test$YIELD, Pred1)  #a value of 0.9 shows yield difference is 0.9 given the range of yield in the dataset
summary(chiapas_test$YIELD) 


library(RWeka) # This can improve the model, reduce the mean absolute error. 
mlb.m5<-M5P(model_stepwise, data=chiapas_train)
summary(mlb.m5)
mlb.p.m5<-predict(mlb.m5, chiapas_test)
cor(mlb.p.m5, chiapas_test$YIELD)  
MAE(chiapas_test$YIELD, mlb.p.m5) # it didn't improve our model



##### NEURAL NETWORKS
library(caret)
chiapas <- read.csv("initial_data.csv", stringsAsFactors = T)
chiapas <-chiapas [-c(5)]

Validation <- read.csv("to_be_solved.csv", stringsAsFactors = T)
Validation <-Validation [-c(4)]

### Convert categorical variables to numeric dummy variables since neural networks only accept numeric
Chiapas_NN  <- sapply(chiapas, unclass)
Validation <- sapply(Validation, unclass)

Chiapas_NN <- data.frame(Chiapas_NN)
Validation <- data.frame(Validation)

## Normalizing the data 
normalize <- function(x) {return((x - min(x)) / (max(x) - min(x)))}
chiapas_norm <-as.data.frame(lapply(Chiapas_NN[2:19], normalize))
Validation_norm  <-as.data.frame(lapply(Validation  [2:18], normalize))

chiapas_norm$System <- as.numeric(chiapas_norm$System )

#Assigning training and testing sets
set.seed(1234)
train_index <- sample(seq_len(nrow(chiapas_norm)), size = 0.75*nrow(chiapas_norm))
chiapas_train<-chiapas_norm[train_index, ]
chiapas_test<-chiapas_norm[-train_index, ]


# install.packages("neuralnet")
library(neuralnet)
model_NN <-neuralnet(System ~ YIELD + Year + Tillage + Planting + Nitrogen + Potassium + 
                       Slope + Clay + CEC + PH + prcp + srad + tmax + tmin + vp, data=chiapas_train)
plot(model_NN)

#Model performance 
NN_Pred <-compute(model_NN, chiapas_test[, c(2:18)])
pred_results <-NN_Pred$net.result
Accuracy <- cor(pred_results, chiapas_test$System)
Accuracy
pred1 <- ifelse(pred_results[,1] > 0.8,"Landrance","Hybrid")
table(pred1)
prop.table(table(pred1))

## Adding hidden layers 
library(neuralnet)
model_NN2 <-neuralnet(System ~ YIELD + Year + Tillage + Planting + Nitrogen + Potassium + 
                        Slope + Clay + CEC + PH + prcp + srad + tmax + tmin + vp, data=chiapas_train, hidden = 3)
plot(model_NN2)

#Model performance 
NN_Pred2 <-compute(model_NN2, chiapas_test[, c(2:18)])
pred_results <-NN_Pred2$net.result
cor(pred_results, chiapas_test$System)
pred1 <- ifelse(pred_results[,1] > 0.6,"Landrance","Hybrid")
table(pred1)
prop.table(table(pred1))

### Validate model using to_be_solved data 

NN_Pred3 <-compute(model_NN, Validation_norm)
pred_results3 <-NN_Pred3$net.result

out.sol3<-cbind(Validation_norm,pred_results3)

out.sol3$System <-ifelse(pred_results3[,1] > 0.35,"Hybrid","Landrace")
out.sol3<-cbind(Validation$Field_ID,out.sol3)
colnames(out.sol3)[1] <- "Field_ID"

outsol4<-out.sol3[,c(1,23)]
write.csv(outsol4, file ="Mujjabi_Test_NN5.csv", row.names = F)