#!/usr/bin/env Rscript

# Example data  ./dataset/data_exp-hs-ddg-allr.dat
# compare this script with cal_aic_fit_dg.R
args <- commandArgs(TRUE)

if (length(args) != 2){
    print ("Usage: cal_cv_caret.R input_dataFL outputModelFL\n")
    quit();
}

library(caret)

inputFL =  '/home/lixue/PSHOMPPIv1.3/uploadData/test/statistics_qryIDpair_wo_sameProt.csv'; #args[1]
data_all <- read.delim(inputFL, sep=",", header=TRUE)
head(data_all) #- visually checking the data

outModelFL = args[2] #-- output file 1
outDIR = dirname(outModelFL)
resultFL = paste(outDIR,"/result.txt", sep='') #-- output file 2

#--- average (CC1 + CC2)
data_all$aveCC = (data_all$CC_1+ data_all$CC_2)/2;
data_all$logEval1 = log(as.numeric(data_all$Eval1 )  )
data_all$logEval2 = log(data_all$Eval2 )


#--- remove columns that we will not use

keep = c(   
 "numInt_Homolog1", "numInt_Homolog2", "logEVal1"         ,  "logEVal2",
  "SID1"          ,  "SID2"          ,  "Positives1"   ,   "Positives2",
  "Similarity1"   ,  "Similarity2"   ,  "hhpredProb1"  ,   "hhpredProb2",
  "LAL_Q1"        ,  "LAL_H1"        ,  "LAL_Q2"       ,   "LAL_H2",
  "frac_LAL1"     ,  "frac_LAL2"   , "aveCC"   );

data=data_all[, names(data_all) %in% keep]



#-- to check the colinearity between variables
#library(car)
#scatterplotMatrix(data,var.labels=colnames(data), diagonal=c("histogram"))


#------ training and predicting

train_control <- trainControl(method = "repeatedcv", number=10, repeats = 5,  returnData = T,  savePredictions = "final")

#---------------------- Train model -----------------------------------#
## check training methods http://topepo.github.io/caret/modelList.html
#---------------------- Train model -----------------------------------#

## Multiple Linear Regression with Stepwise Selection, default both directions.
set.seed(134)
model <- train(aveCC ~ ., data=data  , trControl=train_control, method ="lmStepAIC")

#-- random forest
#train_control <- trainControl(method = "repeatedcv", number=5, repeats = 1,  returnData = T,  savePredictions = "final", search ="grid")
#mtry <- sqrt(ncol(data ) ) # -- suggested value
#tunegrid <- expand.grid(.mtry=c(1:20))
##tunegrid <- expand.grid(.mtry=mtry)
#model <- train(aveCC ~ ., data=data , trControl=train_control, method="rf", tuneGrid=tunegrid ) #-- R-sq = 0.76, mtry = 1 

# Check model content
print(model)
plot(model)
names(model)
summary(model)
names(model$finalModel)
model$resample

#-- plot pred vs obs
obs = model$pred$obs
pred = model$pred$pred
min = min( c(obs, pred) )
max = max( c(obs, pred) )

plot(obs, pred,  xlim = c(min, max), ylim = c(min, max) ) 

#-- check residue plot
resi = pred - obs
plot(pred,resi)

#---
cor=cor(pred, obs) # RF: 0.47 
line = paste(inputFL, as.character(cor))
write(line, file=resultFL, append=T)

#--- save the model
save(model, file=outModelFL)

#--
line = paste (resultFL ,'and', outModelFL, 'generated')
print( line  )



