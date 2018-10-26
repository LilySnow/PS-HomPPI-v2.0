library(caret)


inputFL = "/home/mwalter/statanalysis/compiledData.txt"

data_all <- read.delim(inputFL, sep="\t", header=TRUE)
#head(data_all) #- visually checking the data

#ModelFL = args[2] #-- output file 1
outModelFL = "/home/mwalter/statanalysis/testoutput/result1.model"
outDIR = dirname(outModelFL)
resultFL = paste(outDIR,"/result1.txt", sep='') #-- output file 2

#--- remove homologs that have the same PDB IDs as the query proteins
a=substr(data_all$homologPair,1,4)
b = substr(data_all$qryPairID,1,4)
index = (a!=b)
#which (index == F)
data_all = data_all[index,]

#---

data_all$aveCC = (data_all$CC_1+ data_all$CC_2)/2

data_all$logEVal1 = log(data_all$EVal1)
data_all$logEVal2 = log(data_all$EVal2)
data_all[data_all=="-Inf"]<-"-450"

#> min(data_all$logEVal1)
#[1] -319.3662

#optional: check which gives better r^2:
data_all$logLAL_Q1 = log(data_all$LAL_Q1)
data_all$logLAL_H1 = log(data_all$LAL_H1)
data_all$logLAL_Q2 = log(data_all$LAL_Q2)
data_all$logLAL_H2 = log(data_all$LAL_H2)
data_all[data_all=="-Inf"]<-"-450"


#--- remove the queries with  <=4 interfacial residues (based on numInt_Q1 <= 4 or  numInt_Q2 <= 4).
#data = data [ data_all$Kd != 'na', ]
data_all = data_all[ data_all$numInt_Q1 > 4, ]
data_all = data_all[ data_all$numInt_Q2 > 4, ]



keep1 = c( "logEVal1"       ,  "logEVal2",
          "SID1"           ,  "SID2"          ,  "Positives1"   ,   "Positives2",
          "Similarity1"    ,  "Similarity2"   ,  "hhpredProb1"  ,   "hhpredProb2",
          "LAL_Q1"         ,  "LAL_H1"        ,  "LAL_Q2"       ,   "LAL_H2",
          "frac_LAL1"      ,  "frac_LAL2"     ,  "aveCC"   )


keep2 = c( "logEVal1"       ,  "logEVal2",
          "SID1"           ,  "SID2"          ,  "Positives1"   ,   "Positives2",
          "Similarity1"    ,  "Similarity2"   ,  "hhpredProb1"  ,   "hhpredProb2",
          "logLAL_Q1"      ,  "logLAL_H1"     ,  "logLAL_Q2"    ,   "logLAL_H2",
          "frac_LAL1"      ,  "frac_LAL2"     ,  "aveCC"   )



#keep1 (using the normal LAL values) gives an r^2 of: 0.2571438
#keep2 (using the log -  LAL values) gives an r^2 of: 0.2486042
# therefore use keep1 for now!

data=data_all[, names(data_all) %in% keep1]




#data = data[ , data$Eval1 ] <- log(data[ ,data$Eval1 ])
#data = data[ , data$Eval2 ] <- log(data[ ,data$Eval2 ])
#
#data[data=="-Inf"]<-"-450"

set.seed(134)
#sub_folds <- createFolds(y = prot_names, k=n_prot, list = TRUE, returnTrain =   TRUE)



#in_train <- holdout <- vector(mode = "list", length = length(sub_folds)) #      initialization

#row_index <- 1:nrow(data)


#for(i in seq(along = sub_folds)) {
#  ## Which subjects are in fold i
#  sub_in <- prot_names[sub_folds[[i]]]
#  ## which rows of the data correspond to those subjects
#  in_train[[i]] <- row_index[data$Complex %in% sub_in]
#  holdout[[i]]  <- row_index[!(data$Complex %in% sub_in)]
#}

#names(in_train) <- names(holdout) <- names(sub_folds)

#------ training and predicting
train_control <- trainControl(method = "repeatedcv", number=10, repeats = 5,    returnData = T,  savePredictions = "final")


#------ multiple linear regression
model <- train(aveCC ~ ., data=data  , trControl=train_control, method ="lmStepAIC")

#------ random forest
#-- random forest
#train_control <- trainControl(method = "repeatedcv", number=5, repeats = 1,     returnData = T,  savePredictions = "final", search ="grid")
##mtry <- sqrt(ncol(data ) ) # -- suggested value
#tunegrid <- expand.grid(.mtry=c(1:20))
##tunegrid <- expand.grid(.mtry=mtry)
#model <- train(aveCC ~ ., data=data , trControl=train_control, method="rf",     tuneGrid=tunegrid ) #-- R-sq = 0.76, mtry = 1


# Check model content
print(model)
#plot(model)
#names(model)
summary(model)
#names(model$finalModel)
#model$resample



#-- plot pred vs obs
obs = model$pred$obs
pred = model$pred$pred
min = min( c(obs, pred) )
max = max( c(obs, pred) )

plot(obs, pred,  xlim = c(min, max), ylim = c(min, max) )


#-- use ggplot to plot pred vs obs
library('ggplot2')
df = data.frame(obs,pred) #-- create a data frame
head(df)


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



