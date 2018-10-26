library(caret)


#inputFL = "/home/mwalter/statanalysis/compiledData.txt"
inputFL='statistics_wo_sameProt_hugeCasesExcl.txt';

data_all <- read.delim(inputFL, sep="\t", header=TRUE)
#head(data_all) #- visually checking the data
nrow(data_all) #-  574410
ncol(data_all) #-  45
str(data_all)


#ModelFL = args[2] #-- output file 1
outModelFL = "RFregression.model" #/home/mwalter/statanalysis/testoutput/result1.model"
outDIR = dirname(outModelFL)
resultFL = paste(outDIR,"/result.txt", sep='') #-- output file 2

a=substr(data_all$homologPair,1,4)
b = substr(data_all$qryPairID,1,4)
index = (a!=b)
#which (index == F)
data_all = data_all[index,]
nrow(data_all) # 61653



data_all$aveCC = (data_all$CC_1+ data_all$CC_2)/2

data_all$logEVal1 = log10(data_all$EVal1)
data_all$logEVal2 = log10(data_all$EVal2)

data_all[data_all$logEVal1=="-Inf"]<- -450 #xue
data_all[data_all$logEVal2=="-Inf"]<- -450

data_all$avelogEVal = (data_all$logEVal1 + data_all$logEVal2)/2


#> min(data_all$logEVal1)
#[1] -319.3662

#optional: check which gives better r^2:
data_all$logLAL_Q1 = log10(data_all$LAL_Q1)
data_all$logLAL_H1 = log10(data_all$LAL_H1)
data_all$logLAL_Q2 = log10(data_all$LAL_Q2)
data_all$logLAL_H2 = log10(data_all$LAL_H2)
#data_all[data_all=="-Inf"]<-"-450" # no zeroes expected so this line  unneeded

data_all$avelogLAL_Q = (data_all$logLAL_Q1 + data_all$logLAL_Q2)/2
data_all$avelogLAL_H = (data_all$logLAL_H1 + data_all$logLAL_H2)/2
data_all$avelogLAL = (data_all$avelogLAL_Q + data_all$avelogLAL_H)/2


#--- remove the queries with  <=4 interfacial residues (based on numInt_Q1 <= 4 or  numInt_Q2 <= 4).
#data = data [ data_all$Kd != 'na', ]
data_all = data_all[ data_all$numInt_Q1 > 4, ]
data_all = data_all[ data_all$numInt_Q2 > 4, ]

#data_all$numInt = (data_all$numInt_Q1+data_all$numInt_Q2)/2

data_all$aveSID = (data_all$SID1 + data_all$SID2)/2
data_all$avePositives = (data_all$Positives1 + data_all$Positives2)/2
data_all$aveSimilarity = (data_all$Similarity1+data_all$Similarity2)/2
data_all$avehhpredProb = (data_all$hhpredProb1 + data_all$hhpredProb2)/2
data_all$avefrac_LAL = (data_all$frac_LAL1 + data_all$frac_LAL2)/2



keep1 = c( "avelogEVal"     ,    "aveSID"   ,   "avePositives"  ,   "aveSimilarity"  , 
           "avehhpredProb"  ,    "avelogLAL_Q" ,   "avelogLAL_H"      ,   "avefrac_LAL"    ,   "aveCC"   )

keep_xue = 
c("logEVal1","logEVal2",
"SID1","SID2","Positives1","Positives2",
"Similarity1","Similarity2","hhpredProb1","hhpredProb2",
"hhpred_pVal1", "hhpred_pVal2",
"hhpred_SS1", "hhpred_SS2",
"LAL_Q1","LAL_H1","LAL_Q2","LAL_H2",
"frac_LAL1","frac_LAL2",
"aveCC")


#keep1 = c( "logEVal1"       ,  "logEVal2",
#          "SID1"           ,  "SID2"          ,  "Positives1"   ,   "Positives2",
#          "Similarity1"    ,  "Similarity2"   ,  "hhpredProb1"  ,   "hhpredProb2",
#          "LAL_Q1"         ,  "LAL_H1"        ,  "LAL_Q2"       ,   "LAL_H2",
#          "frac_LAL1"      ,  "frac_LAL2"     ,  "aveCC"   )


#keep2 = c( "logEVal1"       ,  "logEVal2",
#          "SID1"           ,  "SID2"          ,  "Positives1"   ,   "Positives2",
#          "Similarity1"    ,  "Similarity2"   ,  "hhpredProb1"  ,   "hhpredProb2",
#          "logLAL_Q1"      ,  "logLAL_H1"     ,  "logLAL_Q2"    ,   "logLAL_H2",
#          "frac_LAL1"      ,  "frac_LAL2"     ,  "aveCC"   )



#keep1 (using the normal LAL values) gives an r^2 of: 0.2571438
#keep2 (using the log -  LAL values) gives an r^2 of: 0.2486042
# therefore use keep1 for now!

data=data_all[, names(data_all) %in% keep_xue]

#--- visually check each feature vs. aveCC
pairs(data)

#--- plot the medians
ggplot(aes(x=factor(LAL_Q1),y= aveCC), data=data_all) +geom_point(color='grey')+ stat_summary(fun.y = "median", geom='point', color='red') + ggtitle("LAL_Q1 vs aveCC") 

ggplot(aes(x=factor(LAL_Q2),y= aveCC), data=data_all) +geom_point(color='grey')+ stat_summary(fun.y = "median", geom='point', color='red') 

ggplot(aes(x=factor(LAL_H1),y= aveCC), data=data_all) +geom_point(color='grey')+ stat_summary(fun.y = "median", geom='point', color='red') 

ggplot(aes(x=factor(LAL_H2),y= aveCC), data=data_all) +geom_point(color='grey')+ stat_summary(fun.y = "median", geom='point', color='red') 

frame()
ggplot(aes(x=factor(avelogLAL),y= aveCC), data=data_all) +geom_point(color='grey')+ stat_summary(fun.y = "median", geom='point', color='red') 



#data = data[ , data$Eval1 ] <- log(data[ ,data$Eval1 ])
#data = data[ , data$Eval2 ] <- log(data[ ,data$Eval2 ])
#
#data[data=="-Inf"]<-"-450"

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


#------ multiple linear regression
#train_control <- trainControl(method = "repeatedcv", number=5, repeats = 1,    returnData = T,  savePredictions = "final")
#set.seed(134)
#model <- train(aveCC ~ ., data=data  , trControl=train_control, method ="lmStepAIC")

#------ random forest
#-- random forest
train_control <- trainControl(method = "repeatedcv", number=5, repeats = 1,     returnData = T,  savePredictions = "final", search ="grid")
mtry <- sqrt(ncol(data ) ) # -- suggested value. The number of variables for each tree.
mtry = floor(mtry)
#tunegrid <- expand.grid(.mtry=c(1:16))
tunegrid <- expand.grid(.mtry=mtry)
set.seed(134)
model <- train(aveCC ~ ., data=data , trControl=train_control, method="rf", ntree = 200, tuneGrid=tunegrid ) #-- R-sq = 0.76, mtry = 1


# Check model content
model$results
#      mtry     RMSE  Rsquared      RMSESD  RsquaredSD
#1 4.582576 0.150678 0.6643608 0.002304845 0.009767977

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



