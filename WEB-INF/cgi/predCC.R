#!/usr/bin/env Rscript

# Li Xue
# Jan. 14th, 2017
# 
# Given the alignment quality, predict intrface conservatation in terms of CC  


#--
args = commandArgs(TRUE);

if (length(args) != 2){
   print ("Usage: predCC.R statFL_ori statFL_withPredCC");
   quit();
}

library(caret)
load('RFregression.model')

#---- read alignment quality
statFL = args[1]  #-- contains alignment quality
resultFL = args[2] #-- output file, i.e., statFL with predCC
data_ori = read.csv(statFL, comment.char = '#', sep ="\t") 
data = data_ori;


#---- calculate additional features
data$logEVal1 = log10(data$EVal1)
data$logEVal2 = log10(data$EVal2)
data[data$logEVal1=="-Inf"]<- -450 #xue
data[data$logEVal2=="-Inf"]<- -450

data$avelogEVal = (data$logEVal1 + data$logEVal2)/2

data$logLAL_Q1 = log10(data$LAL_Q1)
data$logLAL_H1 = log10(data$LAL_H1)
data$logLAL_Q2 = log10(data$LAL_Q2)
data$logLAL_H2 = log10(data$LAL_H2)

data$avelogLAL_Q = (data$logLAL_Q1 + data$logLAL_Q2)/2
data$avelogLAL_H = (data$logLAL_H1 + data$logLAL_H2)/2
data$avelogLAL = (data$avelogLAL_Q + data$avelogLAL_H)/2


#---- make predictions
pred_CC = predict(final_model, data)

#----
data_ori$pred_CC = sprintf("%.2f", pred_CC);
write.table(data_ori, file = resultFL, sep="\t", row.names = F)

quit()





