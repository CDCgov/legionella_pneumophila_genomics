#!/usr/bin/env Rscript

library(caret)
library(randomForest)

load(file ="/scicomp/home/xxh5/LabGroupData/Serogroup_2018/pipelines/20180904_rFModel_update.RData")

args <- commandArgs(trailingOnly=TRUE)

dataInfo = read.csv(args[1], header=TRUE,stringsAsFactor=FALSE,colClasses =c("character"))

#print(args[1])

dataSet = subset(dataInfo, select = -c(SgType))

predictData = predict(rFModel_update, newdata = dataSet)

outFile = paste0(args[2],"-predictResults.txt")

finalOutFile = paste0(args[3],outFile)
#write.csv(predictData, file="resultsPredict.csv")
write.csv(predictData, file=finalOutFile)
