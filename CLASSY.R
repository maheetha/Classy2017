#!/usr/bin/Rscript


########################################################################
########################################################################
######## This Script is the Classy Pipeline. It takes all the ##########
######## extracted features, combines them, and runs random   ##########
######## random forest to calculate category probabilities    ##########
######## The final output probability for each gene is output ##########
########################################################################
########################################################################
########################################################################
######## Read training data can be derived from PanCan data or #########
######## you can input existing training data. In this script ##########
######## you must provide the flag "-d" at the end of your args ######## 
######## list in order to derive training data. ########################
########################################################################
########################################################################


library('kernlab')
args = commandArgs(trailingOnly = TRUE)

######## Read in all other data ################
traingenes = read.table(args[1], header = T, stringsAsFactors = F)
mafdata = read.table(args[2], header = T, stringsAsFactors = F)
cnadata = read.table(args[3], header = T, stringsAsFactors = F)
expdata = read.table(args[4], header = T, stringsAsFactors = F)


######## Create one data frame. ################
entire_data = merge(mafdata, cnadata, by =c('Gene'), all = TRUE)
entire_data = merge(entire_data, expdata, by = c('Gene'), all = TRUE)
entire_data = entire_data[complete.cases(entire_data),]
write.table(entire_data, 'entire_data', sep = "\t", quote = F, row.names = F)
entire_data = subset(entire_data, entire_data[,14] != entire_data[,15])


######## extract features #######################
trunc = entire_data[,2]/entire_data[,5]
recur = entire_data[,3]/entire_data[,5]
domain = entire_data[,4]/entire_data[,5]
totalsamps =  rowSums(entire_data[6:13]) + 1
cp1 = entire_data[,6]/totalsamps
cp2 = entire_data[,7]/totalsamps
cp3 = entire_data[,8]/totalsamps
cp4 = entire_data[,9]/totalsamps
cp5 = entire_data[,10]/totalsamps
cp6 = entire_data[,11]/totalsamps
cp7 = entire_data[,12]/totalsamps
cp8 = entire_data[,13]/totalsamps
meanmutate = entire_data[,14]/(entire_data[,14] + entire_data[,15] + 1)

######## Formatting Test and Training Data ########

data = as.data.frame(cbind(entire_data$Gene, trunc, recur, domain, cp1, cp2, cp3, cp4, cp5,  cp6, cp7, cp8, meanmutate))
names(data) = c('Gene', 'trunc', 'recur', 'domain', 'cp1', 'cp2', 'cp3', 'cp4', 'cp5','cp6', 'cp7', 'cp8', 'meanmutate')

if (length(args) != 6){
	trainingdata = as.matrix(traingenes[,2:13])
	class = as.numeric(traingenes$class)
} else {
	train = as.data.frame(subset(data, data[,1] %in% traingenes[,1]))
	train = merge(train, traingenes, by = c('Gene'), all = TRUE)
	train = train[complete.cases(train), ]
	trainingdata = as.matrix(train[,2:13])
	class(trainingdata) <- "numeric"
	class = as.numeric(train[,14])
	write.table(train, 'trainingdata', quote = F, sep = "\t", row.names = F)
}

testingdata = as.matrix(data[,2:13])
class(testingdata) <- "numeric"


######## Implement Random Forest ################

library('randomForest')
model = randomForest(trainingdata, class, keep.forest=TRUE)
importance(model)
prediction = predict(model, testingdata, type = "response")
data$randomForest = prediction

######## Output ########
write.table(data, args[5], quote = F, sep = "\t", row.names = F)

