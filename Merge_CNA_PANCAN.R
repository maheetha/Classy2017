#!/usr/bin/Rscript

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### the purpose of this script is to develop a Pan Can Merged Version #### 
#### it takes in all of the *.output files in the cna directory #### #### 
#### and combines the counts to create a PANCAN sample count files #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 


args = commandArgs(trailingOnly = TRUE)
files = list.files(pattern = glob2rx("*.output"))
merged = 1

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### merging all of the files that have *.output ending ##################
##########################################################################
##########################################################################

for (x in 1:length(files)){

	cancer = files[x]
	tomerge  = read.table(cancer, header = F, stringsAsFactors = F)

	one = paste(cancer, "0_05", sep = ".")
	two = paste(cancer, "05_1", sep = ".")
	three = paste(cancer, "1_15", sep = ".")
	four = paste(cancer, "15_2", sep = ".")
	six = paste(cancer, "2_25", sep = ".")
	seven = paste(cancer, "25_3", sep = ".")
	eight = paste(cancer, "3_35", sep = ".")
	nine = paste(cancer, "35_4", sep = ".")

	names(tomerge) = c('Gene', one, two, three, four,  six, seven, eight, nine)
	
	if (merged == 1){
		merged = tomerge
	
	} else {
		merged = merge(merged, tomerge, by = c('Gene'), all = TRUE)

	}

}

merged[is.na(merged)] = 0
write.table(merged, args[1], col.names = F, quote = F, row.names = F, sep = "\t")



#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### merging counts across newly merged PANCAN file  #####################
##########################################################################
##########################################################################

merged = read.table(args[1], header = F, stringsAsFactors = F)
multiplication = 1


##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### This is a file provided by MuSic that helps filter. The idea ##### ##### ##### 
##### is to only count samples for a gene from cancer types only in ##### ##### ##### 
##### which the gene is significantly mutated in. This way we don't ##### ##### ##### 
##### accumulate noise from cancer types where the gene isn't ##### ##### ##### ##### 
##### significantly mutated. This is completely optional. File format ##### ##### ##### 
##### is described in the comments section ##### ##### ##### ##### ##### ##### ##### ## 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 


if (args[3] != 'NULL'){
	cancermatrix = read.table(args[3], header = T, stringsAsFactors = F)
	merged = merged[order(merged$V1), ]
	cancermatrix = cancermatrix[order(cancermatrix$Genes),]
	merged = subset(merged, merged$V1 %in% cancermatrix$Genes)
	cancermatrix = subset(cancermatrix, cancermatrix$Genes %in% merged$V1)
	merged = merged[order(merged$V1), ]
	cancermatrix = cancermatrix[order(cancermatrix$Genes),]
	multiplication = cancermatrix[,2:ncol(cancermatrix)]
}


##### Merging Counts ######3
onevec = rowSums(merged[seq(2, ncol(merged), 8)]*multiplication)
twovec = rowSums(merged[seq(3, ncol(merged), 8)]*multiplication)
threevec = rowSums(merged[seq(4, ncol(merged), 8)]*multiplication)
fourvec = rowSums(merged[seq(5, ncol(merged), 8)]*multiplication)
fivevec = rowSums(merged[seq(6, ncol(merged), 8)]*multiplication)
sixvec = rowSums(merged[seq(7, ncol(merged), 8)]*multiplication)
sevenvec = rowSums(merged[seq(8, ncol(merged), 8)]*multiplication)
eightvec = rowSums(merged[seq(9, ncol(merged), 8)]*multiplication)


##### Writing to Final Table ##### 
finalmerged = as.data.frame(cbind(merged$V1, onevec, twovec, threevec, fourvec, fivevec, sixvec, sevenvec, eightvec))
names(finalmerged) = c('Gene', 'count0_05', 'count05_1', 'count1_15', 'count15_2', 'count2_25', 'count25_3', 'count3_35', 'count35_4')
write.table(finalmerged, args[2], row.names = F, sep = "\t", quote = F)

	


