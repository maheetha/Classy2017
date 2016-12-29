#!/usr/bin/Rscript

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### the purpose of this script is to develop a Pan Can Merged Version #### 
#### it takes in all of the *.rsem files in the maf directory #### #### 
#### and combines the counts to create a PANCAN sample count files #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

args = commandArgs(trailingOnly = TRUE)
files = list.files(pattern = glob2rx("*.rsem"))
merged = 1


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### merging all of the files that have *.rsem ending ####################
##########################################################################
##########################################################################
for (x in 1:length(files)){
	
	cancer = files[x]
	tomerge  = read.table(cancer, header = T, stringsAsFactors = F)

	if (merged == 1){
		merged = tomerge
	} else {
		temp = tomerge[,2:3]
		merged = as.data.frame(cbind(merged,temp)) 
	}


}

merged[is.na(merged)] = 0
write.table(merged, args[1], quote = F, row.names = F, sep = "\t")

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### merging counts across newly merged PANCAN dataframe  ################
##########################################################################
##########################################################################
merged= read.table(args[1], header = F, stringsAsFactors = F)
multiplication = 1


##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### This is a file provided by MuSic that helps filter. The idea ##### ##### ##### 
##### is to only count samples for a gene from cancer types only in ##### ##### ##### 
##### which the gene is significantly mutated in. This way we don't ##### ##### ##### 
##### accumulate noise from cancer types where the gene isn't ##### ##### ##### ##### 
##### significantly mutated. This is completely optional. File format ##### ##### ##### 
##### is described in the comments section ##### ##### ##### ##### ##### ##### ##### ## 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

if (length(args) == 3){
	cancermatrix = read.table(args[3], header = T, stringsAsFactors = F)
	merged = merged[order(merged[,1]), ]
	cancermatrix = cancermatrix[order(cancermatrix$Genes),]
	merged = subset(merged, merged[,1] %in% cancermatrix$Genes)
	cancermatrix = subset(cancermatrix, cancermatrix$Genes %in% merged[,1])
	merged = merged[order(merged[,1]), ]
	cancermatrix = cancermatrix[order(cancermatrix$Genes),]
	multiplication = cancermatrix[,2:ncol(cancermatrix)]
}


one = rowSums(merged[seq(2, ncol(merged), 2)]*multiplication)
two = rowSums(merged[seq(3, ncol(merged), 2)]*multiplication)


finalmerged = as.data.frame(cbind(merged$Gene, one, two))
names(finalmerged) = c('Gene', 'meannonmutate', 'meanmutate')
write.table(finalmerged, args[2], row.names = F, sep = "\t", quote = F)
