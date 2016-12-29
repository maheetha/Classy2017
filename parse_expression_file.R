#!/usr/bin/Rscript

#####################################################################
#####################################################################
### For ease, this script tranposes the expression files so     #####
### that sample names are the rows and gene names are the columns ###
#####################################################################
#####################################################################

args = commandArgs(trailingOnly = TRUE)
gene_expression_table = read.table(args[1], sep = "\t", stringsAsFactor = F)

tranposed_gene_expression_table = t(gene_expression_table)
tranposed_gene_expression_table[,1] = substr(tranposed_gene_expression_table[,1], 1, 12)

#####################################################################
#####################################################################
### Read in the mutated samples that were output from parse_maf.sh ##   
#####################################################################
#####################################################################
samples = read.table(args[2],  stringsAsFactors = F)
samples = as.data.frame(samples)


### Obtain number of colummns and initialize final table ###
size = ncol(tranposed_gene_expression_table)
final_table = c()


#### A for loop that parses file line by line #####

for (x in 3:(size-1)){
  
  ## obtain gene and it's RSEM values for the samples ###
  genetable = as.data.frame(cbind(tranposed_gene_expression_table[,1],tranposed_gene_expression_table[,x]))
  gene = as.vector(tranposed_gene_expression_table[,x])
  
  ### obtaining the name of the gene ####
  gene_line = gene[1]
  index = which(strsplit(gene_line, "")[[1]]=="|")
  gene_id = substr(gene_line, 1, index-1)
  
  
  #subsetting for the gene expression profile for this gene
  gene_expression_profile = subset(samples, (samples$V2 == gene_id & (samples$V3 == "Missense_Mutation" 
                           | samples$V3 == "Nonsense_Mutation" | samples$V3 == "Nonstop_Mutation" | 
                           samples$V3 == "Frame_Shift_Ins" | samples$V3 == "Frame_Shift_Del" | 
                             samples$V3 == "Splice_Site" )))

  # if no  samples for this gene are mutated, then next
  if (length(datasubset[,1]) == 0) {
    line = paste(gene_id , 1, 1, sep = "\t")
    final_table = c(final_table,line)
    next 
  } 
  
  #getting those samples that are from mutated samples
  len = length(datasubset[,1])
  mutated = c()
  nonmutated = c()
  
  ## Filtering for mutated versus nonmutated ####
  
  mutated = subset(genetable, genetable[,1] %in% gene_expression_profile[,1] )
  nonmutated = subset(genetable, !(genetable[,1] %in% gene_expression_profile[,1]))
  
  #Separate mutated and nonmutated
  nonmutated = nonmutated[-1,]
  
  # if there are no mutated samples, then next
  if (length(mutated) == 0){
    line = paste(gene_id , 1, 1, sep = "\t")
    final_table = c(final_table,line)
    next
  }
  
  # calculate median scores
  meannonmutate = median(as.numeric(as.character(nonmutated$V2)))
  meanmutate = median(as.numeric(as.character(mutated$V2)))
  line = paste(gene_id , meannonmutate, meanmutate,  sep = "\t")
  final_table = c(final_table, line)
  
}


## Writing to output ####
nonmutatename = paste(args[3], "meannonmutate", sep = ".")
mutatename = paste(args[3], "meanmutate", sep = ".")
filename = paste(args[3], "rsem", sep = ".")
title = paste("Gene", nonmutatename, mutatename, sep = "\t")
final_table = c(title, final_table)
write.table(final_table, filename, row.names = F, col.names = F, quote = F, sep = "\t")




