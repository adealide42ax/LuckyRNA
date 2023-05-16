library(RNAseqStat2)
library(dplyr)
args <- commandArgs(T)

params <- read.table(args[1], header=T, sep = "," )
row_counts <- read.delim(args[2], header=T, sep = "\t", row.names=1, check.names=FALSE)
GO_anno <- read.table(args[3], sep = "\t",header = T,check.names = F)
GO_description <- read.table(args[4], sep = "\t",header = T,check.names = F, quote = "")
KEGG_anno <- read.table(args[5], sep = "\t",header = T,check.names = F)
KEGG_description <- read.table(args[6], sep = "\t", header = T, check.names = F, quote = "")

group <- colnames(row_counts)
group <- strsplit(group,split = "_")
group <- sapply(group,"[",1)
group_list <- group

run <- function(row_counts, condition1, condition2, FC, FDR, output) {
  
  colnames <- colnames(row_counts)
  split_colnames <- strsplit(colnames, "_")
  # 获取行名中包含 condition1 或 condition2 的列的名称
  cols <- sapply(split_colnames, function(x) x[1] %in% c(condition1, condition2))
  # 使用这些列的名称选择这些列
  row_count <- row_counts[, cols]
  group_list <- group[group %in% c(condition1, condition2)]
  
  data_i <- Create_DEGContainer(expMatrix = row_count,
                                groupInfo = group_list,
                                caseGroup = condition2,
                                species = "unknow",
                                GOTERM2GENE = GO_anno,
                                GOTERM2NAME = GO_description,
                                KEGGTERM2GENE = KEGG_anno,
                                KEGGTERM2NAME = KEGG_description)
  cutFC(data_i) <- as.numeric(FC)
  cutFDR(data_i) <- as.numeric(FDR)
  cutFC(data_i)
  cutFDR(data_i)
  data_o <- runALL(object = data_i,dir = output,GO = TRUE,KEGG = TRUE)
}

for (i in 1:nrow(params)) {
  
  condition1 <- params[i, 1]
  condition2 <- params[i, 2]
  FC <- params[i, 3]
  FDR <- params[i, 4]
  output <- paste(condition1,condition2,FC,FDR, sep = "_")
  run(row_counts,condition1, condition2, FC, FDR,output)
}
