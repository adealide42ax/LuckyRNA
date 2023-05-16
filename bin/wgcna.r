#!/usr/bin/env Rscript
library(WGCNA)
library(reshape2)
library(stringr)
options(stringsAsFactors = FALSE)
args<-commandArgs(T)
exprMat<-args[1]
dataExpr<-read.table(exprMat,sep='\t',row.names=1,header=T,quote="",comment="",check.names=F)
type="unsigned"
corType="pearson"
corFnc=ifelse(corType=="pearson",cor,bicor)
dataExpr_log<-log(dataExpr+1)
input_data <- as.data.frame(t(dataExpr_log))
#查看是否有离群样本
nGenes = ncol(input_data)
nSamples = nrow(input_data)
sampleTree = hclust(dist(input_data), method = "average")
png(file="sampleTree.png")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()
#软阈值计算
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(input_data, powerVector=powers, networkType=type, verbose=5)
cex1 = 0.9;
pdf(file="power.pdf")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
dev.off()
#平均连接度计算
pdf(file="mean_connect.pdf")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
 text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
power=sft$powerEstimate
power
#一步法构建网络
net = blockwiseModules(input_data, power = sft$powerEstimate,TOMType = "unsigned", minModuleSize = 30,reassignThreshold = 0,
                       mergeCutHeight = 0.25,numericLabels = TRUE, pamRespectsDendro = FALSE,saveTOMs = TRUE,
                       saveTOMFileBase = "rattan.TOM",verbose = 3)
table(net$colors)
#标签转化为颜色；
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
#动态剪切后的modules
unmergedColors = labels2colors(net$unmergedColors)
pdf(file="module_combined.pdf", width=12, height=12)
plotDendroAndColors(
  net$dendrograms[[1]],
  cbind(unmergedColors[net$blockGenes[[1]]], mergedColors[net$blockGenes[[1]]]),
  c("Dynamic Tree Cut" , "Merged colors"),
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)
dev.off()
#导出模块里的基因
MEList = moduleEigengenes(input_data, colors = unmergedColors)
MEs = MEList$eigengenes
moduleColors = labels2colors(net$colors)
modNames = substring(names(MEs), 3)
AllGenes =colnames(input_data)
for (module in modNames)
{
    modGenes = (mergedColors==module)
    modLLIDs = AllGenes[modGenes]
    fileName = paste(module, ".txt", sep="")
    write.table(as.data.frame(modLLIDs), file = fileName,row.names = FALSE, col.names = FALSE)
}
#导出TOM矩阵
TOM = TOMsimilarityFromExpr(input_data, power = sft$powerEstimate)
write.table(TOM,file="tom.cor.txt",quot=FALSE,sep="\t", row.names = TRUE,col.names = TRUE)