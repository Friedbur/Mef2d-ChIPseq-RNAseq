library(ensembldb) #Import all required libraries
library(tximport)
library(DESeq2)
library(vsn)
library(apeglm)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ensembldb)
library(tximport)
library(DESeq2)
library(vsn)
library(apeglm)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
organism="org.Mm.eg.db"
library(organism, character.only=TRUE)
library(dplyr)
library(scales)
library(gplots)

### Code from Class##
edb <- EnsDb("/projects/academic/courses/bch519/RNA-seq/Mus_musculus.GRCm38.100.sqlite") #Load Ensembl database of mouse annotation
Tx <- transcripts(edb, return.type="DataFrame")[9:8] # extract columns 8 and 9: transcipt annotations
gene_annotation <- genes(edb, return.type="DataFrame")

files <- c("./WT1/quant.sf","./WT2/quant.sf","./KO1/quant.sf","./KO2/quant.sf") 
names(files) <- paste0("sample", 1:4) 
txi <- tximport(files,type="salmon",tx2gene=Tx,ignoreTxVersion=TRUE) 

sampleTable <- data.frame(condition=factor(rep(c("WT","KO"), each=2)))
rownames(sampleTable) <- colnames(txi$counts)
dds <- DESeqDataSetFromTximport(txi,sampleTable,~condition)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds$condition <- relevel(dds$condition,ref="WT")
dds <- DESeq(dds)
res <- lfcShrink(dds,coef="condition_KO_vs_WT",type="apeglm") #to work remove low levels

row_match <- match(rownames(res),gene_annotation[,1])
res <- as.data.frame(cbind(res, gene_annotation[row_match,2]))
colnames(res)[6] <- "Gene"
##########################################

res <- as.data.frame(cbind(res, txi$abundance[rownames(res),])) #add abundance Values
colnames(res)[7] <- "WT1" #rename to WT1
colnames(res)[8] <- "WT2" #rename to WT2
colnames(res)[9] <- "KO1" #rename to KO1
colnames(res)[10] <- "KO2" #rename to KO2

sig <- subset(res,padj <= 0.01 &abs(log2FoldChange) >= 0.5 ) #filter for only significant genes

####### For Final part?????##

mef2d_targets <- read.csv("Mef2d_targets.csv") ## Read Mef2d target genes file
colnames(mef2d_targets) <- c("Gene", "Peak", "Distance") #assign col names to match
sig <- sig[sig$Gene %in% mef2d_targets$Gene, ] #finds overlap between  expression and binding
row_match <- match(sig$Gene, mef2d_targets$Gene) #Match gene order between sig and mef2d_targets: need to attach dist and peak
sig <- cbind( sig, mef2d_targets[row_match, c("Peak", "Distance")]) #append Peak and distance   from Mef2d
sig <- sig[order(sig$padj), ] # sort the table by the p values
write.csv(sig, file="Mef2d_final.csv") # save output

#for Final Med part ^^^
##########################################

down <- subset(sig, log2FoldChange < 0) #extract downregulated genes
down <- down[order(down$padj), ] #sort them
cat("downregulated genes:") #print the downregulated genes
print(down[1:10, c("Gene","log2FoldChange","padj","WT1","WT2","KO1","KO2")])

up <- subset(sig, log2FoldChange > 0) # now extract upregualted genes
up <- up[order(up$padj), ] # sort them by p value
cat("Upregulated genes:") # print top 10 upregulated genes
print(up[1:10, c("Gene","log2FoldChange","padj","WT1","WT2","KO1","KO2")])

up_genes <- up$Gene #store uop genes for analysis
down_genes <-down$Gene # store list of down genes for analysis
sig_bp_up <- enrichGO(gene=up_genes, universe= res$Gene,OrgDb=organism, keyType="SYMBOL",readable=TRUE, ont="BP", qvalueCutoff=0.1) #GO Biological Process enrichment for upregulated genes
sig_bp_down <- enrichGO(gene=down_genes, universe= res$Gene,OrgDb=organism, keyType="SYMBOL",readable=TRUE, ont="BP", qvalueCutoff=0.1) #GO Biological Process enrichment for downpregulated genes

top_up <- head(sig_bp_up@result[order(sig_bp_up@result$qvalue), ], 5) #get top 5 upregualted go results
top_down <- head(sig_bp_down@result[order(sig_bp_down@result$qvalue), ], 5) # extract top 5 downreg go results
cat("GO terms-Upregulated:") # Print go enrichment for upreg genes (top 5)
print(top_up[, c("ID","Description","qvalue")])
cat("GO terms-Downreulated:") # do the same for downreg go results)
print(top_down[, c("ID","Description","qvalue")])

write.csv(res, file="WT-KO-All.csv") #save full results
write.csv(sig, file="W:WT-KO-Sig.csv") #save significant results

save.image() # idk if I actually need this but have it from class so ;-;.

#DONE and hopefully works
