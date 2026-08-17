library(limma)
library(DESeq2)
library("RColorBrewer")
library("ggrepel")
library(amap)
library(gplots)
library(ggplot2)
library(BiocParallel)
library(pheatmap)
library(dplyr)
#install.packages("tidyr")  # 重新安装
library(tidyr)  

data<-read.table("E:/数据分析/ZZ5-LHB-mix-20260613/RNA-seq/RNA-seq all_counts.txt")

YXD=data.frame(TN1=data[,6],TN2=data[,7],TN3=data[,8],
                 WT1=data[,9],WT2=data[,10],WT3=data[,11])
rownames(YXD)=data[,1]
group=c("TN","TN","TN","WT","WT","WT")
YXD_2 <- YXD[rowSums(YXD)>10,] 
conditions<-group
batch<-factor(c(1,1,1,1,1,1))
sample<-data.frame(conditions,batch)   
rownames(sample)<-colnames(YXD_2)
#产生DESeq数据集并计算标准化因子

ddsFullCountTable <- DESeqDataSetFromMatrix(countData =YXD_2,colData = sample,  design= ~ conditions)
dds <- DESeq(ddsFullCountTable)

normalized_counts <- counts(dds, normalized=TRUE)
#normalized_counts_mad <- apply(normalized_counts, 1, mad) 
#normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
#write.table(normalized_counts, file="E:/数据分析/ZZ5-LHB-mix-20260613/Data.normalized.xls",quote=F,sep="\t", row.names=T, col.names=T)
'''
## log转换后的结果并输出
rld <- rlog(dds, blind=FALSE)
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]
write.table(rlogMat, file="E:/DESeq2.normalized.rlog.xls",quote=F, sep="\t", row.names=T, col.names=T)
'''
######样本相关性热图绘制及PCA分析
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
pearson_cor <- as.matrix(cor(normalized_counts, method="pearson"))
hc <- hcluster(t(normalized_counts), method="pearson")
heatmap.2(pearson_cor, Rowv=as.dendrogram(hc), symm=T, trace="none",
          col=hmcol, margins=c(11,11), main="Pearson correlation")
#pca_data <- plotPCA(dds, intgroup=c("conditions"), returnData=T, ntop=3000)

sampleA = "TN"
sampleB = "WT"
contrastV <- c("conditions", sampleA, sampleB)
res <- results(dds,contrast=contrastV)
#res$padj[is.na(res$padj)] <- 1  #校正后p-value为NA的赋值为1
res <- cbind(ID=rownames(res),normalized_counts ,as.data.frame(res))
res <- res[order(res$padj),]
res_de <- subset(res, res$padj<0.05, select=c('ID','WT1','WT2','WT3',"TN1","TN2","TN3",'log2FoldChange', 'padj'))
res_de_up <- subset(res_de, res_de$log2FoldChange>=1)
res_de_dw <- subset(res_de, res_de$log2FoldChange<=(-1))



write.table(res,"E:/数据分析/ZZ5-LHB-mix-20260613/ALL-data-nor-log.xls",sep="\t", quote=F, row.names=T,col.names = T)
write.table(res_de_up,"E:/数据分析/ZZ5-LHB-mix-20260613/up.xls",sep="\t", quote=F, row.names=T)
write.table(res_de_dw,"E:/数据分析/ZZ5-LHB-mix-20260613/down.xls",sep="\t", quote=F, row.names=T)

cut_off_qvalue = 0.05
cut_off_logFC = 1
res$Sig <- ifelse(res$padj < cut_off_qvalue & 
                    abs(res$log2FoldChange) > cut_off_logFC,
                  ifelse(res$log2FoldChange > cut_off_logFC ,'Enriched','Depleted'),'Similar')
res <- data.frame(res)
tmp <- res %>% drop_na(Sig)
table(res$Sig)
###add name delete #
#gene_tmp <- c("nuoJ","mukF","nlpE","glpK","ssrS","rrsD","rrsG","rrsH","rrsE","rrsB","ssrA","rrsA","rrsC")
#gene_tmp <- data.frame(gene_tmp)
#gene_tmp$geneList <- gene_tmp$gene_tmp
#ID <- res$ID
#tmp <- res %>% left_join(gene_tmp,by = c("ID" = "gene_tmp"))

ggplot(res, aes(x = res$log2FoldChange, y = -log10(res$padj), colour=Sig)) +
  geom_point(alpha=0.4, size=1.5) +
  scale_color_manual(values=c("blue","red","black")) + 
  xlim(c(-10, 10)) + 
  ylim(c(0,150))+ 
  geom_vline(xintercept=c(-cut_off_logFC,cut_off_logFC),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_qvalue),
             lty=4,col="black",lwd=0.8) +
  labs(x="Log2(FC)",
       y="-log10(Padj)")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank() 
  )# +  
#  geom_text_repel(aes(label=geneList), fontface="bold",
#                  color="black", box.padding=unit(0.35, "lines"),
#                  point.padding=unit(5, "lines"),segment.color ="black",max.overlaps=Inf)
#dev.off()

#logCounts <- log2(res$baseMean+1)
#logFC <- res$log2FoldChange
#FDR <- res$padj
#plot(logFC, -1*log10(FDR), col=ifelse(FDR<=0.01, "red", "black"),
#     xlab="logFC", ylab="-log10 Padj", main="Differentially expressed gene",
#     pch=1)

###choose 1
res_de_up_sorted <- res_de_up %>% arrange(padj)
res_de_dw_sorted <- res_de_dw %>% arrange(padj)
res_de_up_top20_id <- as.vector(head(res_de_up_sorted$ID,20))
res_de_dw_top20_id <- as.vector(head(res_de_dw_sorted$ID,20))
#choose 2
#res_de_up_top20_id <- as.vector(head(res_de_up$ID,20))
#res_de_dw_top20_id <- as.vector(head(res_de_dw$ID,20))

res_de_top20 <- c(res_de_up_top20_id, res_de_dw_top20_id)   
res_de_top20_expr <- normalized_counts[rownames(normalized_counts) %in% res_de_top20,]
res_de_top20_expr <-YXD_2[rownames(YXD_2) %in% res_de_top20,]
pheatmap(res_de_top20_expr, cluster_row=T, scale="row", annotation_col=sample)
library(pheatmap)


library(clusterProfiler)
gene2ko <- download_KEGG(organism = "eco", keggType = "KEGG")$KEGG2Gene
colnames(gene2ko) <- c("KO", "gene_id")
# 2. 获取KO-eco通路映射
ko2path <- download_KEGG(organism = "eco", keggType = "KEGG")$Pathway2Gene
colnames(ko2path) <- c("pathway_id", "KO")
# 3. 三表合并：gene-KO-eco通路
gene_ko_path <- merge(gene2ko, ko2path, by = "KO")
# 4. 导出csv本地保存
write.csv(gene_ko_path, "Ecoli_gene_KO_ecoPath.csv", row.names = F)








up_genes <- rownames(res_de_up)
Tn2<-read.table("E:/数据分析/ZZ7-LHB-mix-20260701/转录组/KEGG GO/up-kegg-forcluster.txt")

Tn2<-read.table("E:/数据分析/ZZ7-LHB-mix-20260701/转录组/KEGG GO/down-kegg-forcluster.txt")
#symbols <-bitr(Tn5, fromType = "ALIAS", toType = "ENTREZID",OrgDb = 'org.EcK12.eg.db')

Tn5<-Tn2[,1]
kegg <- enrichKEGG(
  gene = Tn5,# up_genes
  organism = 'eco', 
  pAdjustMethod = 'fdr', 
  pvalueCutoff = 0.95, #0.2  #0.4  #0.3  #0.6
  qvalueCutoff = 0.95, 
)
kegg@result
dotplot(kegg, x = "FoldEnrichment", color = "pvalue") +
  scale_color_distiller(palette = "RdYlBu", direction = -1)
write.csv(kegg,"E:/数据分析/ZZ7-LHB-mix-20260701/转录组/KEGG GO/up-KEGG-GC-v2.csv")

Tn2<-read.table("E:/数据分析/ZZ7-LHB-mix-20260701/转录组/KEGG GO/up-go-forcluster.txt")

Tn2<-read.table("E:/数据分析/ZZ7-LHB-mix-20260701/转录组/KEGG GO/down-go-forcluster.txt")
Tn5<-Tn2[,1]
library(org.EcK12.eg.db)
#up_genes_entrez <- bitr(up_genes, fromType = "SYMBOL",toType = "ENTREZID", OrgDb = org.EcK12.eg.db)
ego <- enrichGO(
  gene          = Tn5,
  keyType = "ENTREZID",
  OrgDb         = org.EcK12.eg.db,
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.95,
  qvalueCutoff  = 0.95,
  readable      = TRUE)
dotplot(ego)#,showCategory = 50)
#barplot(go_up, showCategory =8, split = "ONTOLOGY") +
#  facet_grid(ONTOLOGY ~ ., scales ="free")
write.csv(ego,"E:/数据分析/ZZ7-LHB-mix-20260701/转录组/KEGG GO/up-GO-gc-v2.csv")#GO-result-ALL.csv")


# 1. 定义要筛选的通路名称（确保与结果中的 Description 完全一致）
target_pathways <- c("Starch and sucrose metabolism - Escherichia coli K-12 MG1655","Lysine degradation - Escherichia coli K-12 MG1655","Quorum sensing - Escherichia coli K-12 MG1655","Pentose and glucuronate interconversions - Escherichia coli K-12 MG1655","Microbial metabolism in diverse environments - Escherichia coli K-12 MG1655","Propanoate metabolism - Escherichia coli K-12 MG1655","ABC transporters - Escherichia coli K-12 MG1655","Glycerophospholipid metabolism - Escherichia coli K-12 MG1655","Glycerolipid metabolism - Escherichia coli K-12 MG1655")
target_pathways <- c("Aminoacyl-tRNA biosynthesis - Escherichia coli K-12 MG1655")
target_pathways <- c("quorum sensing","aspartate family amino acid catabolic process","polyol metabolic process","response to pH","glycerol metabolic process","small molecule catabolic process","cellular response to stress","oxidoreductase activity")
target_pathways <- c("monoatomic cation transmembrane transporter activity","active monoatomic ion transmembrane transporter activity","inorganic molecular entity transmembrane transporter activity","monoatomic ion transmembrane transporter activity","inorganic cation transmembrane transporter activity")

# 2. 提取富集结果数据框
res <- kegg@result  ####!!!!!!!!!!!!!!!!!!!!!!!xiugai

res <- ego@result  ####!!!!!!!!!!!!!!!!!!!!!!!xiugai

# 3. 筛选出目标通路（如果某个通路不在结果中，则会被忽略）
sub_res <- res[res$Description %in% target_pathways, ]

# 检查筛选结果
print(sub_res[, c("Description", "pvalue", "GeneRatio", "BgRatio")])

# 如果 sub_res 为空，说明通路名称不匹配，请检查是否有大小写或空格差异
# 可以先用 unique(res$Description) 查看实际名称

# 4. 计算 FoldEnrichment（若没有该列）
if(!"FoldEnrichment" %in% colnames(sub_res)){
  gene_ratio <- as.numeric(sapply(strsplit(sub_res$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))
  bg_ratio <- as.numeric(sapply(strsplit(sub_res$BgRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))
  sub_res$FoldEnrichment <- gene_ratio / bg_ratio
}

# 5. 绘图（使用 ggplot2）
library(ggplot2)

# 确保数值列为 numeric
sub_res$FoldEnrichment <- as.numeric(sub_res$FoldEnrichment)
sub_res$pvalue <- as.numeric(sub_res$pvalue)
sub_res$Count <- as.numeric(sub_res$Count)

# 按 FoldEnrichment 降序排列（使图上通路从高到低显示）
sub_res <- sub_res[order(sub_res$FoldEnrichment, decreasing = TRUE), ]

p4 <- ggplot(sub_res, aes(x = FoldEnrichment, y = reorder(Description, FoldEnrichment))) +
  geom_point(aes(size = Count, color = pvalue)) +
  scale_color_distiller(palette = "RdYlBu", direction = -1) +   # 红-黄-蓝渐变
  labs(
    x = "Fold Enrichment",
    y = "KEGG Pathway",
    color = "P value",
    size = "Gene count"
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10))

print(p1+p2+p3+p4)






go_data<-read.table("E:/数据分析/ZZ7-LHB-mix-20260701/转录组/up-KEGG-GC-v2.txt")#,# header = FALSE,
       # sep = "\t", stringsAsFactors = FALSE, fill = TRUE, quote = " ")
go_data<-read.csv2("E:/数据分析/ZZ7-LHB-mix-20260701/转录组/up-KEGG-GC-v2.csv",# header = FALSE,
                    sep = ",", stringsAsFactors = FALSE, fill = TRUE, quote = " ")
go_data1<-go_data[,-1]
#colnames(go_data1)<-go_data1[1,]
#go_data1<-go_data1[-1,]
#  c("ONTOLOGY", "ID", "Description", "GeneRatio", "BgRatio", "RichFactor", "FoldEnrich", "Score", "pvalue", "p.adjust", "qvalue", "geneID", "Count")
go_data1$FoldEnrich <- as.numeric(go_data1$FoldEnrichment)
go_data1$pvalue <- as.numeric(go_data1$pvalue)
go_data1$Count <- as.numeric(go_data1$Count)

# 按 FoldEnrich 排序（使纵轴从高到低显示）
go_data <- go_data1[order(go_data1$FoldEnrich, decreasing = TRUE), ]

ggplot(go_data, aes(x = FoldEnrich, 
                         y = reorder(Description, FoldEnrich))) +  # 按富集倍数排序
  geom_point(aes(size = Count, color = pvalue)) +  # 点大小表示基因数
  scale_color_gradient(
    low = "#B22222",    # 显著（p小）→ 深红
    high = "#1E3A8A"    # 不显著（p大）→ 深蓝
  ) +
  # 如果您喜欢示例图中的红-黄-蓝渐变，可改为：
  # scale_color_distiller(palette = "RdYlBu", direction = -1) +
  labs(
    x = "Fold Enrichment",
    y = "GO Term",
    color = "P value",
    size = "Gene count"
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10))
p
