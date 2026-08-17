library(Seurat)
library(cowplot)
library(patchwork)
library(dplyr)
library(harmony)
Data1 <- read.table("E:/数据分析/WH-LHB-scRNA-20260602/0427-3.txt", sep="\t")
Data11 <- Data1[-1,-1]   
colnames(Data11) <- Data1[1,-1]
rownames(Data11) <- Data1[-1,1]
Data2 <- read.table("E:/数据分析/WH-LHB-scRNA-20260602/S0427-4.txt", sep="\t")
Data21 <- Data2[-1,-1]   
colnames(Data21) <- Data2[1,-1]
rownames(Data21) <- Data2[-1,1]
# Set up control object
ecoli24h_1 <- CreateSeuratObject(counts = t(Data11), project = "Data1", min.cells = 5)
VlnPlot(ecoli24h_1, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
ecoli24h_1$group <- "replicate1" 
ecoli24h_1_1 <- subset(ecoli24h_1, subset = nFeature_RNA > 30 & nFeature_RNA < 4000)##cell level  >50  <2500
ecoli24h_1_2 <- NormalizeData(ecoli24h_1_1, verbose = FALSE)
ecoli24h_1_3 <- FindVariableFeatures(ecoli24h_1_2, selection.method = "vst", nfeatures = 500)

ecoli24h_2 <- CreateSeuratObject(counts = t(Data21), project = "Data2", min.cells = 5)
VlnPlot(ecoli24h_2, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
ecoli24h_2$group <- "replicate2"
ecoli24h_2_1 <- subset(ecoli24h_2, subset = nFeature_RNA >20 & nFeature_RNA < 4000)##cell level
ecoli24h_2_2 <- NormalizeData(ecoli24h_2_1, verbose = FALSE)
ecoli24h_2_3 <- FindVariableFeatures(ecoli24h_2_2, selection.method = "vst", nfeatures = 500)

biofilm.anchors <- FindIntegrationAnchors(object.list = list(ecoli24h_1_3, ecoli24h_2_3), dims = 1:20)
biofilm.combined <- IntegrateData(anchorset = biofilm.anchors, dims = 1:20)
VlnPlot(biofilm.combined, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

#Perform an integrated analysis
DefaultAssay(biofilm.combined) <- 'integrated'
# Run the standard workflow for visualization and clustering
biofilm.combined_1 <- ScaleData(biofilm.combined, verbose = FALSE)
biofilm.combined_2 <- RunPCA(biofilm.combined_1, npcs = 30, verbose = FALSE)

VizDimLoadings(biofilm.combined_2, dims = 1:4, reduction = 'pca')
DimPlot(biofilm.combined_2, reduction = 'pca')
)#
DimHeatmap(biofilm.combined_2, dims = 1, cells = 500, balanced = TRUE)
P<-DimHeatmap(biofilm.combined_2, dims = 1:9, cells = 500, balanced = TRUE) & guides(fill = guide_colorbar(title = "Expression"))
P + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom",
        legend.box = "horizontal") &  # 横向放置
  guides(fill = guide_colorbar(barwidth = 20, barheight = 1.5, 
                               title.position = "top", 
                               title.hjust = 0.5))
ggsave("E:/数据分析/WH-LHB-scRNA-20260602/DimHeatmap_with_legend.pdf", plot = P, width = 16, height = 8, dpi = 150)
# Determine the 'dimensionality' of the dataset
#More approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
biofilm.combined_2 <- JackStraw(biofilm.combined_2, num.replicate = 100)
pbmc9 <- ScoreJackStraw(biofilm.combined_2, dims = 1:16)

JackStrawPlot(pbmc9, dims = 1:16)
ElbowPlot(pbmc9)

biofilm.combined_3 <- RunHarmony(biofilm.combined_2, group.by.vars = "orig.ident") 
DimPlot(biofilm.combined_3, group.by=NULL,reduction = 'pca',pt.size = 0.85)

# t-SNE and Clustering
biofilm.combined_4 <- RunUMAP(biofilm.combined_3, reduction = 'pca', dims = 1:5)
biofilm.combined_5 <- FindNeighbors(biofilm.combined_4, reduction = 'pca', dims = 1:5)
biofilm.combined_6 <- FindClusters(biofilm.combined_5, resolution = 0.17)   #0.6
# Visualization
DimPlot(biofilm.combined_6, reduction = 'umap', group.by = 'stim')#,cols = c("orange","blue")"red", "green"
DimPlot(biofilm.combined_6, reduction = 'umap', label = TRUE,split.by = "stim")
#
sub1 <- as.matrix(GetAssayData(biofilm.combined_6)[, WhichCells(biofilm.combined_6, ident = "3")])
epithelial_cells <- subset(biofilm.combined_6, subset = stim %in% c("3"))#####修改为用cluster选子群

sub2 <- as.matrix(GetAssayData(biofilm.combined_6, slot = "counts")[, WhichCells(biofilm.combined_6, ident = "2")])
write.table(data.frame(data.frame(sub1)), 
            'G:/YXD-scRNA-seq data/24h cluster2 gene-cell matrix.xlsx',sep = '\t', quote = FALSE, row.names = FALSE)
cells_cluster3 <- WhichCells(biofilm.combined_6, idents = 3)
num_cells <- length(cells_cluster3)
print(num_cells)
pbmc.markers <- FindAllMarkers(biofilm.combined_6, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
##part marker gene
pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top10
DoHeatmap(biofilm.combined_6, features = top10$gene) #+ NoLegend()
#biofilm.combined_6@assays$RNA@scale.data <- scale(biofilm.combined_6@assays$RNA@data, scale = TRUE)
#DoHeatmap(biofilm.combined_6, features = top10$gene,size = 5, angle = -50, hjust=0.8,slot = "scale.data")+ NoLegend()
write.table(data.frame(data.frame(pbmc.markers)), 
            'F:/Data1-2 cluster marker gene.xlsx',sep = '\t', quote = FALSE, row.names = FALSE)

FeaturePlot(biofilm.combined_6, features = c("b3517","b1493"), pt.size = 0.15)
DotPlot(object = biofilm.combined_6,features =c("b3517","b1493"))
FeaturePlot(biofilm.combined_6,features = c("nCount_RNA"), pt.size = 0.85)
RidgePlot(biofilm.combined_6, features = c("b3517", "b1493"), ncol = 1, slot = "data")
FeatureScatter(biofilm.combined_6, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", 
               group.by = "seurat_clusters", shuffle = TRUE, jitter = TRUE)
VlnPlot(biofilm.combined_6, features = c("b3517", "b1493"), pt.size = 0.5, ncol = 2)
VlnPlot(biofilm.combined_6, features = c("b3517", "b1493"), 
        stack = TRUE, fill.by = "ident", flip = TRUE) 
FeaturePlot(biofilm.combined_6,features = c("b3517", "b1493"), pt.size = 0.85,
           # ncol = 5, 
           # pt.size = 0.1, 
            cols = colorRampPalette(c( "#91bfdb", "#ffffbf", "#FC8D62"))(100))# 将 my_colors 赋给 celltype 这一列的颜色属性
FeaturePlot(biofilm.combined_6,features = c("nCount_RNA"), pt.size = 0.85)

# 从 FindAllMarkers 结果统计每个 cluster 的显著基因数
sig_markers <- pbmc.markers %>% 
  filter(p_val_adj < 0.05, avg_log2FC > 0.5)
count_per_cluster <- table(sig_markers$cluster)
barplot(count_per_cluster, xlab = "Cluster", ylab = "# Marker Genes", 
        col = rainbow(length(count_per_cluster)))

save.image(file = "E:/数据分析/WH-LHB-scRNA-20260602/脚本/20260603-R.R")
