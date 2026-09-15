library(Seurat)
library(cowplot)
library(patchwork)
library(dplyr)
library(harmony)
Data1 <- read.table("D:/研究课题/LHB-gadAB-文章材料/0427-3.txt", sep="\t")
Data11 <- Data1[-1,-1]   
colnames(Data11) <- Data1[1,-1]
rownames(Data11) <- Data1[-1,1]
Data2 <- read.table("D:/研究课题/LHB-gadAB-文章材料/0427-4.txt", sep="\t")
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

DefaultAssay(ecoli24h_1_3) <- "RNA"
ecoli24h_1_3 <- ScaleData(ecoli24h_1_3, verbose = FALSE)
ecoli24h_1_3 <- RunPCA(ecoli24h_1_3, npcs = 30, verbose = FALSE)
DefaultAssay(ecoli24h_2_3) <- "RNA"
ecoli24h_2_3 <- ScaleData(ecoli24h_2_3, verbose = FALSE)
ecoli24h_2_3 <- RunPCA(ecoli24h_2_3, npcs = 30, verbose = FALSE)
obj.list <- list(ecoli24h_1_3, ecoli24h_2_3)
## ---------- 2. CCA 整合（主流程 A、B 共用） ----------
anchors.cca <- FindIntegrationAnchors(object.list = obj.list,reduction = "cca", dims = 1:20)
combined.cca <- IntegrateData(anchorset = anchors.cca, dims = 1:20)
DefaultAssay(combined.cca) <- "integrated"
combined.cca <- ScaleData(combined.cca, verbose = FALSE)
combined.cca <- RunPCA(combined.cca, npcs = 30, verbose = FALSE)

library(ggplot2)
library(patchwork)
# A. CCA + Harmony
objA <- RunHarmony(combined.cca, group.by.vars = "orig.ident")
objA <- RunUMAP(objA, reduction = "harmony", dims = 1:5, reduction.name = "umap")
objA <- FindNeighbors(objA, reduction = "harmony", dims = 1:5)
objA <- FindClusters(objA, resolution = 0.17)
# B. CCA 不加 Harmony
objB <- RunUMAP(combined.cca, reduction = "pca", dims = 1:5, reduction.name = "umap")
objB <- FindNeighbors(objB, reduction = "pca", dims = 1:5)
objB <- FindClusters(objB, resolution = 0.17)
## ---------- 5. 对比图 ----------
pA <-DimPlot(objA, reduction = "umap", group.by = "orig.ident") + ggtitle("A: CCA + Harmony (batch)")
pB <-DimPlot(objB, reduction = "umap", group.by = "orig.ident") + ggtitle("B: CCA only (batch)")

pA2 <-DimPlot(objA, reduction = "umap", label = TRUE) + ggtitle("A: CCA + Harmony (clusters)")
pB2 <- DimPlot(objB, reduction = "umap", label = TRUE) + ggtitle("B: CCA only (clusters)")

(pA | pB ) / (pA2 | pB2 )
PPA<-FeaturePlot(objA, features = c("b3517","b1493"), pt.size = 0.15)
PPB<-FeaturePlot(objB, features = c("b3517","b1493"), pt.size = 0.15)
PPA / PPB 

# 样本1独立分析
obj1 <- ecoli24h_1_3
obj1 <- ScaleData(obj1, verbose = FALSE)
obj1 <- RunPCA(obj1, npcs = 30, verbose = FALSE)
obj1 <- FindNeighbors(obj1, dims = 1:10)
obj1 <- FindClusters(obj1, resolution = 0.17)
obj1 <- RunUMAP(obj1, dims = 1:5)
P1a<-DimPlot(obj1, reduction = 'umap', label = TRUE, pt.size = 0.85)
P1<-FeaturePlot(obj1, features = c("b3517","b1493"))

obj2 <- ecoli24h_2_3
obj2 <- ScaleData(obj2, verbose = FALSE)
obj2 <- RunPCA(obj2, npcs = 30, verbose = FALSE)
obj2 <- FindNeighbors(obj2, dims = 1:10)
obj2 <- FindClusters(obj2, resolution = 0.17)
obj2 <- RunUMAP(obj2, dims = 1:5)
P2a<-DimPlot(obj2, reduction = 'umap', label = TRUE, pt.size = 0.85)
P2<-FeaturePlot(obj2, features = c("b3517","b1493"))
P1a | P2a
P1 / P2


