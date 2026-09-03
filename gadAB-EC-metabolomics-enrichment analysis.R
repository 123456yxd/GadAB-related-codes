install.packages("KEGGREST")
BiocManager::install('KEGGREST')
library(KEGGREST)
library(clusterProfiler)

library(clusterProfiler)
library(dplyr)
eco_pathways <- keggList("pathway", "eco")
get_pathway_compounds <- function(pathway_id) {
  pathway_info <- keggGet(pathway_id)
  compounds <- pathway_info[[1]]$COMPOUND
  if (!is.null(compounds)) {
    return(data.frame(
      pathway = pathway_id,
      compound = names(compounds),
      stringsAsFactors = FALSE
    ))
  } else {
    return(NULL)
  }
}

all_compounds_list <- lapply(names(eco_pathways), get_pathway_compounds)
kegg_eco_compound_pathway <- do.call(rbind, all_compounds_list)
head(kegg_eco_compound_pathway)


sig_metabolites<-read.table("E:/数据分析/ZZ7-LHB-mix-20260701/代谢组/cpd-ID-down-forkegg.txt")
all_metabolites<-read.table("E:/数据分析/ZZ7-LHB-mix-20260701/代谢组/cpd_ID-all.txt")#: 字符向量，包含所有检测到的代谢物ID (作为背景集)

sig_metabolites<-sig_metabolites[,1]
all_metabolites<-all_metabolites[,1]
enrich_result <- enricher(
  gene = sig_metabolites,          
  universe = all_metabolites,      
  TERM2GENE = kegg_eco_compound_pathway[, c("pathway", "compound")], # 通路-代谢物对应表
  pvalueCutoff = 0.75,             
  qvalueCutoff = 0.2,              
  minGSSize = 3                    # 通路中最少包含的代谢物数
)

# --- 查看与可视化结果 ---
head(enrich_result@result)
enrich_result@result
barplot(enrich_result,showCategory = 20,color = "pvalue")#, by = "ID")#, showCategory = 10, title = "KEGG Pathway Enrichment")
dotplot(enrich_result, showCategory = 10, title = "KEGG Pathway Enrichment")
write.csv(enrich_result,"E:/数据分析/ZZ7-LHB-mix-20260701/代谢组/代谢组做kegg go/cpd-ID-up-kegg-result.csv")#GO-result-ALL.csv")


Tn2<-read.table("E:/数据分析/ZZ7-LHB-mix-20260701/up-go.txt")
Tn5<-Tn2[,1]
ego <- enrichGO(
  gene          = Tn5,
  keyType = "ENTREZID",
  OrgDb         = org.EcK12.eg.db,
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.5,
  readable      = TRUE)
dotplot(ego)#,showCategory = 50)
write.csv(ego,"E:/数据分析/ZZ7-LHB-mix-20260701/KEGG GO/up-GO-gc.csv")#GO-result-ALL.csv")

Tn2<-read.table("E:/数据分析/ZZ7-LHB-mix-20260701/down-kegg.txt")
Tn5<-Tn2[,1]
kegg <- enrichKEGG(
  gene = Tn5,# up_genes
  organism = 'eco', 
  pAdjustMethod = 'fdr', 
  pvalueCutoff = 0.05, #0.2  #0.4  #0.3  #0.6
  qvalueCutoff = 0.5, 
)
dotplot(kegg)
write.csv(kegg,"E:/数据分析/ZZ7-LHB-mix-20260701/KEGG GO/dw-KEGG-GC.csv")

