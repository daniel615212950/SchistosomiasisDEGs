
#################################
##clusterProfiler包进行富集分析##
#################################
#相关依赖R包安装
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
library(clusterProfiler)
BiocManager::install("topGO")
library("topGO")
BiocManager::install("clusterProfiler")
library(clusterProfiler)
BiocManager::install("Rgraphviz")
library("Rgraphviz")
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(stringr)
library(enrichplot)


# Go富集
#go富集分析--耗费时间灰常长，很正常

rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load(file = "step4output.Rdata")


library(dplyr)
deg = nrDEG
deg <- mutate(deg,symbol = rownames(deg))
head(deg)


logFC_t=1.5
P.Value_t = 0.05
k1 = (deg$P.Value < P.Value_t)&(deg$logFC < -logFC_t)
k2 = (deg$P.Value < P.Value_t)&(deg$logFC > logFC_t)
change = ifelse(k1,"down",ifelse(k2,"up","stable"))
deg <- mutate(deg,change)

head(deg)
table(deg$change)

###3.加ENTREZID列，后面富集分析要用
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
s2e <- bitr(unique(deg$symbol), fromType = "SYMBOL",  #ID转换核心函数bitr
            toType = c( "ENTREZID"),
            OrgDb = org.Hs.eg.db)
head(s2e)
head(deg)
deg <- inner_join(deg,s2e,by=c("symbol"="SYMBOL"))

head(deg)
table(deg$change)
save(deg,file = "deg.Rdata")


rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)

load(file = "deg.Rdata")
write.csv(deg,file = "degGO.csv")


##############################################################
rm(list = ls()) 
load(file = "deg.Rdata")
library(clusterProfiler)
#输入数据
gene_up= deg[deg$change == 'up','ENTREZID'] 
gene_down=deg[deg$change == 'down','ENTREZID'] 
gene_diff=c(gene_up,gene_down)
gene_all = deg[,'ENTREZID']

#**GO分析三大块总体分析**
#*
#*
#*

go <- enrichGO(gene_diff, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                 qvalueCutoff = 0.2,keyType = 'ENTREZID')

head(go)

write.csv(go,file = "go.csv")

dim(go)
dim(go[go$ONTOLOGY=='BP',])
dim(go[go$ONTOLOGY=='CC',])
dim(go[go$ONTOLOGY=='MF',])
barplot(go,showCategory=20,drop=T)
dotplot(go,showCategory=20)

dotplot(go, split = "ONTOLOGY", font.size = 10, 
        showCategory = 5) + facet_grid(ONTOLOGY ~ ., scale = "free") + 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 45))

#**GO分析分组分析**
#*
#*
#*


#细胞组分
ego_CC <- enrichGO(gene = gene_diff,
                   OrgDb= org.Hs.eg.db,
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.01,
                   readable = TRUE)
#生物过程
ego_BP <- enrichGO(gene = gene_diff,
                   OrgDb= org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.01,
                   readable = TRUE)
#分子功能：
ego_MF <- enrichGO(gene = gene_diff,
                   OrgDb= org.Hs.eg.db,
                   ont = "MF",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.01,
                   readable = TRUE)
save(ego_CC,ego_BP,ego_MF,file = "GO_saj.Rdata")
rm(list = ls())  ## 魔幻操作，一键清空~   
load(file = "GO_saj.Rdata")
#作图
#第一种，条带图，按p从小到大排的
barplot(ego_CC, showCategory=20,title="EnrichmentGO_CC")
barplot(ego_BP, showCategory=20,title="EnrichmentGO_CC")
#如果运行了没出图，就dev.new()
#第二种，点图，按富集数从大到小的
dotplot(ego_CC,title="EnrichmentGO_BP_dot")
#细胞通路图
go.BP <- enrichGO(gene_diff, OrgDb = org.Hs.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
plotGOgraph(go.BP, firstSigNodes = 10, useInfo = "all",sigForAll = TRUE,useFullNames = TRUE)

#展示top通路的共同基因，要放大看。
#Gene-Concept Network
cnetplot(go,categorySize="pvalue", foldChange=gene_diff,colorEdge = TRUE)
cnetplot(go, showCategory = 5,foldChange=gene_diff, circular = TRUE, colorEdge = TRUE)

dev.off()


#KEGG富集分析
#（1）富集分析整体
kegg <- enrichKEGG(gene_diff, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', 
                   minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
head(kegg)
dim(kegg)
dotplot(kegg, showCategory=30)

gene_up= deg[deg$change == 'up','ENTREZID'] 
gene_down=deg[deg$change == 'down','ENTREZID'] 
gene_diff=c(gene_up,gene_down)
gene_all = deg[,'ENTREZID']


#通过ENSEMBL，获取ENTREZID，GENENAME


enrich.kegg <- enrichKEGG(gene =deg$ENTREZID,
                          organism ="hsa",
                          keyType = "kegg",
                          pvalueCutoff = 1,
                          pAdjustMethod = "BH",
                          minGSSize = 10,
                          maxGSSize = 500,
                          qvalueCutoff = 1,
                          use_internal_data =FALSE)
dim(enrich.kegg)

sig.kegg<-filter(enrich.kegg,pvalue<.05,qvalue<0.2)
dim(sig.kegg)
kegg.bar<-barplot(sig.kegg,showCategory=20,color = "pvalue")
kegg.dot<-dotplot(sig.kegg,showCategory=20,color = "pvalue")
db<-plot_grid(kegg.bar,kegg.dot,ncol=2)
db
ggsave("kegg_bar_dot1.tiff",db,width=14,height=6)


#########################
###GSEA富集分析##########
library(dplyr)
geneList <- select(degenes, Entrez.ID, Fold.Change)







#########################
###pathview包##########
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pathview")

rm(list = ls())  
load(file = "step1output.Rdata")
load(file = "step2output.Rdata")
length(unique(ids$symbol))
tail(sort(table(ids$symbol)))
table(rownames(expsaj) %in% ids$probe_id)
expsaj = expsaj[rownames(expsaj) %in% ids$probe_id,]
ids=ids[match(rownames(expsaj),ids$probe_id),]
tmp = by(expsaj,
         ids$symbol,
         function(x) rownames(x)[which.max(rowMeans(x))])

probes = as.character(tmp)
expsaj = expsaj[rownames(expsaj) %in% probes,]
dim(expsaj)

rownames(expsaj)=ids[match(rownames(expsaj),ids$probe_id),2]

class(expsaj)

pathview(gene.data = NULL, cpd.data = NULL, pathway.id,
         species = "hsa", kegg.dir = ".", cpd.idtype = "kegg", gene.idtype =
           "entrez", gene.annotpkg = NULL, min.nnodes = 3, kegg.native = TRUE,
         map.null = TRUE, expand.node = FALSE, split.group = FALSE, map.symbol =
           TRUE, map.cpdname = TRUE, node.sum = "sum", discrete=list(gene=FALSE,
                                                                     cpd=FALSE), limit = list(gene = 1, cpd = 1), bins = list(gene = 10, cpd
                                                                                                                              = 10), both.dirs = list(gene = T, cpd = T), trans.fun = list(gene =
                                                                                                                                                                                             NULL, cpd = NULL), low = list(gene = "green", cpd = "blue"), mid =
           list(gene = "gray", cpd = "gray"), high = list(gene = "red", cpd =
                                                            "yellow"), na.col = "transparent", ...)


