#富集分析

#富集分析准备工作：

##首先对差异表达矩阵nrDEG，进行加工
###1.把行名变成SYMBOL列
rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load(file = "DEGoutput.Rdata")
library(dplyr)
deg = nrDEG
deg <- mutate(deg,symbol = rownames(deg))
head(deg)

###2.加change列：上调或下调，火山图要用


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

save(deg,file = "enrich_input.Rdata")

#####################
######富集分析#######
#####################

rm(list = ls()) 
options(stringsAsFactors = F)
load(file = 'enrich_input.Rdata')

## 1.KEGG pathway analysis
#上调、下调、差异、所有基因

#clusterProfiler作kegg富集分析：
library(clusterProfiler)
gene_up= deg[deg$change == 'up','ENTREZID'] 
gene_down=deg[deg$change == 'down','ENTREZID'] 
gene_diff=c(gene_up,gene_down)
gene_all = deg[,'ENTREZID']
kk.up <- enrichKEGG(gene         = gene_up,
                    organism     = 'hsa',
                    universe     = gene_all,
                    pvalueCutoff = 0.9,
                    qvalueCutoff =0.9)
head(kk.up)[,1:6]
dim(kk.up)
kk.down <- enrichKEGG(gene         =  gene_down,
                      organism     = 'hsa',
                      universe     = gene_all,
                      pvalueCutoff = 0.9,
                      qvalueCutoff =0.9)
head(kk.down)[,1:6]
dim(kk.down)
kk.diff <- enrichKEGG(gene         = gene_diff,
                      organism     = 'hsa',
                      pvalueCutoff = 0.05)
head(kk.diff)[,1:6]

class(kk.diff)
#提取出数据框
kegg_diff_dt <- kk.diff@result

#根据pvalue来选,用于可视化
down_kegg <- kk.down@result %>%
  filter(pvalue<0.05) %>%
  mutate(group=-1)

up_kegg <- kk.up@result %>%
  filter(pvalue<0.05) %>%
  mutate(group=1)

#可视化
kegg_plot <- function(up_kegg,down_kegg){
  dat=rbind(up_kegg,down_kegg)
  colnames(dat)
  dat$pvalue = -log10(dat$pvalue)
  dat$pvalue=dat$pvalue*dat$group 
  
  dat=dat[order(dat$pvalue,decreasing = F),]
  
  g_kegg<- ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 
    geom_bar(stat="identity") + 
    scale_fill_gradient(low="blue",high="red",guide = FALSE) + 
    scale_x_discrete(name ="Pathway names") +
    scale_y_continuous(name ="log10P-value") +
    coord_flip() + theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
    ggtitle("Pathway Enrichment") 
}

g_kegg <- kegg_plot(up_kegg,down_kegg)
g_kegg

ggsave(g_kegg,filename = 'kegg_up_down.png')


#gsea作kegg富集分析：

data(geneList, package="DOSE")
head(geneList)
length(geneList)
names(geneList)
boxplot(geneList)
boxplot(deg$logFC)

geneList=deg$logFC
names(geneList)=deg$ENTREZID
geneList=sort(geneList,decreasing = T)

kk_gse <- gseKEGG(geneList     = geneList,
                  organism     = 'hsa',
                  nPerm        = 1000,
                  minGSSize    = 120,
                  pvalueCutoff = 0.9,
                  verbose      = FALSE)
head(kk_gse)[,1:6]
gseaplot(kk_gse, geneSetID = rownames(kk_gse[1,]))

down_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore < 0,];down_kegg$group=-1
up_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore > 0,];up_kegg$group=1

gse_kegg=kegg_plot(up_kegg,down_kegg)
print(gse_kegg)
ggsave(gse_kegg,filename ='kegg_up_down_gsea.png')


### 2.GO database analysis 
load(file = "step1output.Rdata")
#go富集分析
library(clusterProfiler)
#输入数据
gene_up= deg[deg$change == 'up','ENTREZID'] 
gene_down=deg[deg$change == 'down','ENTREZID'] 
gene_diff=c(gene_up,gene_down)
head(deg)

#**GO分析三大块**
# 1.GO 富集分析-----

#(1)输入数据
gene_up = deg[deg$change == 'up','ENTREZID'] 
gene_down = deg[deg$change == 'down','ENTREZID'] 
gene_diff = c(gene_up,gene_down)
gene_all = deg[,'ENTREZID']

#(2)富集
#以下步骤耗时很长，设置了存在即跳过
if(!file.exists(paste0(gse_number,"_GO.Rdata"))){
  ego <- enrichGO(gene = gene_diff,
                  OrgDb= org.Hs.eg.db,
                  ont = "ALL",
                  readable = TRUE)
  #ont参数：One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
  save(ego,file = paste0(gse_number,"_GO.Rdata"))
}
load(paste0(gse_number,"_GO.Rdata"))

#(3)可视化
#条带图
barplot(ego)
#气泡图
dotplot(ego)
#个性化设置
dotplot(ego, split = "ONTOLOGY", font.size = 10, 
        showCategory = 5) + facet_grid(ONTOLOGY ~ ., scale = "free") + 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 45))

#geneList 用于设置下面图的颜色
geneList = deg$logFC
names(geneList)=deg$ENTREZID
geneList = sort(geneList,decreasing = T)

#(3)展示top通路的共同基因，要放大看。
#Gene-Concept Network
cnetplot(ego,categorySize="pvalue", foldChange=geneList,colorEdge = TRUE)
cnetplot(ego, showCategory = 3,foldChange=geneList, circular = TRUE, colorEdge = TRUE)
#Enrichment Map,这个函数最近更新过，版本不同代码会不同
Biobase::package.version("enrichplot")

if(T){
  emapplot(pairwise_termsim(ego)) #新版本
}else{
  emapplot(ego)#老版本
}

#(4)展示通路关系 https://zhuanlan.zhihu.com/p/99789859
#goplot(ego)
#(5)Heatmap-like functional classification
heatplot(ego,foldChange = geneList,showCategory = 8)

