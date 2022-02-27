
rm(list = ls())  
load(file = "step1output.Rdata")
load(file = "step2output.Rdata")


Group=c(rep("control",times=4),
        rep("Schistosomiasis",times=6))

# 1.PCA ---
dat=t(expsaj)
class(dat)
dat=as.data.frame(dat)
library(FactoMineR)
library(factoextra) 
dat.pca <- PCA(dat, graph =TRUE)
pca_plot <- fviz_pca_ind(dat.pca,
                         geom.ind = "point", # show points only (nbut not "text")
                         col.ind = Group, # color by groups
                         palette = c("#00AFBB", "#E7B800","#E7B000"),
                         addEllipses = TRUE, # Concentration ellipses
                         legend.title = "Groups"
)
pca_plot
ggsave(plot = pca_plot,filename = paste0(gse_number,"_PCA.png"))
save(pca_plot,exp,dat,pd,file = "pca_plot.Rdata")

dev.off()

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

save(exp,expsaj,Group,ids,gse_number,file = "step3output.Rdata")

