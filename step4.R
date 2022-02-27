rm(list = ls()) 



load(file = "step3output.Rdata")


library(reshape2)


head(expsaj)



expl=expsaj
exp_L = melt(expl)

head(exp_L)
colnames(exp_L)=c('symbol','sample','value')
head(exp_L)

library(stringr)
group_list = Group
group_list
exp_L$group = rep(group_list,each=nrow(expsaj))
head(exp_L)

# ggplot2
library(ggplot2)
p = ggplot(exp_L,
           aes(x=sample,y=value,fill=group))+geom_boxplot()
print(p)

##boxplot
p=ggplot(exp_L,aes(x=sample,y=value,fill=group))+geom_boxplot()
p=p+stat_summary(fun.y="mean",geom="point",shape=23,size=3,fill="red")
p=p+theme_set(theme_set(theme_bw(base_size=20)))
p=p+theme(text=element_text(face='bold'),axis.text.x=element_text(angle=30,hjust=1),axis.title=element_blank())
print(p)

head(exp_L)
table(exp_L[,2])


library(limma) 
exp2 = normalizeBetweenArrays(exp_L)
p=ggplot(exp_L,aes(x=sample,y=value,fill=group))+geom_violin()
print(p)

p=ggplot(exp_L,aes(value,fill=group))+geom_histogram(bins = 200)+facet_wrap(~sample, nrow = 4)
print(p)

p=ggplot(exp_L,aes(value,col=group))+geom_density()+facet_wrap(~sample, nrow = 4)
print(p)

p=ggplot(exp_L,aes(value,col=group))+geom_density() 
print(p)

# 
head(expsaj)
colnames(expsaj) = paste(group_list,1:9,sep='')
head(expsaj)
# 定义nodePar
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.7, col = "blue")
# 聚类
hc=hclust(dist(t(expsaj)))
par(mar=c(5,5,5,10)) 
# 绘图
plot(as.dendrogram(hc), nodePar = nodePar,  horiz = TRUE)






library(limma)
# 做分组矩阵 
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(expsaj)
design  #得到的分组矩阵

contrast.matrix<-makeContrasts(paste0(c("control","Schistosomiasis"),collapse = "-"),levels = design)

##step1
fit <- lmFit(expsaj,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix) ##这一步很重要，大家可以自行看看效果
fit2 <- eBayes(fit2)  ## default no trend !!!
##eBayes() with trend=TRUE
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 
#write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
head(nrDEG)
save(nrDEG,file = "DEGoutput.Rdata")


save(Group,group_list,exp,nrDEG,file = "step4output.Rdata")

write.csv(nrDEG,file = "nrDEG.csv")
















