rm(list = ls()) 
options(stringsAsFactors = F)
load(file = "step4output.Rdata")

load(file = "step3output.Rdata")

library(pheatmap)
choose_gene = head(rownames(nrDEG),)
choose_matrix = exp[choose_gene,]
choose_matrix = t(scale(t(choose_matrix)))
pheatmap(choose_matrix)
colnames(nrDEG)
plot(nrDEG$logFC,-log10(nrDEG$P.Value))



dat=expsaj
dat[1:4,1:4] 


# 2.top 1000 sd 热图---- 
cg=names(tail(sort(apply(dat,1,sd)),1000))
n=dat[cg,]
library(pheatmap)
annotation_col=data.frame(group=Group)#注解
rownames(annotation_col)=colnames(n) 
# 用标准化的数据画热图，两种方法的比较：https://mp.weixin.qq.com/s/jW59ujbmsKcZ2_CM5qRuAg，效果是一样的
## 1.使用热图参数
pheatmap(n,
         show_colnames =F,
         show_rownames = F,
         annotation_col=annotation_col,
         scale = "row",
         breaks = seq(-3,3,length.out = 100)
) 
h1
#breaks 参数解读在上面链接
dev.off()
## 2.自行标准化再画热图
n2 = t(scale(t(n)))
pheatmap(n2,
         show_colnames =F,
         show_rownames = F,
         annotation_col=annotation_col,
         breaks = seq(-3,3,length.out = 100)
)
dev.off()


deg=nrDEG

#为deg数据框添加几列
#1.加probe_id列，把行名变成一列
library(dplyr)
deg <- mutate(deg,probe_id=rownames(deg))
head(deg)




#3.加change列,标记上下调基因
logFC_t=1.5
P.Value_t = 0.05
k1 = (deg$P.Value < P.Value_t)&(deg$logFC < -logFC_t)
k2 = (deg$P.Value < P.Value_t)&(deg$logFC > logFC_t)
change = ifelse(k1,"down",ifelse(k2,"up","stable"))
deg <- mutate(deg,change)

#1.火山图-----
library(dplyr)
library(ggplot2)
dat  = deg
#ggplot是可以赋值的
p <- ggplot(data = dat, 
            aes(x = logFC, 
                y = -log10(P.Value))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  theme_bw()
p

#加label，设置想要加上的label用if(T)
if(T){
  #自选基因
  for_label <- dat%>% 
    filter(symbol %in% c("HADHA","LRRFIP1")) 
}
if(F){
  #p值最小的10个
  for_label <- dat %>% head(10)
}

if(F) {
  #p值最小的前3下调和前3上调
  x1 = dat %>% 
    filter(change == "up") %>% 
    head(5)
  x2 = dat %>% 
    filter(change == "down") %>% 
    head(5)
  for_label = rbind(x1,x2)
}
#实现在图上加label


volcano_plot <- p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = probe_id),
    data = for_label,
    color="black"
  )
volcano_plot
ggsave(plot = volcano_plot,filename = paste0(gse_number,"_volcano.png"))

write.csv(deg,file = "deg.csv")
