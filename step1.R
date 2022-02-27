getwd()

setwd()

rm(list = ls())
options(stringsAsFactors = F)
library(GEOquery)
gse_number = "GSE61376"
eSet <- getGEO(gse_number, 
               destdir = '.', 
               getGPL = F)
class(eSet)
length(eSet)
eSet = eSet[[1]]
#(1)提取表达矩阵exp
exp <- exprs(eSet)
exp[1:4,1:4]#查看数据，判断需不需要取log2
exp <- log2(exp+1)
exp[1:4,1:4]
boxplot(exp)#箱线图是用来查看数据是否有异常值
#(2)提取临床信息
pd <- pData(eSet)
#(3)调整pd的行名顺序与exp列名完全一致
p = identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]
#(4)提取芯片平台编号
gpl_number <- eSet@annotation
expsaj=exp[,1:10]

save(gse_number,pd,exp,expsaj,gpl_number,file = "step1output.Rdata")



write.csv(pd,"pd.csv")
write.csv(exp,"exp.csv")
write.csv(expsaj,"expsaj.csv")
