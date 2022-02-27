rm(list = ls())  
load(file = "step1output.Rdata")
library(stringr)
# 1.Group-----

Group = pd$`infection:ch1`
Group=c(rep("control",times=4),
        rep("Schistosomiasis",times=6))




a = getGEO(gpl_number,destdir = ".")
b = a@dataTable@table
colnames(b)
ids2 = b[,c("ID","Symbol")]
colnames(ids2) = c("probe_id","symbol")
ids2 = ids2[ids2$symbol!="" & !str_detect(ids2$symbol,"///"),]

ids=ids2
save(exp,expsaj,Group,ids,gse_number,file = "step2output.Rdata")
write.csv(ids,"ids.csv")
