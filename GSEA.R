

#引用包
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(tidyverse)
#读取文件,并对输入文件进行整理
rt=read.table("TCGA_all_TPM.txt", header=T, sep="\t", check.names=F)
#以01A和11A分组，提取肿瘤数据
rt_T=rt%>% dplyr::select(str_which(colnames(.), "-01A"))
data=avereps(rt_T)


#高低表达组比较，得到logFC
gene="CFTR"              #基因名称
dataL=data[,data[gene,]<=median(data[gene,]),drop=F]
dataH=data[,data[gene,]>median(data[gene,]),drop=F]
meanL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=log2(meanH)-log2(meanL)
logFC=sort(logFC,decreasing=T)
genes=names(logFC)


#GO
#读入基因集文件
gmt=read.gmt("c5.go.v7.4.symbols.gmt")

#富集分析
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$pvalue<0.05,]
write.table(kkTab,file=paste0(gene,"_GSEA.result-GO.txt"),sep="\t",quote=F,row.names = F)

#输出富集的图形
termNum=5    #展示前5个通路
if(nrow(kkTab)>=termNum){
  showTerm=row.names(kkTab)[1:termNum]
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title="")
  gseaplot[[1]]=gseaplot[[1]]+theme(legend.position="right",legend.direction="vertical")
  pdf(file=paste0(gene,"_GSEA-GO.pdf"), width=10, height=8)
  print(gseaplot)
  dev.off()
}
if(nrow(kkTab)<termNum){
  showTerm=row.names(kkTab)[1:nrow(kkTab)]
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title="")
  gseaplot[[1]]=gseaplot[[1]]+theme(legend.position="right",legend.direction="vertical")
  pdf(file=paste0(gene,"_GSEA-GO.pdf"), width=10, height=8)
  print(gseaplot)
  dev.off()
}

#KEGG
#读入基因集文件
gmt=read.gmt("c2.cp.kegg.v7.4.symbols.gmt")

#富集分析
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$pvalue<0.05,]
write.table(kkTab,file=paste0(gene,"_GSEA.result-KEGG.txt"),sep="\t",quote=F,row.names = F)

#输出富集的图形
termNum=5    #展示前5个通路
if(nrow(kkTab)>=termNum){
  showTerm=row.names(kkTab)[1:termNum]
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title="")
  gseaplot[[1]]=gseaplot[[1]]+theme(legend.position="right",legend.direction="vertical")
  pdf(file=paste0(gene,"_GSEA-KEGG.pdf"), width=10, height=8)
  print(gseaplot)
  dev.off()
}
if(nrow(kkTab)<termNum){
  showTerm=row.names(kkTab)[1:nrow(kkTab)]
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title="")
  gseaplot[[1]]=gseaplot[[1]]+theme(legend.position="right",legend.direction="vertical")
  pdf(file=paste0(gene,"_GSEA-KEGG.pdf"), width=10, height=8)
  print(gseaplot)
  dev.off()
}
