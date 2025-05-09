#加载R包
library(limma)
library(future.apply)
library(pbapply)
library(parallel)
library(doParallel)
#加载所需函数
{
merge_TCGA <- function(metadata, path, data.type, mRNA_expr_type="STAR", symbol=T, RNA_type=T){

  filenames <- file.path(path, metadata$file_id, metadata$file_name, 
                         fsep = .Platform$file.sep)
  if (data.type=='RNAseq') {
    message ('###############    正在进行提取，请稍后   ################')
    if(mRNA_expr_type=="STAR"){
      column=4
    }else if(mRNA_expr_type=="TPM"){
      column=7
    }else if(mRNA_expr_type=="FPKM"){
      column=8
    }else if(mRNA_expr_type=="FPKM_UQ"){
      column=9
    }
    plan(multisession)
    rnaMatrix <- do.call("cbind", future_lapply(filenames, function(fl)
      read.table(fl,skip=6,sep="\t")[,column]))
    ensembl <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V1
    gene_symbol <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V2
    type <- read.table(filenames[1],skip=6,sep="\t",stringsAsFactors = F)$V3
    index=grepl("^\\d+$",sapply(strsplit(ensembl, '.', fixed=TRUE), '[',2))
    rnaMatrix=rnaMatrix[index,]
    rownames(rnaMatrix) <- sapply(strsplit(ensembl[index], '.', fixed=TRUE), '[',1)
    gene_symbol=gene_symbol[index]
    type=type[index]
    colnames(rnaMatrix) <- metadata$sample
    nSamples = ncol(rnaMatrix)
    nGenes = nrow(rnaMatrix)
    if(RNA_type){
      rnaMatrix=data.frame(type,rnaMatrix,stringsAsFactors = F,check.names = F)
    }
    
    if(symbol){
      rnaMatrix=data.frame(gene_symbol,rnaMatrix,stringsAsFactors = F,check.names = F)
    }
    message (paste('Number of samples: ', nSamples, '\n', sep=''),
             paste('Number of genes: ', nGenes, '\n', sep=''))
    #返回最后的基因表达矩阵
    return (rnaMatrix)
    
  }else if (data.type=='miRNAs') { 
    message ('############### Merging miRNAs data ###############\n')
    mirMatrix <- future_lapply(filenames, function(fl) filtermir(fl))
    mirs <- mirbase$V1
    mirMatrix <- do.call('cbind', lapply(mirMatrix, 
                                         function(expr) expr[mirs]))
    rownames(mirMatrix) <- mirbase$V2
    colnames(mirMatrix) <- metadata$sample
    mirMatrix[is.na(mirMatrix)] <- 0
    nSamples = ncol(mirMatrix)
    nGenes = nrow(mirMatrix)
    message (paste('Number of samples: ', nSamples, '\n', sep=''),
             paste('Number of miRNAs: ', nGenes, '\n', sep=''))
    return (mirMatrix)
  }else{ 
    stop('data type error!')
  }
}

filtermir <- function(fl) {
  expr <- read.table(fl, header=TRUE, stringsAsFactors = FALSE)
  expr <- expr[startsWith(expr$miRNA_region, "mature"),]
  expr <- aggregate(expr$read_count, list(expr$miRNA_region), sum)
  mirs <- sapply(strsplit(expr$Group.1, ',', fixed=TRUE),'[',2)
  expr <- expr[,-1]
  names(expr) <- mirs
  return(expr)
}

FilterDuplicate <- function(metadata) {
  filter <- which(duplicated(metadata[,'sample']))
  if (length(filter) != 0) {
    metadata <- metadata[-filter,]
  }
  message (paste('Removed', length(filter), 'samples', sep=' '))
  return (metadata)
}
FilterSampleType <- function(metadata) {
  filter <- which(! metadata$sample_type %in% 
                    c('PrimaryTumor', 'SolidTissueNormal'))
  if (length(filter) != 0) {
    metadata <- metadata[-filter,]
  }
  message (paste('Removed', length(filter), 'samples', sep=' '))
  return (metadata)
}
metaMatrix.RNA=read.table("sheet.tsv",sep="\t",header=T)
names(metaMatrix.RNA)=gsub("sample_id","sample",gsub("\\.","_",tolower(names(metaMatrix.RNA))))
metaMatrix.RNA$sample_type=gsub(" ","",metaMatrix.RNA$sample_type)
metaMatrix.RNA <- FilterDuplicate(metaMatrix.RNA)
metaMatrix.RNA <- FilterSampleType(metaMatrix.RNA)
}

#TPM
RNA_TPM=merge_TCGA(metadata=metaMatrix.RNA, 
                           path="RNAseq", 
                           data.type="RNAseq",
                           mRNA_expr_type="TPM",
                           symbol = T,
                           RNA_type=T
)
#提取所有表达数据
RNA_TPM_all=as.matrix(RNA_TPM)
rownames(RNA_TPM_all)=RNA_TPM_all[,1]
RNA_TPM_all_exp=RNA_TPM_all[,3:ncol(RNA_TPM_all)]
RNA_TPM_all_dimnames=list(rownames(RNA_TPM_all_exp),colnames(RNA_TPM_all_exp))
RNA_TPM_all_data=matrix(as.numeric(as.matrix(RNA_TPM_all_exp)),nrow=nrow(RNA_TPM_all_exp),dimnames=RNA_TPM_all_dimnames)
RNA_TPM_all_data=avereps(RNA_TPM_all_data)
write.table(file="TCGA_all_TPM.txt",RNA_TPM_all_data,sep="\t",quote=F)
#提取lncRNA
RNA_TPM_lnc=RNA_TPM[RNA_TPM$type=="lncRNA",]
RNA_TPM_lnc=as.matrix(RNA_TPM_lnc)
rownames(RNA_TPM_lnc)=RNA_TPM_lnc[,1]
RNA_TPM_lnc_exp=RNA_TPM_lnc[,3:ncol(RNA_TPM_lnc)]
RNA_TPM_lnc_dimnames=list(rownames(RNA_TPM_lnc_exp),colnames(RNA_TPM_lnc_exp))
RNA_TPM_lnc_data=matrix(as.numeric(as.matrix(RNA_TPM_lnc_exp)),nrow=nrow(RNA_TPM_lnc_exp),dimnames=RNA_TPM_lnc_dimnames)
RNA_TPM_lnc_data=avereps(RNA_TPM_lnc_data)
write.table(file="TCGA_lnc_TPM.txt",RNA_TPM_lnc_data,sep="\t",quote=F)

#count
RNA_STAR_Counts=merge_TCGA(metadata=metaMatrix.RNA, 
                           path="RNAseq", 
                           data.type="RNAseq",
                           mRNA_expr_type="STAR",
                           symbol = T,
                           RNA_type=T
)
#提取所有
RNA_STAR_Counts_all=as.matrix(RNA_STAR_Counts)
rownames(RNA_STAR_Counts_all)=RNA_STAR_Counts_all[,1]
RNA_STAR_Counts_all_exp=RNA_STAR_Counts_all[,3:ncol(RNA_STAR_Counts_all)]
RNA_STAR_Counts_all_dimnames=list(rownames(RNA_STAR_Counts_all_exp),colnames(RNA_STAR_Counts_all_exp))
RNA_STAR_Counts_all_data=matrix(as.numeric(as.matrix(RNA_STAR_Counts_all_exp)),nrow=nrow(RNA_STAR_Counts_all_exp),dimnames=RNA_STAR_Counts_all_dimnames)
RNA_STAR_Counts_all_data=avereps(RNA_STAR_Counts_all_data)
write.table(file="TCGA_all_counts.txt",RNA_STAR_Counts_all_data,sep="\t",quote=F)

#提取lncRNA
RNA_STAR_Counts_lnc=RNA_STAR_Counts[RNA_STAR_Counts$type=="lncRNA",]
RNA_STAR_Counts_lnc=as.matrix(RNA_STAR_Counts_lnc)
rownames(RNA_STAR_Counts_lnc)=RNA_STAR_Counts_lnc[,1]
RNA_STAR_Counts_lnc_exp=RNA_STAR_Counts_lnc[,3:ncol(RNA_STAR_Counts_lnc)]
RNA_STAR_Counts_lnc_dimnames=list(rownames(RNA_STAR_Counts_lnc_exp),colnames(RNA_STAR_Counts_lnc_exp))
RNA_STAR_Counts_lnc_data=matrix(as.numeric(as.matrix(RNA_STAR_Counts_lnc_exp)),nrow=nrow(RNA_STAR_Counts_lnc_exp),dimnames=RNA_STAR_Counts_lnc_dimnames)
RNA_STAR_Counts_lnc_data=avereps(RNA_STAR_Counts_lnc_data)
write.table(file="TCGA_lnc_counts.txt",RNA_STAR_Counts_lnc_data,sep="\t",quote=F)

#FPKM
RNA_FPKM=merge_TCGA(metadata=metaMatrix.RNA, 
                    path="RNAseq", 
                    data.type="RNAseq",
                    mRNA_expr_type="FPKM",
                    symbol = T,
                    RNA_type=T
)
#提取所有
RNA_FPKM_all=as.matrix(RNA_FPKM)
rownames(RNA_FPKM_all)=RNA_FPKM_all[,1]
RNA_FPKM_all_exp=RNA_FPKM_all[,3:ncol(RNA_FPKM_all)]
RNA_FPKM_all_dimnames=list(rownames(RNA_FPKM_all_exp),colnames(RNA_FPKM_all_exp))
RNA_FPKM_all_data=matrix(as.numeric(as.matrix(RNA_FPKM_all_exp)),nrow=nrow(RNA_FPKM_all_exp),dimnames=RNA_FPKM_all_dimnames)
RNA_FPKM_all_data=avereps(RNA_FPKM_all_data)
write.table(file="TCGA_all_FPKM.txt",RNA_FPKM_all_data,sep="\t",quote=F)
#提取lncRNA
RNA_FPKM_lnc=RNA_FPKM[RNA_FPKM$type=="lncRNA",]
RNA_FPKM_lnc=as.matrix(RNA_FPKM_lnc)
rownames(RNA_FPKM_lnc)=RNA_FPKM_lnc[,1]
RNA_FPKM_lnc_exp=RNA_FPKM_lnc[,3:ncol(RNA_FPKM_lnc)]
RNA_FPKM_lnc_dimnames=list(rownames(RNA_FPKM_lnc_exp),colnames(RNA_FPKM_lnc_exp))
RNA_FPKM_lnc_data=matrix(as.numeric(as.matrix(RNA_FPKM_lnc_exp)),nrow=nrow(RNA_FPKM_lnc_exp),dimnames=RNA_FPKM_lnc_dimnames)
RNA_FPKM_lnc_data=avereps(RNA_FPKM_lnc_data)
write.table(file="TCGA_lnc_FPKM.txt",RNA_FPKM_lnc_data,sep="\t",quote=F)