#########################
## Author: Xi Wang
## Date: 2018/03/01
## Title: Add genename to vep table
#########################
library(dplyr)
library(biomaRt)#get gene name

anndir <- "/data1/workspace/DCI/Kirsch/Chang-Lung.Lee/WES_2018/Results/VEP/vep_output/"
outdir <- "/data1/workspace/DCI/Kirsch/Chang-Lung.Lee/WES_2018/Analysis/Xi_VEP/"

files <- list.files(path=anndir, pattern = '.tsv') #VEP table genereted by Dadong
mouse <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

addgene <-function(file,mart=mouse){
    vepfiltfile <- paste0(anndir, file)
    tools::md5sum(vepfiltfile)
    vep <- read.table(vepfiltfile, sep=' ', header=FALSE, stringsAsFactors=FALSE)
    colnames(vep) <- c("chr", "start", "end", "ALT", 
                        "EnsGene", "EnsFeature", 
                        "Annotation", "Impact")
    genesym <- unique(vep$EnsGene)

    res <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"), 
             filters="ensembl_gene_id", values=genesym, mart=mart)

    vepgene <- left_join(vep, res,  by=c("EnsGene"="ensembl_gene_id"))
    #dim(vepgene) == dim(vep)
    colnames(vepgene)[9]='Gene'
    write.csv(vepgene,file=paste0(anndir,file),sep='\t',col.names=FALSE,row.names=FALSE)
}


for (i in 1:length(files)){
    addgene(files[i])
}



