#########################
## Author: Xi Wang
## Date: 2018/02/21
## Title: Find Trp53 
##     mutation in vep table
#########################
library(dplyr)
library(biomaRt)#get gene name

#trp53 location in NCBI: chr11 69580359..69591873

anndir <- "/data1/workspace/DCI/Kirsch/Chang-Lung.Lee/WES_2018/Results/VEP/vep_output/"
outdir <- "/data1/workspace/DCI/Kirsch/Chang-Lung.Lee/WES_2018/Analysis/Xi_VEP/"


files <- list.files(path=anndir, pattern = '.tsv') #VEP table genereted by Dadong
file2 <- '/home.local/jsg57/kirsch_project/group3_table.tsv' #VEP table denerated by Jeremy

get_trp53mut <-function(file,chrnum=11,trp_start = 69580359,trp_end = 69591873 ){
    vepfiltfile <- paste0(anndir, file)
    tools::md5sum(vepfiltfile)
    vep <- read.table(vepfiltfile, sep=' ', header=FALSE, stringsAsFactors=FALSE)
    colnames(vep) <- c("chr", "start", "end", "ALT", 
                        "EnsGene", "EnsFeature", 
                        "Annotation", "Impact")
    vep['sample']=substr(file,1,3)
    trp_mut <- vep %>% filter(chr==chrnum, (start>=trp_start & start<=trp_end) | (end>=trp_start & end<=trp_end))
    return(trp_mut)
}

trp_mut_df={}
for (i in 1:length(files)){
    trp_mut_df = rbind(trp_mut_df,get_trp53mut(files[i]))
}
dim(trp_mut_df)
#0 9

#reproduce the result using file2
get_trp53mut2 <-function(file,chrnum=11,trp_start = 69580359,trp_end = 69591873 ){
    vepfiltfile <- file
    tools::md5sum(vepfiltfile)
    vep <- read.table(vepfiltfile, sep=' ', header=FALSE, stringsAsFactors=FALSE)
    colnames(vep) <- c("sample","chr", "start", "end", "ALT", 
                        "EnsGene", "EnsFeature", 
                        "Annotation", "Impact")
    trp_mut <- vep %>% filter(chr==chrnum, (start>=trp_start & start<=trp_end) | (end>=trp_start & end<=trp_end))
    #trp_mut <- vep %>% filter(sample %in% c('S28','S29','S31','S32','S33'))
    return(trp_mut)
}
trp_mut_df2 <- get_trp53mut2(file2)
dim(trp_mut_df2)
#write.csv(trp_mut_df2,file=paste0(outdir,'trp53_mut.csv'))
sum(trp_mut_df2[,-1]!=trp_mut_df[,-9])
#0: the two vep generate the same results.


####################################################
#get the total number of gene mutated in each sample
####################################################
get_genecount <- function(file){
    vepfiltfile <- paste0(anndir, file)
    tools::md5sum(vepfiltfile)
    vep <- read.table(vepfiltfile, sep=' ', header=FALSE, stringsAsFactors=FALSE)
    colnames(vep) <- c("chr", "start", "end", "ALT", 
                        "EnsGene", "EnsFeature", 
                        "Annotation", "Impact")
    vep_gene <- vep %>% dplyr:: select(-EnsFeature,-Impact) %>% filter(EnsGene!='-')
    genesym <- unique(vep_gene$EnsGene)
    #length(genesym) 
    mouse <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
    res <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"), 
             filters="ensembl_gene_id", values=genesym, mart=mouse)
    #dim(res)
    #setdiff(genesym, res[,1])
    uvepgene <- left_join(vep_gene, res,  by=c("EnsGene"="ensembl_gene_id"))
    VEPAnno <- aggregate(external_gene_name ~ chr+start+end+ALT, 
               data=uvepgene, FUN=function(x) paste(unique(x), collapse="&"))
    colnames(VEPAnno)[5] <- "Gene" #rename 'external_gene_name'
    count <- data.frame(sample=substr(file,1,3),
               Ensgene_n = length(genesym),
               gene_n = dim(unique(res['external_gene_name']))[1],
               agg_gene_n = length(unique(VEPAnno$Gene)))
    return(count)
}

 

count={}
for (i in 1:length(files)){
    count = rbind(count,get_genecount(files[i]))
}
count





#  sample Ensgene_n gene_n agg_gene_n
#1    S28     10183  10167       9626
#2    S29      1979   1977       1426
#3    S31      1973   1972       1413
#4    S32      2299   2296       1619
#5    S33      2007   2002       1417
# Ensgene_n is the number of unique EnsGene ID 
# gene_n is the number of unique gene names after assigning the gene name from BioMart database to each EnsGene ID.
# agg_gene_n is the number of unique gene names after aggregating the same mutation to one gene name. (because same genetic position might be assigned to different gene names)

#saved in gene_count.txt under the same folder.










