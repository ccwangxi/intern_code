#########################
## Author: XXXXXXXXX
## Date: 2018/02/21
## Title: Find AAAAA gene 
##     mutation in vep table
#########################
library(dplyr)
library(biomaRt)#get gene name

#tAAAA location in NCBI: chrn, start,end

anndir <- ###########
outdir <- ###########


files <- list.files(path=anndir, pattern = '.tsv') #VEP table genereted by person 1
file2 <- ########## #VEP table denerated by person 2

get_AAAmut <-function(file,chrnum=chrn,AA_start = start,AA_end = end ){
    vepfiltfile <- paste0(anndir, file)
    tools::md5sum(vepfiltfile)
    vep <- read.table(vepfiltfile, sep=' ', header=FALSE, stringsAsFactors=FALSE)
    colnames(vep) <- c("chr", "start", "end", "ALT", 
                        "EnsGene", "EnsFeature", 
                        "Annotation", "Impact")
    vep['sample']=substr(file,1,3)
    AA_mut <- vep %>% filter(chr==chrnum, (start>=AA_start & start<=AA_end) | (end>=AA_start & end<=AA_end))
    return(AA_mut)
}

AA_mut_df={}
for (i in 1:length(files)){
    AA_mut_df = rbind(AA_mut_df,get_AAAmut(files[i]))
}
dim(AA_mut_df)
#0 9

#reproduce the result using file2
get_AAAmut2 <-function(file,chrnum=chrn,AA_start = start,AA_end = end){
    vepfiltfile <- file
    tools::md5sum(vepfiltfile)
    vep <- read.table(vepfiltfile, sep=' ', header=FALSE, stringsAsFactors=FALSE)
    colnames(vep) <- c("sample","chr", "start", "end", "ALT", 
                        "EnsGene", "EnsFeature", 
                        "Annotation", "Impact")
    AA_mut <- vep %>% filter(chr==chrnum, (start>=AA_start & start<=AA_end) | (end>=AA_start & end<=AA_end))
    #AA_mut <- vep %>% filter(sample %in% XXXXXXXXXX)
    return(AA_mut)
}
AA_mut_df2 <- get_AAAmut2(file2)
dim(AA_mut_df2)
#write.csv(AA_mut_df2,file=paste0(outdir,'AAA_mut.csv'))
sum(AA_mut_df2[,-1]!=AA_mut_df[,-9])
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


# Ensgene_n is the number of unique EnsGene ID 
# gene_n is the number of unique gene names after assigning the gene name from BioMart database to each EnsGene ID.
# agg_gene_n is the number of unique gene names after aggregating the same mutation to one gene name. 
                         #(because same genetic position might be assigned to different gene names)










