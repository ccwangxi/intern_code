#########################
## Author: XXXXXX
## Date: 2018/03/14
## Title: Tables of counts of variants 
## for each experimental group/subgroup (excluding KRAS unirradiated)
## Recurring genelist for each subgroup
#########################
library(dplyr)
library(biomaRt)#get gene name
library(readxl)

outdir1 <- ############
count_outdir <- ###########

######WES_2018 (Sarcoma & lymphoma)
anndir1 <- #############
######Sarcoma
anndir2 <- #############
######lymphoma group3
anndir3 <- #############

## Select terms based on preferred VEP terms 
CLannlevs1 <- ###########
#############################    VEP    ########################################

get_veppass <-function(anndir,file){
    vepfiltfile <- paste0(anndir, file)
    vep <- read.table(vepfiltfile, sep=' ', header=FALSE, stringsAsFactors=FALSE)
    colnames(vep) <- c("seqnames", "start", "end", "ALT", 
                        "EnsGene", "EnsFeature", 
                        "Annotation", "Impact")
    vep['id']=strsplit(file,'\\-')[[1]][1]
    return(vep)
}
get_vep <- function(anndir){
    files <- list.files(path=anndir, pattern = '.tsv')
    veppass = {}
    for (i in 1:length(files)){
        veppass = rbind(veppass,get_veppass(anndir,files[i]))
    }
    return(veppass)
}


get_anno <- function(veppass, CLannlevs = CLannlevs1){
    ## Drop rows differing by transcript ID
    uvep <- veppass %>% 
                dplyr::select(-EnsFeature,-Impact) %>% 
                filter(EnsGene!='-') %>%       #Filter out those without Gene ID. 
                unique

    ## Reorder mutations by chr and pos
    uvep$chrn <- case_when(uvep$seqnames %in% c("M","MT")~22,
                          uvep$seqnames=="Y"~21,
                          uvep$seqnames=="X"~20,
                          TRUE~as.numeric(uvep$seqnames))
    uvep <- arrange(uvep, id, chrn, start, end)

    ## Expand comma-separated terms
    multi <- grep(",", uvep$Annotation, fixed=T)
    #if Annotation has ",", then there is more than one annotation. grep search for "," pattern.(if pattern shows, =1;otherwise,=0)
    uvepSingle <- uvep[-multi,] # With Single Annotation
    uvepMulti <- uvep[multi,] # With Multiple Annotation
    
    splitAnn <- strsplit(uvepMulti$Annotation, split=',') 
    uvepMultiExpanded <- {}
    for(i in 1:nrow(uvepMulti)){  # Abuse merge to duplicate uvep rows
      uvepMultiExpanded <- rbind.data.frame(uvepMultiExpanded, 
            merge(x=uvepMulti[i,-6], y=splitAnn[[i]]), 
            stringsAsFactors=FALSE)
    }
    colnames(uvepMultiExpanded)[8] <- "Annotation"
    uvepExp <- rbind.data.frame(uvepSingle, uvepMultiExpanded,
            stringsAsFactors=FALSE)
    uvepExp <- unique(uvepExp)

    uvep <- filter(uvepExp, Annotation %in% CLannlevs)

    ## Select most deleterious annotation within Gene
    uvep$AnnoN <- as.numeric(factor(uvep$Annotation, 
        levels=CLannlevs, ordered=TRUE))
    
    uvep <- data.frame(uvep %>% group_by(id, seqnames, start, end, ALT, EnsGene) %>% 
                  filter(AnnoN == min(AnnoN)) %>%
                  unique, stringsAsFactors=FALSE)
    return (uvep)
}


add_genename <- function(uvep,name, outdir=outdir1, CLannlevs = CLannlevs1){
    ## Get Gene names
    genesym <- unique(uvep$EnsGene)
    mouse <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
    res <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"), 
                 filters="ensembl_gene_id", values=genesym, mart=mouse)

    uvepgene <- inner_join(uvep, res,  by=c("EnsGene"="ensembl_gene_id"))  #remove those mutations without a gene name

    ## Aggregrate repeated rows (Gene and annotation)
    VEPAnno <- aggregate(cbind(external_gene_name, Annotation) ~ id+chrn+seqnames+start+end+ALT, 
      data=uvepgene, FUN=function(x) paste(unique(x), collapse="&"))
    colnames(VEPAnno)[7] <- "Gene"

    VEPAnno <- arrange(VEPAnno,id,chrn,start,end)

    ## Add annotation levels for reference
    VEPAnno$AnnoFact <- factor(VEPAnno$Annotation, levels=CLannlevs, ordered=TRUE)
    ## Save results
    vepout <- paste0(outdir, name,".RData")
    save(VEPAnno, file=vepout)
    return (VEPAnno)
}


veppass1 <- get_vep(anndir1)   #EnsGene: Ensemble for gene; Ensfeature: Ensemble for transcript
uvep1 <- get_anno(veppass1)
VEPAnno1 <- add_genename(uvep1,'VEP_WES2018')

veppass2 <- get_vep(anndir2)   
uvep2 <- get_anno(veppass2)
VEPAnno2 <- add_genename(uvep2,'Mouse_Sarcoma')

veppass3 <- get_vep(anndir3)   
uvep3 <- get_anno(veppass3)
VEPAnno3 <- add_genename(uvep3,'Lymphoma_Group3')


#######################################  Phenotype  ######################################
ssplit<-function(s, sp="\\-", n=1){
    sapply(s, function(x) strsplit(x, sp)[[1]][n])
}

pfiledir1 <- ###########
pfile1 <- file.path(pfiledir1, #########)
tools::md5sum(pfile1)
pdat1 <- data.frame(read_excel(pfile1))%>%
    mutate(sample1=ssplit(Code, n=2), 
           sample2=ssplit(Sample.type, sp="\\_", n=2), 
           type=factor(ifelse(grepl("tail", Sample.type), "tail",
                       ifelse(grepl("liver", Sample.type), "liver", Sample.type)), 
                       levels=c("lymphoma", "Sarcoma","tail", "liver")))%>%
    mutate(id=ifelse(is.na(sample2), sample1, sample2))%>% 
    dplyr::select(id, Genotype, type) %>%
    unique


pfiledir2 <- #############
pfile2 <- paste0(pfiledir2, ###########)
pdat2 <- data.frame(read_excel(pfile2, col_names=TRUE)) %>%
    mutate(id = substring(SAMPLE_ID,1, 4),
           Genotype = Sarcoma.genotype,
           type = 'Sarcoma') %>%
    dplyr::select(id, Genotype, type) %>%
    unique

pfiledir3 <- ###############
pfile3 <- paste0(pfiledir3, ###########)
pdat3 <- data.frame(read_excel(pfile3, col_names=TRUE)) %>%
         filter(Group_num == 'Group3') %>%
         mutate(id = as.character(Patient),
                Genotype = Groups,
                type='lymphoma') %>%
         dplyr::select(id, Genotype, type) %>%
         unique

#############################  merge  #####################################################
dat1 <- right_join(pdat1, VEPAnno1, by="id")
dat1_l <- filter(dat1, type == 'lymphoma')  
dat1_s <- filter(dat1, type == 'Sarcoma') 

dat2 <- right_join(pdat2, VEPAnno2, by="id")
dat3 <- right_join(pdat3, VEPAnno3, by="id")

#############################  Genelist for lymphoma & Sarcoma  ###########################
dat_s <- rbind(dat1_s, dat2)
dat_l <- rbind(dat1_l, dat3)
unique(dat_s$type);unique(dat_l$type)

gene_s <- dat_s %>% dplyr::select(id, Gene) %>% unique %>% 
          group_by(Gene) %>% count(Gene) %>% 
          arrange(desc(n)) %>% filter(n>1)
gene_s <- data.frame(gene_s)
#write.csv(gene_s, file.path(count_outdir, paste0("Sarcoma_genelist.csv")), quote=FALSE) #9344 genes

gene_l <- dat_l %>% dplyr::select(id, Gene) %>% unique %>% 
          group_by(Gene) %>% count(Gene) %>% 
          arrange(desc(n)) %>% filter(n>1)
gene_l <- data.frame(gene_l)
#write.csv(gene_l, file.path(count_outdir, paste0("lymphoma_genelist.csv")), quote=FALSE) #48

#############################  Count mutations  ###########################################
count_s <- dat_s %>% dplyr::select(id, Genotype) %>% 
          group_by(Genotype) %>% count %>% arrange(desc(n)) %>% mutate(type='Sarcoma')
count_s <- data.frame(count_s)

count_l <- dat_l %>% dplyr::select(id, Genotype,type) %>% 
          group_by(Genotype) %>% count %>% arrange(desc(n)) %>% mutate(type='lymphoma')
count_l <- data.frame(count_l)

counts <- rbind(count_s,count_l)
#write.csv(counts, file.path(count_outdir, paste0("mutation_counts.csv")), quote=FALSE)



