#########################
## Author: Xi Wang
## Date: 2018/03/01
## Title: Read in VEP tables and
##        generate mutation heatmap
#########################

library(dplyr)
library(biomaRt)#get gene name


#anndir <- "/data1/workspace/DCI/Kirsch/Chang-Lung.Lee/WES_2018/Results/VEP/vep_output/"
outdir <- "/data1/workspace/DCI/Kirsch/Chang-Lung.Lee/WES_2018/Analysis/Xi_Heatmap/"

## Function to return all duplicate rows based on a subset of columns
# Collapse rows to strings
squish <- function(x){
  y <- paste(x, collapse="")
  y <- gsub(" ", "", y, fixed=TRUE)
  return(y)
}
# Compare strings
rowcheck  <- function(df1, df2){
  xx <- apply(df1, 1, squish)
  yy <- apply(df2, 1, squish)
  zz <- xx %in% yy
  return(zz)
}
# Return matching rows
dfdup <- function(dat, cols){
  dat <- dat[,cols]
  dupdf <- dat[duplicated(dat),]
  return(rowcheck(dat, dupdf))
}


###################
##  ##  VEP  ##  ##
#Add Annotation for gene name and mutation types. 
#filter annotation types. order mutations.
###################

## Read in Jeremy's combined results
vepfiltfile <- '/home.local/jsg57/kirsch_project/group3_table.tsv'
tools::md5sum(vepfiltfile)
veppass <- read.table(vepfiltfile, sep=' ', header=FALSE, stringsAsFactors=FALSE)
dim(veppass)
# 1112395       9
colnames(veppass) <- c("id", "seqnames", "start", "end", "ALT", 
                        "EnsGene", "EnsFeature", 
                        "Annotation", "Impact")
head(veppass)
#EnsGene: Ensemble for gene; Ensfeature: Ensemble for transcript

## Drop rows differing by transcript ID
vepgen <- veppass %>% 
            dplyr::select(-id,-EnsFeature,-Impact) %>% 
            unique
dim(vepgen)
# 491732      6
uvep <- unique(filter(vepgen, EnsGene!='-'))
dim(uvep)
# 491221      6    Filter out those without Gene ID. 

## Reorder mutations by chr and pos
uvep$chrn <- case_when(uvep$seqnames %in% c("M","MT")~22,
                      uvep$seqnames=="Y"~21,
                      uvep$seqnames=="X"~20,
                      TRUE~as.numeric(uvep$seqnames))
table(uvep$chrn, uvep$seqnames, exclude=NULL)

snpord <- order(uvep$chrn, uvep$start, uvep$end)
uvep <- uvep[snpord,]
dim(uvep)
# 491221      7

## Expand comma-separated terms
sort(unique(uvep$Annotation))

multi <- grep(",", uvep$Annotation, fixed=T)
#if Annotation has ",", then there is more than one annotation. grep search for "," pattern.(if pattern shows, =1;otherwise,=0)
uvepSingle <- uvep[-multi,]
dim(uvepSingle)
#397734  With Single Annotation
uvepMulti <- uvep[multi,]
dim(uvepMulti)
#93487   With Multiple Annotation

splitAnn <- strsplit(uvepMulti$Annotation, split=',') 
uvepMultiExpanded <- NULL
for(i in 1:nrow(uvepMulti)){ 
  # Abuse merge to duplicate uvep rows
  uvepMultiExpanded <- rbind.data.frame(uvepMultiExpanded, 
        merge(x=uvepMulti[i,-6], y=splitAnn[[i]]), 
        stringsAsFactors=FALSE)
}
dim(uvepMultiExpanded)

colnames(uvepMultiExpanded)[7] <- "Annotation"
uvepExp <- rbind.data.frame(uvepSingle, uvepMultiExpanded,
            stringsAsFactors=FALSE)
dim(uvepExp)
# 588891  7

uvepExp <- unique(uvepExp)
dim(uvepExp)
# 517616

## Select terms based on Chang-Lung's preferred VEP terms (email 7/27/17)
CLannlevs <- c('transcript_ablation', 'frameshift_variant', 
  'start_lost', 'stop_gained', 
  'splice_region_variant',
  'splice_acceptor_variant',
  'splice_donor_variant', 
  'incomplete_terminal_codon_variant',
  'stop_lost', 'missense_variant',
  'inframe_deletion', 'inframe_insertion')

plotMuts$Type <- case_when(plotMuts$Annotation %in% 'frameshift_variant' ~ "Frame shift",
    plotMuts$Annotation %in% c('splice_acceptor_variant',
'splice_donor_variant', 'splice_region_variant') ~ "Splice site",
    plotMuts$Annotation %in% c('stop_gained', 'incomplete_terminal_codon_variant') ~ "Nonsense",
    plotMuts$Annotation %in% 'inframe_insertion' ~ "In frame",
    plotMuts$Annotation %in% 'missense_variant' ~ "Missense")


uvep <- filter(uvepExp, Annotation %in% CLannlevs)
dim(uvep)
# 63211     7

## Select most deleterious annotation within Gene
uvep$AnnoN <- as.numeric(factor(uvep$Annotation, 
        levels=CLannlevs, ordered=TRUE))
tab <- table(uvep$Annotation, uvep$AnnoN, useNA="always")
#rownames(tab)[which(tab[,ncol(tab)]!=0)] #chracter(0)

uvep <- data.frame(uvep %>% group_by(seqnames, start, end, ALT, EnsGene) %>% 
              filter(AnnoN == min(AnnoN)), stringsAsFactors=FALSE)
nrow(uvep)
# 61147
uvep <- unique(uvep);dim(uvep)
#61147
sum(dfdup(uvep, 1:5))
#0
sum(dfdup(uvep, 1:4))
#860

## Get Gene names
genesym <- unique(uvep$EnsGene)
length(genesym)
#14972

mouse <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
res <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"), 
             filters="ensembl_gene_id", values=genesym, mart=mouse)
dim(res)
#14940 

setdiff(genesym, res[,1])#Not the same
#[1] "ENSMUSG00000112508" "ENSMUSG00000094798" "ENSMUSG00000095982"
#[4] "ENSMUSG00000110859" "ENSMUSG00000044060" "ENSMUSG00000058923"
#[7] "ENSMUSG00000093785" "ENSMUSG00000093635" "ENSMUSG00000093692"
#[10] "ENSMUSG00000058731" "ENSMUSG00000093523" "ENSMUSG00000093411"
#[13] "ENSMUSG00000068426" "ENSMUSG00000093679" "ENSMUSG00000058588"
#[16] "ENSMUSG00000093501" "ENSMUSG00000093451" "ENSMUSG00000094206"
#[19] "ENSMUSG00000093764" "ENSMUSG00000057612" "ENSMUSG00000093632"
#[22] "ENSMUSG00000045713" "ENSMUSG00000057161" "ENSMUSG00000004613"
#[25] "ENSMUSG00000091095" "ENSMUSG00000057799" "ENSMUSG00000094600"
#[28] "ENSMUSG00000060490" "ENSMUSG00000071295" "ENSMUSG00000113525"
#[31] "ENSMUSG00000066944" "ENSMUSG00000063732"

#uvepgene <- left_join(uvep, res,  by=c("EnsGene"="ensembl_gene_id")) #61147
#remove those mutations without gene name
uvepgene <- inner_join(uvep, res,  by=c("EnsGene"="ensembl_gene_id"))
dim(uvepgene)
# 60953
#194 mutations are removed because their EnsGene do not have gene names.

sum(dfdup(uvepgene, 1:5))#0
sum(dfdup(uvepgene, 1:4))#814
#uvepgene[dfdup(uvepgene, 1:4),]

# Check if overlapping genes also occur outside of overlap
#which(uvepgene$external_gene_name %in% uvepgene[dfdup(uvepgene, 1:4),9])


## Aggregrate repeated rows (Gene and annotation)
# VEPAnno <- aggregate(uvepgene$external_gene_name, by=uvep[,c(1:4,6:8)], 
#                       paste, collapse="&")
VEPAnno <- aggregate(cbind(external_gene_name, Annotation) ~ chrn+seqnames+start+end+ALT, 
  data=uvepgene, FUN=function(x) paste(unique(x), collapse="&"))
dim(VEPAnno)
# 60481

colnames(VEPAnno)[6] <- "Gene"
length(grep("&", VEPAnno$Gene, fixed=T, value=T))#
grep("&", VEPAnno$Annotation, fixed=T, value=T)

snpord <- order(VEPAnno$chrn, VEPAnno$start, VEPAnno$end)
VEPAnno <- VEPAnno[snpord,]
dim(VEPAnno)
#60481
dim(unique(VEPAnno[,1:4]))
#60209 
length(unique(VEPAnno$Gene))
# 14884


## Add annotation levels for reference
VEPAnno$AnnoFact <- factor(VEPAnno$Annotation, levels=CLannlevs, ordered=TRUE)
levels(VEPAnno$AnnoFact)
table(VEPAnno$Annotation, as.numeric(VEPAnno$AnnoFact), exclude=NULL)
head(VEPAnno)

## Save results
vepout <- paste0(outdir, "VEP.RData")
save(VEPAnno, file=vepout)
tools::md5sum(vepout)
#   "07ef4040b4c9ea3b2136721c41248808"


################################
#########VCF ###################
################################
library(VariantAnnotation)#preprocess
#change sampdir to the Sarcoma_WES mutect2 directory.
sampdir <- "/data1/workspace/DCI/Kirsch/Chang-Lung.Lee/WES_2018/mutect2/output/"
anndir2 <- "/data1/workspace/DCI/Kirsch/Chang-Lung.Lee/WES_2018/Analysis/Xi_Heatmap/"

## Get list of VCF files
# Each file is in separate subdir
samples <- list.dirs(path=sampdir, full.names=F, recursive=F)
samples
length(samples)#19 samples
vcfs <- paste0(sampdir,samples, "/MuTect/", samples, "-MuTect-1.1.7.vcf")
vcfs[1]
all(file.exists(vcfs))

## Read in all mutations from vcf files
allMutations <- data.frame(NULL, stringsAsFactors=FALSE)
for(i in 1:length(samples)){
   vcf2 <- readVcf(vcfs[i])
   keep <- which(info(vcf2)$STR==FALSE) #STR:insertion
   muts <- as.data.frame(rowRanges(vcf2[keep,]), row.names=NULL)[,-c(5,6,9)] #remove strand paramRangeID QUAL
   muts$id <- samples[i] 
   allMutations <- rbind.data.frame(allMutations, filter(muts, FILTER=="PASS"& width==1))
}
dim(allMutations)

#aa<-unstrsplit(allMutations$ALT, sep=",")
allMutations$ALT <- unlist(lapply(allMutations$ALT, as.character))
allMutations$seqnames <- as.character(allMutations$seqnames)

table(allMutations$width, exclude=NULL)
table(allMutations$seqnames)
table(allMutations$id)

## Make chromosome:position key
allMutations$chrn <- case_when(allMutations$seqnames %in% c("M","MT")~22,
                             allMutations$seqnames=="Y"~21,
                             allMutations$seqnames=="X"~20,
                             TRUE~as.numeric(allMutations$seqnames))
table(allMutations$chrn, allMutations$seqnames, exclude=NULL)
allMutations$cpos <- paste0(allMutations$chrn,":",allMutations$start)
length(unique(allMutations$cpos))
#50481

##### Read in Annotation data #######
annfile <- paste0(anndir2, "VEP.RData") 
tools::md5sum(annfile)# the same
load(annfile) #load VEPAnno dataset

dim(VEPAnno)
#60481

CLannlevs2 <- levels(VEPAnno$AnnoFact)
CLannlevs2 == CLannlevs #all TRUE

VEPAnno$chrn <- case_when(VEPAnno$seqnames %in% c("M","MT")~22,
                            VEPAnno$seqnames=="Y"~21,
                            VEPAnno$seqnames=="X"~20,
                            TRUE~as.numeric(VEPAnno$seqnames))
table(VEPAnno$chrn, VEPAnno$seqnames, exclude=NULL)

VEPAnno$cpos <- paste0(VEPAnno$chrn,":",VEPAnno$start)

length(unique(VEPAnno$cpos))
# 60063

# Dadong didnt filter by PASS, 23511 was removed
length(setdiff(VEPAnno$cpos, unique(allMutations$cpos)))#41400
# Not annotated or non-protein altering mutations
length(setdiff(unique(allMutations$cpos), VEPAnno$cpos)) #31678

dim(unique(VEPAnno[,c("cpos","Gene")])) #60063
length(unique(VEPAnno$"cpos"))          #60063


## Merge mutations and function annotations
allMuts2 <- inner_join(allMutations[,c("id", "cpos", "start", "end", 
                          "REF", "ALT")], 
                       VEPAnno[,c("cpos","start", "end", "ALT","Gene","Annotation","AnnoFact")])
dim(allMuts2)
# 18562
length(unique(allMuts2$Gene))
# 10022

## Some results (insertions) not matching end/ALT
# VEP counts REF in insertions, readVcf does not
# The investigator is also interested at the insertions, so we keep them.
close <- filter(allMutations, allMutations$cpos %in% VEPAnno$cpos & 
                             !allMutations$cpos %in% allMuts2$cpos)
dim(close)#153 10
all(close$end-close$start == 0)#TRUE
length(unique(close[,10]))#151

caa<- filter(VEPAnno, cpos %in% close$cpos)
closeVEP <- filter(VEPAnno, cpos %in% close$cpos)
dim(closeVEP)#160. This is to verify that all the mutations are from VEP
all(closeVEP$end-closeVEP$start == 1)#FALSE there are insertions/deletions

close$end <- close$end+1
close$ALT <- substring(close$ALT,2)

closeMuts <- inner_join(close[,c("id", "cpos", "start", "end", "REF", "ALT")], 
                     closeVEP[,c("cpos","start", "end", "ALT","Gene","Annotation","AnnoFact")], by = c("cpos", "start", "end", "ALT"))
dim(closeMuts)
#149 9

allMuts2 <- rbind.data.frame(allMuts2, closeMuts)
dim(allMuts2)
# 18711

head(allMuts2)
levels(allMuts2$AnnoFact)

## Save annotated mutations
mutsfile <- paste0(anndir2, "ProtAltMuts.RData")
save(allMuts2, file=mutsfile)
tools::md5sum(mutsfile)
#"d5535269dec56107dd60c77607fd4efe"

################################
#make heatmaps
################################
library(openxlsx)
library(Matrix)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

mandir <- "/data1/workspace/DCI/Kirsch/Chang-Lung.Lee/WES_2018/Data/"
outdir <- "/data1/workspace/DCI/Kirsch/Chang-Lung.Lee/WES_2018/Analysis/Xi_Heatmap/"

datafile <- paste0(outdir, "ProtAltMuts.RData")
tools::md5sum(datafile)
#"d5535269dec56107dd60c77607fd4efe"
load(datafile)
dim(allMuts2)
# 18711

## Number of protein-altering mutations per sample
palt <- table(allMuts2$id)
#matrix(palt, dimnames=list(names(palt),NULL))
#     S1 S10 S11 S12 S13 S2 S27  S28 S29 S3 S30 S31 S32 S33  S34  S35  S36  S37
#[1,] 31  20  26  25  19 25  25 2448   5 24  18  39  40  19 2203 2644 2243 1653
#      S38  S39 S4 S40 S41 S42 S43 S44 S45 S5 S6 S7 S8 S9
#[1,] 3113 3863 13  29  18  12  12  15   5 17 22 37 31 17


## List of Genes mutated in more than one sample
tab2 <- table(allMuts2$Gene, allMuts2$id)

# Use >0 since some pts have more than one mutation in a single gene. At least 2 samples need to have the mutation 
multimuts2 <- rownames(tab2)[rowSums(tab2>0)>1]
count <- rowSums(tab2>0)[rowSums(tab2>0)>1]
multimuts2 <- cbind(multimuts2,count)
multimuts2 <- multimuts2[order(count,decreasing=TRUE),]
dim(multimuts2) 
#3820. the list of genes
#write.csv(multimuts2,file=paste0(outdir,"multimuts_genelist.csv"),quote=FALSE,row.names=FALSE,col.names=TRUE)





## Plot genes mutated in more than one (Group 1+3) sample + requesed genes from CL
plotgenes <- c("Ddx3x", "Gm16602", "mt-Nd5", "Naip5", "Notch1", "Ttn", "Zfp987", 
               "Ikzf1", "mt-Cytb", "mt-Nd2", "Fbxw7", "Naip6") 

counts2 <- filter(allMuts2, Gene %in% plotgenes)
dim(counts2)#125 9
table(counts2$id, counts2$Gene)



## Combine groups
allMuts <- allMuts2
dim(allMuts) #
CLannlevs <- levels(allMuts$AnnoFact)
CLannlevs

# ## Spreadsheet of mutations for CL
# ssmuts <- filter(allMuts, Gene %in% plotgenes)
# dim(ssmuts)
# ssmuts$chr <- as.numeric(sub(":.*", "", ssmuts$cpos))
# ssmuts <- ssmuts[order(ssmuts$id, ssmuts$chr, ssmuts$start, ssmuts$Gene, ssmuts$AnnoFact),]
#  # change chr to X, Y, M
# ssmuts$chr <- c(1:19,"X","Y","M")[ssmuts$chr]
# # Add type
# ssmuts$impact <- ifelse(as.numeric(ssmuts$AnnoFact) <= 9, "High", "Moderate")
# table(ssmuts$Annotation, ssmuts$impact, exclude=NULL)
# outcols <- c("id","chr","start","end","REF","ALT","Gene","Annotation","impact")
# dim(ssmuts[,outcols])
# sheet <- paste0(outdir, "Mutations_in_plotted_genes.xlsx")
# write.xlsx(list("Sarcoma"=ssmuts[,outcols]), file=sheet, rowNames=FALSE)
# tools::md5sum(sheet)#73d53ee93c9fb42db0c72be02fbd6217

## Collapse Multiple positions into single genes for plot
plotMuts <- unique(filter(allMuts, Gene %in% plotgenes)[, c("id","Gene","Annotation","AnnoFact")])
dim(plotMuts)
# 37 4

## Select most deleterious annotation within Gene
plotMuts$AnnoN <- as.numeric(plotMuts$AnnoFact)
plotMuts <- data.frame(plotMuts %>% group_by(Gene, id) %>% 
              filter(AnnoN == min(AnnoN)), stringsAsFactors=FALSE)
dim(plotMuts)
# 27 5

# ## KO QC: Are all mutations SNPs? YES. all equal 1
# table(nchar(allMuts$REF), exclude=NULL)
 
# snpsQ <- unique(filter(allMuts, Gene %in% plotgenes)[, c("id","Gene","Annotation","ALT")])
# dim(snpsQ)
# # 63  4
# 
# snpsQ <- inner_join(snpsQ, plotMuts,by = c("id", "Gene", "Annotation"))
# dim(snpsQ)
# #37 6
# table(nchar(snpsQ$ALT), exclude=NULL)
# #  1
# # 37 
# # 37 SNPs 
# ## KO QC: REF==ALT? No.
# any(allMuts$REF == allMuts$ALT) #FALSE


## Collapse annotation terms per Chang-Lung's preferences
#paste(CLannlevs, collapse="','")

#plotMuts$Type <- case_when(plotMuts$Annotation %in% 'frameshift_variant' ~ "Frame shift",
#    plotMuts$Annotation %in% c('splice_acceptor_variant',
#'splice_donor_variant', 'splice_region_variant') ~ "Splice site",
#    plotMuts$Annotation %in% c('stop_gained', 'incomplete_terminal_codon_variant') ~ "Nonsense",
#    plotMuts$Annotation %in% 'inframe_insertion' ~ "In frame",
#    plotMuts$Annotation %in% 'missense_variant' ~ "Missense")
#table(plotMuts$Annotation, plotMuts$Type, exclude=NULL)


plotMuts$Type <- case_when(
  plotMuts$Annotation %in% c('transcript_ablation','frameshift_variant',
    'start_lost','stop_gained','splice_region_variant','splice_acceptor_variant',
    'splice_donor_variant','incomplete_terminal_codon_variant',
    'stop_lost') ~ "High impact",
  plotMuts$Annotation %in% c('missense_variant','inframe_deletion',
    'inframe_insertion') ~ "Moderate impact")

table(plotMuts$Annotation, plotMuts$Type, useNA='always')
#                        High impact Moderate impact <NA>
#  missense_variant                0              17    0
#  splice_donor_variant            1               0    0
#  splice_region_variant          level 3          0    0
#  stop_gained                     6               0    0
#  <NA>                            0               0    0


## Make sparse matrix of mutations by pts
# to order genes by chr:pos
genedf <- unique(filter(allMuts[,c("cpos","start","Gene")], Gene %in% plotgenes))
genedf <- genedf[!duplicated(genedf$Gene),]
genedf$chr <- sub(":.*", "", genedf$cpos)
genedf <- genedf[order(as.numeric(genedf$chr), as.numeric(genedf$start)),]

genefactorlevels <- genedf$Gene
genefactorlevels #"Notch1" "Ttn"    "Fbxw7"  "Ikzf1"  "Naip5"  "Naip6"  "Ddx3x"  "mt-Nd2"

patfact <- as.factor(plotMuts[,"id"])
snpfact <- factor(plotMuts[,"Gene"], levels=genefactorlevels, ordered=TRUE)
mutfact <- factor(plotMuts[,"Type"], levels=c("High impact","Moderate impact"), ordered=TRUE)
#   levels=c("Frame shift", "Splice site", "Nonsense", "In frame", "Missense"), ordered=TRUE)

## Some samples have no mutations in the selected genes. "4017" "4018" "4019" "4020" "4021"
missamp <- setdiff(unique(allMuts$id),levels(patfact))
missamp

i <- as.numeric(snpfact)
j <- as.numeric(patfact)+length(missamp) # Add enmpty columns
x <- as.numeric(mutfact)

mapmat <- sparseMatrix(i=i,j=j,x=x) # SNPs in rows to make room for legend
rownames(mapmat) <- levels(snpfact)
colnames(mapmat) <- c(missamp,levels(patfact))
levels(mutfact) #"High impact"     "Moderate impact"

mapmat <- mapmat[,order(colnames(mapmat))]

dim(mapmat) # 8 19
sum(!mapmat==0) #27


## Read in sample subgroups
manfile <- paste0(mandir, "Lee_WES_samplelist_Group4_20180202.xlsx")
man <- read.xlsx(manfile, sheet=1, startRow=1, colNames=TRUE, rowNames=FALSE)
dim(man)#64 4
man$Patient <- substring(man$Code,6, 9) #extract patient ID
man$Groups <- man$Genotype

sampdir <- "/data1/workspace/DCI/Kirsch/Chang-Lung.Lee/WES_2018/mutect2/output/"
samples <- list.dirs(path=sampdir, full.names=F, recursive=F)
samples

man <- filter(man, Patient %in% samples)
dim(man)#32
ptgrps <- unique(man[,c("Patient","Groups")])
dim(ptgrps)#32

ptgrps <- ptgrps[ptgrps$Patient %in% colnames(mapmat),]
table(ptgrps$Groups, exclude=NULL)
#Kras mut + p53 mut      MCA + p53 mut       MCA + p53 WT 
#                 7                  5                  7 

##############Find top 10 recurring genes for each group###################
tab2 <- as.data.frame.matrix(tab2)
allpts <- as.data.frame(ptgrps)
allpts$Groups <- as.factor(allpts$Groups)
allgenes <- tab2[rowSums(tab2>0)>1,] #show in more than one sample

getgenes <- function(level){
grp1pt <- allpts$Patient[allpts$Groups==level]
group1 <- allgenes %>% dplyr::select(grp1pt)
grp1gene <- sort(rowSums(group1),decreasing=TRUE)
}

group1<-getgenes(level="Kras mut + p53 mut")
as.data.frame(group1[1:15])
group2<-getgenes(level="MCA + p53 mut")
as.data.frame(group2[1:15])
group3<-getgenes(level="MCA + p53 WT")
as.data.frame(group3[1:15])

getgenes2 <- function(level){
grp1pt <- allpts$Patient[allpts$Groups==level]
group1 <- allgenes %>% dplyr::select(grp1pt)
grp1gene <- sort(rowSums(group1>0),decreasing=TRUE)
}

group11<-getgenes2(level="Kras mut + p53 mut")
as.data.frame(group11[1:15])
group22<-getgenes2(level="MCA + p53 mut")
as.data.frame(group22[1:15])
group33<-getgenes2(level="MCA + p53 WT")
as.data.frame(group33[1:15])



##############Total mutation number for each sample########################
#######not filter mutations for more than 1 sample
mutsum1 <- colSums(tab2)
#######filter mutations for more than 1 sample
mutsum2 <- colSums(allgenes)
colnames(mutsum1)==colnames(mutsum2)
mutsum <- rbind(mutsum1,mutsum2)
rownames(mutsum) <- c("not_filtered","filtered")
mutsum <- t(mutsum)
mutsum2 <- cbind(Patient=rownames(mutsum),mutsum)
mutsum3 <- merge(mutsum2,ptgrps)
write.xlsx(mutsum3,file=paste0(outdir,"mutsum.xlsx"),quote=FALSE,row.names=TRUE,col.names=TRUE)

## Reorder patients by Group
ptgrpord <- c("MCA + p53 WT","MCA + p53 mut","Kras mut + p53 mut")
idord <- order(factor(ptgrps$Groups, levels=ptgrpord, ordered=T), ptgrps$Patient)
ptgrps <- ptgrps[idord,]
ptgrps

annpts <- data.frame(Group=ptgrps$Groups)
rownames(annpts) <- ptgrps$Patient
annpts

## Generate Heatmap ##
## Colors for Mutations
levels(mutfact)

# red = frame shift mutation
# orange = splice site
# yellow = nonsense
# Green = in frame
# Blue = missense

cols1 <- brewer.pal(9,"Set1") 
# plot(1:9,col=cols1,pch=16)
nun0white <- "#FFFFFF"
# nablack <- "#000000"
frm1 <- cols1[1] #red
# spl2 <- cols1[5] #orange
# non3 <- cols1[6] #yellow
# inf4 <- cols1[3] #green
mis5 <- cols1[2] #blue

# mapcols <- c(nun0white, frm1, spl2, non3, inf4, mis5)
mapcols <- c(nun0white, frm1, mis5)
cbind(c("",levels(mutfact)), mapcols)

## Colors for pt groups
cols2 <- brewer.pal(12,"Paired")
#  plot(1:8,col=cols2,pch=16) 
ptcols <- c(cols2[5], cols2[3], cols2[7])
names(ptcols) <- ptgrpord
ptcols

library(grid)

pngname <- paste0(outdir,"GenesOfInterest.png")
pngw <- 4
pngh <- 3
png(pngname, width=pngw, height=pngh, units="in", res=300)
# use pheatmap(filename) first, then measure results to get png dimensions
pheatmap(t(mapmat[,idord]), cluster_rows=FALSE, cluster_cols=FALSE, 
  annotation_row=annpts,fontsize=6.5,
  #main="Genes with mutations in \n≥ two samples per group", # Warnings from use of "≥"
  main="Mutations in selected genes",
  cellwidth=9, cellheight=9, 
  color=mapcols, legend=FALSE, 
  annotation_colors=list(Group=ptcols)
#   , filename=pngname
)

# Add mutation type legend
# zoom in to check alignment/sizing
ledgX <- 0.554
ledgY <- 0.73
boxhight <- 0.045
boxwidth <- boxhight*pngh/pngw
textsize <- 0.55
grid.text("Variant", ledgX, ledgY, just='left', vjust=1,
  gp=gpar(fontface="bold",cex=textsize))
ledgY <- ledgY - 0.02 # Spacing between legend name and boxes
for(m in 1:(length(mapcols)-1)){
  grid.rect(x=ledgX, y=ledgY-boxhight*m, width=boxwidth, height=boxhight,
    just='left', vjust=0.5,
    gp=gpar(col="darkgrey",fill=mapcols[m+1]))
  grid.text(levels(mutfact)[m], ledgX+boxwidth+.01, ledgY-boxhight*m,
    just='left', vjust=0.5, gp=gpar(cex=textsize))
}
graphics.off()

tools::md5sum(pngname)
#71747ad50cd0fb6afba1d463ebcc2e01
