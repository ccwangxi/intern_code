#########################
## Author: XXXXXX
## Date: 2018/03/01
## Title: Read in VEP tables and
##        generate mutation heatmap
#########################

library(dplyr)
library(biomaRt)#get gene name


#anndir <- ##########
outdir <- ###########

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

## Read in combined results
vepfiltfile <- XXXXXXXXXXXXX
tools::md5sum(vepfiltfile)
veppass <- read.table(vepfiltfile, sep=' ', header=FALSE, stringsAsFactors=FALSE)
dim(veppass)
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

uvep <- unique(filter(vepgen, EnsGene!='-'))
dim(uvep)
#   Filter out those without Gene ID. 

## Reorder mutations by chr and pos
uvep$chrn <- case_when(uvep$seqnames %in% c("M","MT")~22,
                      uvep$seqnames=="Y"~21,
                      uvep$seqnames=="X"~20,
                      TRUE~as.numeric(uvep$seqnames))
table(uvep$chrn, uvep$seqnames, exclude=NULL)

snpord <- order(uvep$chrn, uvep$start, uvep$end)
uvep <- uvep[snpord,]
dim(uvep)


## Expand comma-separated terms
sort(unique(uvep$Annotation))

multi <- grep(",", uvep$Annotation, fixed=T)
#if Annotation has ",", then there is more than one annotation. grep search for "," pattern.(if pattern shows, =1;otherwise,=0)
uvepSingle <- uvep[-multi,]
dim(uvepSingle)
#With Single Annotation
uvepMulti <- uvep[multi,]
dim(uvepMulti)
#With Multiple Annotation

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

## Select terms based on preferred VEP terms 
CLannlevs <- XXXXXXXXXXXX

plotMuts$Type <- case_when(plotMuts$Annotation %in% 'AA' ~ "aa",
    plotMuts$Annotation %in% c('BB',
'CC', 'DD') ~ "bcd")


uvep <- filter(uvepExp, Annotation %in% CLannlevs)
dim(uvep)

## Select most deleterious annotation within Gene
uvep$AnnoN <- as.numeric(factor(uvep$Annotation, 
        levels=CLannlevs, ordered=TRUE))
tab <- table(uvep$Annotation, uvep$AnnoN, useNA="always")
#rownames(tab)[which(tab[,ncol(tab)]!=0)] #chracter(0)

uvep <- data.frame(uvep %>% group_by(seqnames, start, end, ALT, EnsGene) %>% 
              filter(AnnoN == min(AnnoN)), stringsAsFactors=FALSE)
nrow(uvep)

uvep <- unique(uvep);dim(uvep)

sum(dfdup(uvep, 1:5))

sum(dfdup(uvep, 1:4))


## Get Gene names
genesym <- unique(uvep$EnsGene)
length(genesym)


mouse <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
res <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"), 
             filters="ensembl_gene_id", values=genesym, mart=mouse)
dim(res)


setdiff(genesym, res[,1])#Not the same


#uvepgene <- left_join(uvep, res,  by=c("EnsGene"="ensembl_gene_id")) #61147
#remove those mutations without gene name
uvepgene <- inner_join(uvep, res,  by=c("EnsGene"="ensembl_gene_id"))
dim(uvepgene)
# mutations are removed because their EnsGene do not have gene names.

sum(dfdup(uvepgene, 1:5))#
sum(dfdup(uvepgene, 1:4))#
## Aggregrate repeated rows (Gene and annotation)

VEPAnno <- aggregate(cbind(external_gene_name, Annotation) ~ chrn+seqnames+start+end+ALT, 
  data=uvepgene, FUN=function(x) paste(unique(x), collapse="&"))
dim(VEPAnno)

colnames(VEPAnno)[6] <- "Gene"
length(grep("&", VEPAnno$Gene, fixed=T, value=T))#
grep("&", VEPAnno$Annotation, fixed=T, value=T)

snpord <- order(VEPAnno$chrn, VEPAnno$start, VEPAnno$end)
VEPAnno <- VEPAnno[snpord,]
dim(VEPAnno)

dim(unique(VEPAnno[,1:4]))

length(unique(VEPAnno$Gene))



## Add annotation levels for reference
VEPAnno$AnnoFact <- factor(VEPAnno$Annotation, levels=CLannlevs, ordered=TRUE)
levels(VEPAnno$AnnoFact)
table(VEPAnno$Annotation, as.numeric(VEPAnno$AnnoFact), exclude=NULL)
head(VEPAnno)

## Save results
vepout <- paste0(outdir, "VEP.RData")
save(VEPAnno, file=vepout)
tools::md5sum(vepout)



################################
#########VCF ###################
################################
library(VariantAnnotation)#use readVcf
#change sampdir to the Sarcoma_WES mutect2 directory.
sampdir <- XXXXXXXXXXXXXXXXXXX
anndir2 <- XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

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


##### Read in Annotation data #######
annfile <- paste0(anndir2, "VEP.RData") 
tools::md5sum(annfile)# the same
load(annfile) #load VEPAnno dataset

dim(VEPAnno)


CLannlevs2 <- levels(VEPAnno$AnnoFact)
CLannlevs2 == CLannlevs #all TRUE

VEPAnno$chrn <- case_when(VEPAnno$seqnames %in% c("M","MT")~22,
                            VEPAnno$seqnames=="Y"~21,
                            VEPAnno$seqnames=="X"~20,
                            TRUE~as.numeric(VEPAnno$seqnames))
table(VEPAnno$chrn, VEPAnno$seqnames, exclude=NULL)

VEPAnno$cpos <- paste0(VEPAnno$chrn,":",VEPAnno$start)

length(unique(VEPAnno$cpos))


# didnt filter by PASS
length(setdiff(VEPAnno$cpos, unique(allMutations$cpos)))#
# Not annotated or non-protein altering mutations
length(setdiff(unique(allMutations$cpos), VEPAnno$cpos)) #

dim(unique(VEPAnno[,c("cpos","Gene")])) #
length(unique(VEPAnno$"cpos"))          #


## Merge mutations and function annotations
allMuts2 <- inner_join(allMutations[,c("id", "cpos", "start", "end", 
                          "REF", "ALT")], 
                       VEPAnno[,c("cpos","start", "end", "ALT","Gene","Annotation","AnnoFact")])
dim(allMuts2)
# 
length(unique(allMuts2$Gene))
# 

## Some results (insertions) not matching end/ALT
# VEP counts REF in insertions, readVcf does not
# The investigator is also interested at the insertions, so we keep them.
close <- filter(allMutations, allMutations$cpos %in% VEPAnno$cpos & 
                             !allMutations$cpos %in% allMuts2$cpos)
dim(close)#
all(close$end-close$start == 0)#TRUE
length(unique(close[,10]))#

caa<- filter(VEPAnno, cpos %in% close$cpos)
closeVEP <- filter(VEPAnno, cpos %in% close$cpos)
dim(closeVEP)# This is to verify that all the mutations are from VEP
all(closeVEP$end-closeVEP$start == 1)#FALSE there are insertions/deletions

close$end <- close$end+1
close$ALT <- substring(close$ALT,2)

closeMuts <- inner_join(close[,c("id", "cpos", "start", "end", "REF", "ALT")], 
                     closeVEP[,c("cpos","start", "end", "ALT","Gene","Annotation","AnnoFact")], by = c("cpos", "start", "end", "ALT"))
dim(closeMuts)

allMuts2 <- rbind.data.frame(allMuts2, closeMuts)
dim(allMuts2)

head(allMuts2)
levels(allMuts2$AnnoFact)

## Save annotated mutations
mutsfile <- paste0(anndir2, "ProtAltMuts.RData")
save(allMuts2, file=mutsfile)
tools::md5sum(mutsfile)

################################
#make heatmaps
################################
library(openxlsx)
library(Matrix)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

mandir <- XXXXXXXXXXXXX
outdir <- XXXXXXXXXXXXXX

datafile <- paste0(outdir, "ProtAltMuts.RData")
tools::md5sum(datafile)

load(datafile)
dim(allMuts2)


## Number of protein-altering mutations per sample
palt <- table(allMuts2$id)
#matrix(palt, dimnames=list(names(palt),NULL))


## List of Genes mutated in more than one sample
tab2 <- table(allMuts2$Gene, allMuts2$id)

# Use >0 since some pts have more than one mutation in a single gene. At least 2 samples need to have the mutation 
multimuts2 <- rownames(tab2)[rowSums(tab2>0)>1]
count <- rowSums(tab2>0)[rowSums(tab2>0)>1]
multimuts2 <- cbind(multimuts2,count)
multimuts2 <- multimuts2[order(count,decreasing=TRUE),]
dim(multimuts2) 
#the list of genes
#write.csv(multimuts2,file=paste0(outdir,"multimuts_genelist.csv"),quote=FALSE,row.names=FALSE,col.names=TRUE)





## Plot genes mutated in more than one (Group 1+3) sample + requesed genes from investigator
plotgenes <- Xxxxxxx

counts2 <- filter(allMuts2, Gene %in% plotgenes)
dim(counts2)#125 9
table(counts2$id, counts2$Gene)



## Combine groups
allMuts <- allMuts2
dim(allMuts) #
CLannlevs <- levels(allMuts$AnnoFact)
CLannlevs


## Collapse Multiple positions into single genes for plot
plotMuts <- unique(filter(allMuts, Gene %in% plotgenes)[, c("id","Gene","Annotation","AnnoFact")])
dim(plotMuts)


## Select most deleterious annotation within Gene
plotMuts$AnnoN <- as.numeric(plotMuts$AnnoFact)
plotMuts <- data.frame(plotMuts %>% group_by(Gene, id) %>% 
              filter(AnnoN == min(AnnoN)), stringsAsFactors=FALSE)
dim(plotMuts)

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

## Some samples have no mutations in the selected genes. 
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

dim(mapmat) # 
sum(!mapmat==0) #


## Read in sample subgroups
manfile <- paste0(mandir, "ZXX.xlsx")
man <- read.xlsx(manfile, sheet=1, startRow=1, colNames=TRUE, rowNames=FALSE)

man$Patient <- substring(man$Code,6, 9) #extract patient ID
man$Groups <- man$Genotype

sampdir <- XXXXXXXXXXXXX
samples <- list.dirs(path=sampdir, full.names=F, recursive=F)
samples

man <- filter(man, Patient %in% samples)
dim(man)#
ptgrps <- unique(man[,c("Patient","Groups")])
dim(ptgrps)#

ptgrps <- ptgrps[ptgrps$Patient %in% colnames(mapmat),]
table(ptgrps$Groups, exclude=NULL)


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

# red = 
# orange = 
# yellow = 
# Green =
# Blue = 

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

