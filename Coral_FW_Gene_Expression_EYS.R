# source("http://bioconductor.org/biocLite.R")
# biocLite("DESeq2")

#set your working directory
setwd("~/Documents/GitHub/Assign2_Coral_FW/") #you will need to change to your own directory

###Install packages
BiocManager::install("DESeq2")
#BiocManager::install("phangorn")
install.packages("backports")
install.packages("caTools")
BiocManager::install("affycoretools")
BiocManager::install("arrayQualityMetrics")
install.packages("pheatmap")

###conduct array quality metrics to detect and remove outliers
library(DESeq2)
library(affycoretools)
library(arrayQualityMetrics)
library(genefilter)
library(Biobase)
#read in counts 
c
head(countData)
length(countData[,1])
names(countData)=c("fav_recoveryA",	"fav_recoveryB", "fav_recoveryC","fav_stressA", "fav_stressB", "fav_stressC")
row.names(countData)=sub("", "isogroup", rownames(countData))
head(countData)


setwd("~/Documents/GitHub/Assign2_Coral_FW/")
v=setwd("~/Documents/GitHub/Assign2_Coral_FW/")
# # # #look for outliers
treat=c("fav_recovery",	"fav_recovery", "fav_recovery","fav_stress", "fav_stress", "fav_stress")
genotype = c("A", "B", "C", "A", "B", "C")
g = data.frame(cbind(treat, genotype))
g
colData= g

dds=DESeqDataSetFromMatrix(countData=countData,
                           colData = g,
                           design = ~treat)

vsd.ge=assay(vst(dds))
rl=vst(dds)
e=ExpressionSet(assay(rl), AnnotatedDataFrame(as.data.frame(colData(rl))))
arrayQualityMetrics(e,outdir=v,intgroup=c("treat"),force=T)
# dev.off() for windows only
# double-click index.html

#####you only ever need to run the above code once. Outliers are decided at the beginning. 
## I like to close R and restart with packages etc
##So, please save your script and restart R
#no outliers detected!
# if there's outliers we need to remove them

setwd("~/Documents/GitHub/Assign2_Coral_FW/")
library("DESeq2")
library("ggplot2")

#read in counts 
countData <- read.table("allcounts_Harvey_fav.txt")
head(countData)
length(countData[,1])
#19937

names(countData)=c("fav_recoveryA",	"fav_recoveryB", "fav_recoveryC","fav_stressA", "fav_stressB", "fav_stressC")
row.names(countData)=sub("", "isogroup", rownames(countData))
head(countData)

totalCounts=colSums(countData)
totalCounts
barplot(totalCounts, col=c("coral", "coral", "coral", "red", "red", "red", "blue", "blue", "blue"), ylab="raw counts")

min(totalCounts) #458772
max(totalCounts)  #6596244

treat=c("fav_recovery",	"fav_recovery", "fav_recovery","fav_stress", "fav_stress", "fav_stress")
replicate = c("A", "B", "C", "A", "B", "C")
g = data.frame(cbind(treat, replicate))
g
colData<- g

dds<-DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~treat+genotype) #can only test for the main effects of site, pco2, temp

#one step DESeq
dds<-DESeq(dds)
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing

head(dds)
res<- results(dds)

#Look at dispersions plot
plotDispEsts(dds, main="Dispersion plot Snails")

####################stress vs. recovery pairwise comparisons
colData$treat<-factor(colData$treat, levels=c("fav_stress", "fav_recovery"))
##second term is the "control"
resstress <- results(dds, contrast=c("treat","fav_stress", "fav_recovery"))
#how many FDR < 10%?
table(resstress$padj<0.1)
# 0.1=46
# 0.05=28
# 0.01=19
summary(resstress)

nrow(resstress[resstress$padj<0.05 & !is.na(resstress$padj),])  # Num significantly differentially expressed genes excluding the no/low count genes   #28

dev.off()
plotMA(resstress, main="stress vs recovery")
plotMA(resstress, main="stress vs recovery", ylim=c(-2,2))

results <- as.data.frame(resstress)
head(results)

nrow(resstress[resstress$padj<0.1 & resstress$log2FoldChange > 0 & !is.na(resstress$padj),])
nrow(resstress[resstress$padj<0.1 & resstress$log2FoldChange < 0 & !is.na(resstress$padj),])
#UP in stress 31
#DOWN in stress 15

write.table(resstress, file="stress.txt", quote=F, sep="\t")

cd <- read.table("stress.txt")
head(cd)

##make the GO table for MWU
head(cd)

library(dplyr)
cd
go_input_stress = cd %>%
  tibble::rownames_to_column(var = "iso") %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(iso, mutated_p_updown)
head(go_input_stress)
colnames(go_input_stress) <- c("gene", "pval")
head(go_input_stress)
write.csv(go_input_stress, file="stress_GO.csv", quote=F, row.names=FALSE)
# gene name is isogroupisogroup (has one extra isogroup)


###############################################################################################
##############################################################################
#--------------get pvals
valstress=cbind(resstress$pvalue, resstress$padj)
head(valstress)
colnames(valstress)=c("pval.stress", "padj.sterss")
length(valstress[,1])
table(complete.cases(valstress))


######-------------make rlogdata and pvals table
rlog=rlogTransformation(dds, blind=TRUE) 
rld=assay(rlog)
head(rld)
colnames(rld)=paste(colData$treat)
head(rld)
length(rld[,1])

rldpvals=cbind(rld,valstress)
head(rldpvals)
dim(rldpvals)
#[1] 19937     8
table(complete.cases(rldpvals))
#FALSE  TRUE 
#15353  4584 

write.csv(rldpvals, "stress_RLDandPVALS.csv", quote=F)

colnames(rld)=paste(colData$treat)
head(rld)

library(RColorBrewer)
# Sample distance heatmap
sampleDists <- as.matrix(dist(t(rld)))
library(gplots)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          margin=c(10, 10), main="Sample Distance Matrix")



#################################################################################
# VENN Diagram to include both up and down regulated genes in common for PC02
library(VennDiagram)

stress_up=row.names(resstress[resstress$padj<0.1 & !is.na(resstress$padj) & resstress$log2FoldChange>0,])
length(p76_up) #66
p76_down=row.names(resstress[resstress$padj<0.1 & !is.na(resstress$padj) & resstress$log2FoldChange<0,])
length(p76_down) #420
p75_up=row.names(resph75[resph75$padj<0.1 & !is.na(resph75$padj) & resph75$log2FoldChange>0,])
length(p75_up) #38
p75_down=row.names(resph75[resph75$padj<0.1 & !is.na(resph75$padj) & resph75$log2FoldChange<0,])
length(p75_down) #80

p76=row.names(resstress[resstress$padj<0.1 & !is.na(resstress$padj),])
p75=row.names(resph75[resph75$padj<0.1 & !is.na(resph75$padj),])

#UP
pdegs05_up=union(p76_up,p75_up)
length(pdegs05_up)
#93

#DOWN
pdegs05_down=union(p76_down,p75_down)
length(pdegs05_down)
#432

#ALL
pdegs05=union(p76,p75)
length(pdegs05)
#524

###do UP, DOWN, ALL
candidates=list("7.6"=p76, "7.5"=p75)
quartz()
prettyvenn=venn.diagram(
  x = candidates,
  filename=NULL,
  col = "transparent",
  fill = c("coral2", "forestgreen"),
  alpha = 0.5,
  # label.col = c("darkred", "white", "darkgreen", "white", "white", "white", "blue4"),
  cex = 2.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col = c("darkred", "darkgreen"),
  cat.cex = 2.5,
  cat.fontfamily = "sans",
  cat.dist = c(0.08, 0.08),
  cat.pos = 1
);
grid.draw(prettyvenn)

###########################heat map of sample distances for pco2
rldpvals <- read.csv(file="Crep2016_RLDandPVALS.csv", row.names=1)
head(rldpvals)
rld=rldpvals[,1:9]
head(rld)

sampleDists <- dist(t(rld))
sampleDistMatrix <- as.matrix( sampleDists )
treat=c( "pH7.5", "pH7.5", "pH7.5", "pH7.6", "pH7.6", "pH7.6", "pH8", "pH8",  "pH8")
colnames(sampleDistMatrix)=paste(treat)
rownames(sampleDistMatrix)=paste(treat)

library(pheatmap)
heat.colors = colorRampPalette(rev(c("blue","yellow","red")),bias=0.3)(100)
pheatmap(sampleDistMatrix,color = heat.colors,cex=0.9,border_color=NA,cluster_rows=T,cluster_cols=T)

library(vegan)
library(ggplot2)
library(ggrepel)
library(tidyverse)

rld_t=t(rld)
pca <- prcomp(rld_t,center = TRUE, scale. = TRUE)
head(pca)
li <- pca$sdev^2 / sum(pca$sdev^2)
pc1v <- round(li[1] * 100, 1)
pc2v <- round(li[2] * 100, 1)
pca_s <- as.data.frame(pca$x)
head(pca_s)
pca_s <- pca_s[,c(1,2)]
pca_s$Samples = row.names(pca_s)
pca_s$treat=colData$treat
head(pca_s)

cbPalette <- c("darkgoldenrod2",  "darkolivegreen3", "dodgerblue3")
ggplot(pca_s, aes(PC1, PC2, color = treat, pch = treat)) +
  geom_point(size=3) +
  #  geom_text_repel(aes(label=Samples)) +
  scale_colour_manual(values=cbPalette)+
  theme_bw() +
  # geom_density2d(alpha=.5)+
  geom_polygon(alpha=.2)+
  xlab(paste0("PC1: ",pc1v,"% variance")) +
  ylab(paste0("PC2: ",pc2v,"% variance")) 
head(pca)
adonis(pca$x ~ treat, data = pca_s, method='eu', na.rm = TRUE)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# treat      2     40140 20069.9  13.048 0.81306  0.003 **
# Residuals  6      9229  1538.1         0.18694          
# Total      8     49369                 1.00000 


###################################heatmaps for genes NS vs FR
rldpvals <- read.csv(file="Crep2016_RLDandPVALS.csv", row.names=1)
head(rldpvals)
rld_site= rldpvals[,1:9]
head(rld_site)
gg=read.table("Crep454_iso2gene.tab",sep="\t", quote="", row.names=1)
head(gg)

nrow(rldpvals[rldpvals$padj.76<0.01& !is.na(rldpvals$padj.76),])
#242

topnum= 100 # number of DEGS
head(rldpvals)
top100=head(rldpvals[order(rldpvals$padj.76), ],topnum)
head(top100)
length(top100[,1])
summary(top100)
###
library(pheatmap)
head(top100)
p.val=0.1 # FDR cutoff
conds=top100[top100$padj.76<=p.val & !is.na(top100$padj.76),]
length(conds[,1])

exp=conds[,1:9] # change numbers to be your vsd data columns
means=apply(exp,1,mean) # means of rows
explc=exp-means # subtracting them
head(explc)

ccol=colorRampPalette(rev(c("red","chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
col0=colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)

pheatmap(explc,cluster_cols=T,scale="row",color=col0, show_rownames = F)

###################################Heatmap for the genes in common
rldpvals <- read.csv(file="Crep2016_RLDandPVALS.csv", row.names=1)
head(rldpvals)
p.val=0.1 # FDR cutoff
conds=rldpvals[rldpvals$padj.76<=p.val & !is.na(rldpvals$padj.76) & rldpvals$padj.75<=p.val & !is.na(rldpvals$padj.75),]
rld_data= conds[,c(1:9)]
head(rld_data)
nrow(rld_data)
gg=read.table("Crep454_iso2gene.tab",sep="\t", row.names=1)
library(pheatmap)
means=apply(rld_data,1,mean) # means of rows
explc=rld_data-means # subtracting them

ccol=colorRampPalette(rev(c("red","chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
col0=colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)

pheatmap(explc,cluster_cols=T,scale="row",color=col0, show_rownames = F)

# Make annotation table for pheatmap
ann = data.frame(cond = c('7.5', '7.5', '7.5', '7.6', '7.6', '7.6', '8', '8', '8'))
rownames(ann) <- names(explc)

# Set colors
Var1        <- c("darkgoldenrod2",  "darkolivegreen3", "dodgerblue3")
names(Var1) <- c("7.5", "7.6", "8")
anno_colors <- list(cond = Var1)

pheatmap(as.matrix(explc),annotation_col=ann,annotation_colors=anno_colors,cex=1.2,color=col0,border_color=NA,clustering_distance_rows="correlation",clustering_distance_cols="correlation", show_rownames=T)
