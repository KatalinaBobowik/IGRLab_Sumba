### File created 07.05.17
### The objective of this script is to create an RNASeq pipeline for the Yamagishi et al 2014 reads from the paper: Interactive transcriptome analysis of malaria patients and infecting Plasmodium falciparum
### Script written following the "RNA-seq analysis is easy as 1-2-3 with limma, Glimma and edgeR [version 2; referees: 3 approved]" pipeline
### PCA data closely modelled after Irene Gallego Romero's script on crab-eating macaque's: https://bitbucket.org/ee_reh_neh/tomoki_rna/src/507726b71d0c20ed3a202b611a770ce81575ed30/expression_plots_and_variance.r?at=master&fileviewer=file-view-default


#############################################################################################
### 0. Prepare Workspace ###
#############################################################################################

library(Rsubread)
library(RColorBrewer)
library(edgeR)
library(Homo.sapiens)
library(limma)
library(tximport)
library(tximportData)
library(ensembldb)
library(EnsDb.Hsapiens.v75)
library(rhdf5)
library(readr)
library(AnnotationDbi)
library(org.Pf.plasmo.db)
library(vegan)

#############################################################################################
### 1. Check quality of reads ###
#############################################################################################

# Mapped reads from from SRA Run selector
readSummary=read.table("/home/users/ntu/kbobowik/Sumba/Human_BWA_AlignmentStats/humanAlignment_all_YamagishiReads_PysamSummaryStatistics.txt", as.is=T, header=F, sep=" ")
colnames(readSummary) <- c("SRA_Accession_Number", "Total", "Mapped", "Unmapped", "QC_Fail", "Mapping_quality_lessThan_30")

# Organise df by accession number
readSummary=readSummary[order(readSummary$SRA_Accession_Number),]
pdf("readSummary_MappedvsUnmapped_BWA.pdf", height=10, width=30)
barplot(t(readSummary[,3:4])*1e-6, col=c("blue","red"), ylim=c(0, max(readSummary$Total)*1e-6+20), names.arg=readSummary$SRA_Accession_Number, las=3, cex.names=0.75, ylab="Total Mapped/Unmapped")
mtext("Reads", side=1, line=5)
legend("topright", legend=c("Mapped","Unmapped"), fill=c("blue", "red"))
dev.off()

#############################################################################################
### 2. Read in Human count data and organise sample information ###
#############################################################################################

# Reading in human count data with Tximport, using Kallisto-alligned files. Ensemble human annotation based on genome_build: GRCh37
edb <- EnsDb.Hsapiens.v75
Tx.ensemble <- transcripts(edb, columns = c("tx_id", "gene_id", "gene_name"), return.type = "DataFrame")
tx2gene.hg38 <- Tx.ensemble[,c(1,2)]
dir.hg38="/home/users/ntu/kbobowik/scratch/Kallisto/Human_PF_Combined"
samples.hg38=list.files(dir.hg38, pattern="DRR")
hg38.files <- file.path(dir.hg38, samples.hg38, "noDecimal_abundance.tsv")
names(hg38.files) <- samples.hg38
all(file.exists(hg38.files))
tx.hg38 <- tximport(hg38.files, type = "kallisto", tx2gene = tx2gene.hg38, reader = read_tsv)
dge_hg38 <- DGEList(tx.hg38$counts)

# Organising gene annotation using bioconductor's Homo.sapiens package, version 1.3.1 (Based on genome:  hg19)
geneid <- rownames(dge_hg38)
genes <- select(Homo.sapiens, keys=geneid, columns=c("SYMBOL", "TXCHROM"), keytype="ENSEMBL")

# Check for and remove duplicated gene IDs, then add genes dataframe to DGEList object
genes <- genes[!duplicated(genes$ENSEMBL),]
dge_hg38$genes <- genes

# Visualise library size
cols <- brewer.pal(12,"Set3")
pdf("kallisto_Total_n_Reads_ENSEMBLGenes.pdf", height=10, width=30)
barplot(dge_hg38$samples$lib.size*1e-6, ylab="Library size (millions)", cex.names=0.75, col=cols, names=colnames(dge_hg38), las=3, ylim=c(0,max(dge_hg38$samples$lib.size*1e-6)+10))
dev.off()

# Total number of genes
pdf("kallisto_Total_n_Genes_ENSEMBL.pdf", height=10, width=30)
barplot(apply(dge_hg38$counts, 2, function(c)sum(c!=0)), ylab="n Genes", cex.names=0.75, col=cols, names=colnames(dge_hg38), las=3)
dev.off()

#############################################################################################
### 3. Read in PF count data ###
#############################################################################################

dir.pf="/home/users/ntu/kbobowik/scratch/Kallisto/Human_PF_Combined"
samples.pf=list.files(dir.pf, pattern="DRR")
pf.files <- file.path(dir.pf, samples.pf, "noDecimal_abundance.tsv")
names(pf.files)<- samples.pf
all(file.exists(pf.files))
# readDGE(pf.files, columns=c(1,3))
pf <- org.Pf.plasmo.db
k <- keys(pf, keytype="SYMBOL")
tx2gene.pf <- select(pf, keys = k, keytype = "SYMBOL", columns = c("ORF", "SYMBOL"))
tx.pf <- tximport(pf.files, type = "kallisto", tx2gene =tx2gene.pf, reader = read_tsv)
dge_pf <- DGEList(tx.pf$counts)
geneid <- rownames(dge_pf)
genes <- select(org.Pf.plasmo.db, keys=geneid, columns=c("ORF", "GENENAME"), keytype="SYMBOL")
genes <- genes[!duplicated(genes$SYMBOL),]
dge_pf$genes <- genes
# dge_hg38=dge_hg38[,which(dge_hg38$samples$lib.size >= 10000000)]
cpm.pf <- cpm(dge_pf) 
lcpm.pf <- cpm(dge_pf, log=TRUE)

# Remove genes that are lowly expressed- a gene is only retained if it is expressed at log-CPM > 1 in at least half of the libraries
keep.exprs.pf <- rowSums(lcpm.pf>1) >= (nrow(dge_pf$samples)*0.5)
dge_pf <- dge_pf[keep.exprs.pf,, keep.lib.sizes=FALSE]

```
# Duplicates in PF file
a=read.delim(kallisto.files[1])
idx <- duplicated(a[,1]) | duplicated(a[,1], fromLast = TRUE)
x=a[idx,]
n_occur <- data.frame(table(a[,1]))
n=a[a[,1] %in% n_occur$Var1[n_occur$Freq > 1],]
x <- [!duplicated(genes$ENSEMBL),]
```

#############################################################################################
### 4. Data Pre-processing ###
#############################################################################################

# Load in Yamagishi et al 2014 supplementary table 11 (with patient information) and SRA run table
sup11=read.delim("/home/users/ntu/kbobowik/Sumba/Yamagishi2014_Supplemental_Table_11.txt", as.is=T, strip.white=T)
sra=read.delim("/home/users/ntu/kbobowik/Sumba/SraRunTable.txt", as.is=T, strip.white=T)

# Create empty matrix and fill in with patient names, Age, Sex, and SRA Run ID
sample_summary=matrix("NA",nrow=nrow(sup11),ncol=5)
sample_summary[,1:3]=c(sup11$Patient_Name,sup11$Age, sup11$Sex)
colnames(sample_summary)=c("Patient_Name", "Age", "Sex", "SRA_ID", "sample_ID")

# Discrepancy in SRA table and Supplementary 11 table (i.e., one says "malaria7#09" and one says "malaria7#009"
sample_summary[103,1]="malaria7#009"

# Match patient names with the patient names in the SRA table, then grab the corresponding SRA run IDs
sample_summary[,4]=sra[match(sample_summary[,1], sra[,4]),7]
sample_summary[,5]=paste("aligned_human_",sample_summary[,4],".sam", sep="")

# Assign gender and age
dge_hg38$samples$group <- as.factor(sample_summary[match(colnames(dge_hg38), sample_summary[,4]),3])
dge_hg38$samples$age <- as.factor(sample_summary[match(colnames(dge_hg38), sample_summary[,4]),2])

# Assign number of reads mapped to PF genome
#PF.reads <- dge_pf$samples$lib.size
dge_hg38$samples$PF.lib.size <- dge_pf$samples$lib.size

# Filter out samples with library size <10 million and samples with no gender information and assign covariate names
dge_hg38=dge_hg38[,which(dge_hg38$samples$lib.size >= 10000000)]
dge_hg38=dge_hg38[,which(dge_hg38$samples$group != "NA")]
age <- dge_hg38$samples$age
group <- dge_hg38$samples$group
PF.lib.size <- dge_hg38$samples$PF.lib.size

# Visualise library size after filtering
cols <- brewer.pal(12,"Set3")
pdf("kallisto_afterFilter_Total_n_Reads_ENSEMBLGenes.pdf", height=10, width=30)
barplot(dge_hg38$samples$lib.size*1e-6, ylab="Library size (millions)", cex.names=0.75, col=cols, names=colnames(dge_hg38), las=3, ylim=c(0,max(dge_hg38$samples$lib.size*1e-6)+10))
dev.off()

pdf("kallisto_afterFilter_Total_n_Genes_ENSEMBL.pdf", height=10, width=30)
barplot(apply(dge_hg38$counts, 2, function(c)sum(c!=0)), ylab="n Genes", cex.names=0.75, col=cols, names=colnames(dge_hg38), las=3)
dev.off()

# Transform from the raw scale to CPM and log-CPM values
# Prior count for logCPM = 0.25
cpm <- cpm(dge_hg38) 
lcpm <- cpm(dge_hg38, log=TRUE)

# Remove genes that are lowly expressed- a gene is only retained if it is expressed at log-CPM > 1 in at least half of the libraries
keep.exprs <- rowSums(lcpm>1) >= (nrow(dge_hg38$samples)*0.5)
dge_hg38 <- dge_hg38[keep.exprs,, keep.lib.sizes=FALSE]

# Compare library sizes before and after removing lowly-expressed genes
nsamples <- ncol(dge_hg38)
col <- colorRampPalette(brewer.pal(nsamples, "Paired")) (nsamples)
pdf("kallisto_librarySize_beforeANDafterLowlyExpressedGenes_ENSEMBL.pdf")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,max(density(lcpm)$x)+.2), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
}

lcpm <- cpm(dge_hg38, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.31), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
}
dev.off()

# Normalise gene expression distributions (i.e., no bias introduced during sample preparation/sequencing)
dge_hg38 <- calcNormFactors(dge_hg38, method = "TMM")

# Duplicate data, set normalisation back to 1, and plot difference between normalised and non-normalised data
y2 <- dge_hg38
y2$samples$norm.factors <- 1
pdf("kallisto_TMM_normalisedVSnonnormalised_ENSEMBL.pdf", height=10, width=30)
lcpm <- cpm(y2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
y2 <- calcNormFactors(y2)
lcpm <- cpm(y2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")
dev.off()

# MDS
lcpm <- cpm(dge_hg38, log=TRUE)
col.group <- group
levels(col.group) <- c("#FF0000", "#00FF00")
col.group <- as.character(col.group)
col.age=cut(as.numeric(as.character(age)), c(0,15,30,70))
levels(col.age) <- brewer.pal(nlevels(col.age),"Set1")
col.age <- as.character(col.age)
pdf("kallisto_MDS_F1000Pipeline_ENSEMBL_GenderAge.pdf")
plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Sample Gender")
plotMDS(lcpm, labels=age, col=col.age)
title(main="A. Sample Age")
dev.off()

# PCA plotting function
plot.pca <- function(dataToPca, speciesCol, namesPch, sampleNames){
    pca <- prcomp(t(dataToPca), scale=T, center=T)
    pca.var <- pca$sdev^2/sum(pca$sdev^2)
    plot(pca$x[,1], pca$x[,2], col=speciesCol, pch=namesPch, cex=2, xlab=paste("PC1 (", round(pca.var[1]*100, digits=2), "% of variance)", sep=""), ylab=paste("PC2 (", round(pca.var[2]*100, digits=2), "% of variance)", sep=""), main="PCA")
    plot(pca$x[,2], pca$x[,3], col=speciesCol, pch=namesPch, cex=2, xlab=paste("PC2 (", round(pca.var[2]*100, digits=2), "% of variance)", sep=""), ylab=paste("PC3 (", round(pca.var[3]*100, digits=2), "% of variance)", sep=""))
    plot(pca$x[,3], pca$x[,4], col=speciesCol, pch=namesPch, cex=2, xlab=paste("PC3 (", round(pca.var[3]*100, digits=2), "% of variance)", sep=""), ylab=paste("PC4 (", round(pca.var[4]*100, digits=2), "% of variance)", sep=""))

    return(pca)
}

# PCA association function
pc.assoc <- function(pca.data){
    all.pcs <- data.frame()
    for (i in 1:ncol(pca.data$x)){
        all.assoc <- vector()
        for (j in 1:ncol(all.covars.df)){
            test.assoc <- anova(lm(pca.data$x[,i] ~ all.covars.df[,j]))[1,5]
            all.assoc <- c(all.assoc, test.assoc)
        }
        single.pc <- c(i, all.assoc)
        all.pcs <- rbind(all.pcs, single.pc)
    }
    names(all.pcs) <- c("PC", colnames(all.covars.df))

    print ("Here are the relationships between PCs and some possible covariates")
    print (all.pcs)
    return (all.pcs)
}

#Prepare covariate matrix
all.covars.df <- dge_hg38$samples[,c(1,2,4,5)] 
all.covars.df$group <- factor(all.covars.df$group)
all.covars.df$age <- factor(all.covars.df$age)

# Plot PCA
pdf(file="pca_clean.pdf")
pcaresults <- plot.pca(lcpm, col.group, 19, colnames(dge_hg38))
dev.off()

all.pcs <- pc.assoc(pcaresults)
all.pcs

# Write out the covariates:
write.table(all.pcs, file="pca_covariates.txt", col.names=T, row.names=F, quote=F, sep="\t")


#############################################################################################
### 5. Differential expression analysis ###
#############################################################################################

# Make a design matrix. The whole purpose of teh deisgn matrix is to set up your parameters for estimating the relationships among your  variables (linear regression)
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
contr.matrix <- makeContrasts(MvsF=M-F, levels=colnames(design))

# Remove heteroscedascity from count data
pdf("kallisto_Voom_F1000Pipeline_ENSEMBL.pdf")
v <- voom(dge_hg38, design, plot=TRUE)
dev.off()

# Fit linear models of gender data to human DGE data
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
pdf("kallisto_VoomWeights_F1000Pipeline_ENSEMBL.pdf")
plotSA(efit)
dev.off()

# Summarise the number of significantly up- and down-regulated genes
topTable(efit)
dt <- decideTests(efit)
summary(dt)


# Examine top DE genes
M.vs.F <- topTreat(efit, coef=1, n=Inf)

pdf("kallisto_PlotMD_DeGenes_F1000Pipeline_ENSEMBL.pdf")
plotMD(efit, column=1, status=dt[,1], main=colnames(efit), xlim=c(-8,13))
dev.off()
