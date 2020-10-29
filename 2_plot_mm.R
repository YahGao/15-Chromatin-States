#========#
# Fig. 1 #
#========#

# Fig. 1b
All_peaks <- read.table("All_peaks.bed", header = T, sep = "\t")
unique(All_peaks$Sample)

library(GenomicRanges)
source("Bed2GRange.R")
library(corrplot)

All_peaks_1 <- All_peaks
All_peaks_1$Chromosome <- gsub("chr","",All_peaks_1$Chromosome)
All_peaks_autosome <- subset(All_peaks_1,Chromosome%in%paste(1:29))
Smaples <- unique(All_peaks_autosome$Sample)

All_peaks_autosome_GR <- makeGRangesFromDataFrame(All_peaks_autosome)
All_peaks_autosome_merged <- reduce(All_peaks_autosome_GR)
All_peaks_autosome_merged_DF <- Grange2df(All_peaks_autosome_merged)

Overlaps_hits <- findOverlaps(All_peaks_autosome_GR,All_peaks_autosome_merged,select = "all")
Overlaps_hits_DF <- as.data.frame(Overlaps_hits)
Overlaps_hits_DF$Samples <- All_peaks_autosome$Sample
All_Samples <- unique(Overlaps_hits_DF$Samples)

MyCorrelations <- array(NA,dim = c(length(All_Samples),length(All_Samples)))
colnames(MyCorrelations) <- All_Samples
row.names(MyCorrelations) <- All_Samples

for(i in 1:length(All_Samples)){
  Sample_A <- subset(All_peaks_autosome,Sample==All_Samples[i])
  Sample_A_GR <- makeGRangesFromDataFrame(Sample_A)
  A_length <- dim(Sample_A)[1]
  for(j in 1:length(All_Samples)){
    Sample_B <- subset(All_peaks_autosome,Sample==All_Samples[j])
    Sample_B_GR <- makeGRangesFromDataFrame(Sample_B)
    B_length <- dim(Sample_B)[1]
    Overlaps_hits_A_B <- findOverlaps(Sample_A_GR,Sample_B_GR,select = "all")
    head(Overlaps_hits_A_B)
    Overlaps_hits_A_B_DF <- as.data.frame(Overlaps_hits_A_B)    
    A_B <- length(unique(Overlaps_hits_A_B_DF$queryHits))/A_length
    B_A <- length(unique(Overlaps_hits_A_B_DF$subjectHits))/B_length  
    MyCorrelations[i,j] <- A_B*100
    MyCorrelations[j,i] <- B_A*100
  }
}

tiff(file = "Fig. 1b.tiff", res = 300, width = 2400, height = 2400,compression = "lzw")
col3 <- colorRampPalette(c("red", "lightyellow", "royalblue")) 
corrplot(MyCorrelations/100,method = "color",mar = c(0, 1, 0, 2),col = col3(200), cl.lim = c(0, 1),is.corr = FALSE,tl.col="black",order = "hclust")
dev.off()


# Fig. 1c
aw1 <- read.csv("AW1.tsv", header = T, sep = "\t")
aw2 <- read.csv("AW2.tsv", header = T, sep = "\t")
aw3 <- read.csv("AW3.tsv", header = T, sep = "\t")
bw1 <- read.csv("BW1.tsv", header = T, sep = "\t")
bw2 <- read.csv("BW2.tsv", header = T, sep = "\t")
bw3 <- read.csv("BW3.tsv", header = T, sep = "\t")

aw1 <- data.frame(aw1$Gene.ID, aw1$TPM)
colnames(aw1) <- c("Gene", "AW1")
aw2 <- data.frame(aw2$Gene.ID, aw2$TPM)
colnames(aw2) <- c("Gene", "AW2")
aw2$Ind <- "AW2"
aw3 <- data.frame(aw3$Gene.ID, aw3$TPM)
colnames(aw3) <- c("Gene", "AW3")
aw3$Ind <- "AW3"
bw1 <- data.frame(bw1$Gene.ID, bw1$TPM)
colnames(bw1) <- c("Gene", "BW1")
bw1$Ind <- "BW1"
bw2 <- data.frame(bw2$Gene.ID, bw2$TPM)
colnames(bw2) <- c("Gene", "BW2")
bw2$Ind <- "BW2"
bw3 <- data.frame(bw3$Gene.ID, bw3$TPM)
colnames(bw3) <- c("Gene", "BW3")
bw3$Ind <- "BW3"

tpm_gene <- merge(aw1, aw2)
tpm_gene <- tpm_gene[,-4]
tpm_gene <- merge(tpm_gene, aw3)
tpm_gene <- tpm_gene[,-5]
tpm_gene <- merge(tpm_gene, bw1)
tpm_gene <- tpm_gene[,-6]
tpm_gene <- merge(tpm_gene, bw2)
tpm_gene <- tpm_gene[,-7]
tpm_gene <- merge(tpm_gene, bw3)
tpm_gene <- tpm_gene[,-8]
rownames(tpm_gene) <- tpm_gene$Gene
tpm_gene <- tpm_gene[,-1]

aw1 <- tpm_gene$AW1
aw1 <- as.data.frame(aw1)
aw1$Ind <- "AW1"
colnames(aw1) <- c("TPM", "Ind")
aw2 <- tpm_gene$AW2
aw2 <- as.data.frame(aw2)
aw2$Ind <- "AW2"
colnames(aw2) <- c("TPM", "Ind")
aw3 <- tpm_gene$AW3
aw3 <- as.data.frame(aw3)
aw3$Ind <- "AW3"
colnames(aw3) <- c("TPM", "Ind")
bw1 <- tpm_gene$BW1
bw1 <- as.data.frame(bw1)
bw1$Ind <- "BW1"
colnames(bw1) <- c("TPM", "Ind")
bw2 <- tpm_gene$BW2
bw2 <- as.data.frame(bw2)
bw2$Ind <- "BW2"
colnames(bw2) <- c("TPM", "Ind")
bw3 <- tpm_gene$BW2
bw3 <- as.data.frame(bw3)
bw3$Ind <- "BW3"
colnames(bw3) <- c("TPM", "Ind")

color <- c("#386CB0", "#F0027F", "#66A61E", "#1B9E77", "#FDC086", "#0F27F0", "#386CB0", "#F0027F", "#66A61E", "#1B9E77", "#FDC086", "#0F27F0")

tpm_gene2 <- cor(tpm_gene, method = "spearman")
tpm_gene2 <- 1 - tpm_gene2
tpm_gene2 <- as.matrix(tpm_gene2)

library(ggplot2)
library("MASS")
library(RColorBrewer)
tpm_gene2 <- isoMDS(tpm_gene2,k=5)

pc1 <- tpm_gene2$points[,1]
pc2 <- tpm_gene2$points[,2]
pc3 <- tpm_gene2$points[,3]
pc4 <- tpm_gene2$points[,4]
pc5 <- tpm_gene2$points[,5]

pc1_p <- sum((tpm_gene2$points[,1])^2)/sum((tpm_gene2$points)^2)
pc1_p <- pc1_p*100
pc1_p <- paste("Coordinate 1 (",sprintf("%.2f",pc1_p),"%",")",sep="")
pc2_p <- sum((tpm_gene2$points[,2])^2)/sum((tpm_gene2$points)^2)
pc2_p <- pc2_p*100
pc2_p <- paste("Coordinate 2 (",sprintf("%.2f",pc2_p),"%",")",sep="")
pc3_p <- sum((tpm_gene2$points[,3])^2)/sum((tpm_gene2$points)^2)
pc3_p <- pc3_p*100
pc3_p <- paste("Coordinate 3 (",sprintf("%.2f",pc3_p),"%",")",sep="")
pc4_p <- sum((tpm_gene2$points[,4])^2)/sum((tpm_gene2$points)^2)
pc4_p <- pc4_p*100
pc4_p <- paste("Coordinate 4 (",sprintf("%.2f",pc4_p),"%",")",sep="")
pc5_p <- sum((tpm_gene2$points[,5])^2)/sum((tpm_gene2$points)^2)
pc5_p <- pc5_p*100
pc5_p <- paste("Coordinate 5 (",sprintf("%.2f",pc5_p),"%",")",sep="")

all <- data.frame(pc1, pc2, pc3, pc4, pc5)
all1 <- data.frame(all, color)
all_p <- data.frame(pc1_p, pc2_p, pc3_p, pc4_p, pc5_p)


tiff(file="Fig. 1c.tiff", width = 5, height = 6, units = "cm", res = 600, pointsize = 8,compression= "lzw")
plot(all$pc1, all$pc2, type = "n", pch =19, col=all1$color, xlab = all_p$pc1_p, 
     ylab = all_p$pc2_p, xlim = c(-0.025, 0.035))
text(all[1,1], all[1,2], "AW1", col = "#386CB0", cex = .8)
text(all[2,1], all[2,2], "AW2", col = "#F0027F", cex = .8)
text(all[3,1], all[3,2], "AW3", col = "#66A61E", cex = .8)
text(all[4,1], all[4,2], "BW1", col = "#1B9E77", cex = .8)
text(all[5,1], all[5,2], "BW2", col = "#FDC086", cex = .8)
text(all[6,1], all[6,2], "BW3", col = "#0F27F0", cex = .8)
dev.off()



#========#
# Fig. 2 #
#========#
# Fig. 2c 
library(gplots)
library(RColorBrewer)
library(data.table)

emissions <- read.table("emissions_15.txt",header = T,stringsAsFactors = F,sep="\t")
emissions <- emissions[,c("H3K4me3","H3K4me1","H3K27Ac","ATAC","CTCF","H3K27me3")]
emissions_order <- emissions[c(8,10,9,7,5,13,6,4,12,14,11,15,3,2,1),]
my_palette <- colorRampPalette(c("white","darkblue"))(n = 199)
col_breaks = c(seq(0,0.5,length=100), seq(0.51,1,length=100))

tiff(file="Fig. 2c .tiff", width = 6, height = 16.7, units = "cm", res = 600, pointsize = 8,compression= "lzw")
heatmap.2(as.matrix(emissions_order), notecol="black", density.info="none", trace="none", margins =c(7,2), key.par=list(mar=c(3.5,1,5,0)), col=my_palette,
          breaks=col_breaks, dendrogram="none", labRow = NA, labCol = colnames(emissions), cexRow=1.5, cexCol=1.5, key=T, keysize=1.2,key.title=NA, Colv=F, Rowv=F,
          symm=T,symkey=F,symbreaks=T, scale="none",srtCol=90, adjCol = c(0.95,0), key.xlab=expression(paste(italic(Emissions),sep=)))
dev.off()


# Fig. 2d
Overlaps_ATAC <- read.table("BW_15_overlap.txt",header = T,stringsAsFactors = F,sep="\t")
Overlaps_ATAC_order <- Overlaps_ATAC[c(8,10,9,7,5,13,6,4,12,14,11,15,3,2,1),]
colnames(Overlaps_ATAC_order) <- gsub(".bed.gz","",colnames(Overlaps_ATAC_order))
Genome_coverage <- Overlaps_ATAC_order[,c("state..Emission.order.","Genome..")]
my_palette <- colorRampPalette(c("white","darkblue"))(n = 199)
col_breaks = c(seq(0,2,length=100), seq(2.1,4,length=100))

tiff(file="Fig. 2d.tiff", width = 6, height = 16.7, units = "cm", res = 600, pointsize = 8,compression= "lzw")
heatmap.2(as.matrix(Genome_coverage[,c(2,2)]), notecol="black", density.info="none", trace="none", margins =c(7,2), key.par=list(mar=c(3.5,1,5,0)), col=my_palette,
          breaks=col_breaks, dendrogram="none", labRow = NA, labCol = NA, cexRow=1.5, cexCol=1.5, key=F, keysize=1.2, key.title=NA, Colv=F, Rowv=F,
          symm=T,symkey=F,symbreaks=T, scale="none",srtCol=90, adjCol = c(0.95,0), key.xlab=expression(paste(italic(Emissions),sep=)))
dev.off()


# Fig. 2e
Overlaps_ATAC_order_gene <- Overlaps_ATAC_order[,c("state..Emission.order.","G_CpG_island","G_Genic",
                                                   "G_CDS","G_Introns","G_FiveUTR","G_Promoter","G_TSS",
                                                   "G_TES","G_ThreeUTR","E_Expr","E_Expr_TSS","E_Expr_TES",
                                                   "E_Repr","E_Repr_TSS","E_Repr_TES","G_ZNF","C_TF","REPC_sg",
                                                   "R_LINE","R_SINE","R_Satellite","R_Simple_repeat","R_LTR")]
Mynames <- c("CpG_island","Genes","CDS","Introns","5'UTR","Promoter","TSS","TES",
             "3'UTR","Expr_genes","Expr_TSS","Expr_TES","Repr_genes","Repr_TSS","Repr_TES",
             "ZNF_genes","TF","REPC_SG","LINE","SINE","Satellite","Simple_repeat","LTR")
my_palette <- colorRampPalette(c("white","red"))(n = 199)
col_breaks = c(seq(0,5,length=100), seq(5.1,10,length=100))

tiff(file="Fig. 2e.tiff", width = 24, height = 16.7, units = "cm", res = 600, pointsize = 8,compression= "lzw")
heatmap.2(as.matrix(Overlaps_ATAC_order_gene[,-c(1)]), cellnote = round(as.matrix(Overlaps_ATAC_order_gene[,-c(1)]),digits = 2), notecol="black", density.info="none", trace="none",
          margins =c(10,2), key.par=list(mar=c(3.5,1,5,0)), col=my_palette, breaks=col_breaks, dendrogram="none", labRow = NA, labCol = Mynames, cexRow=1.5,cexCol=1.5,
          key=F, keysize=0.3,key.title=NA,Colv=F,Rowv=F, symm=T,symkey=F,symbreaks=T, scale="none",srtCol=90, adjCol = c(0.95,0), key.xlab=expression(paste(italic(Overlap),sep=""))) 
dev.off()


# Fig. 2f
emision <- read.table("emissions_15_new.txt", header = T, sep = "\t")
colnames(emision) <- c("State", "New")
Mycolors <- c("#ff0000", "#ff5050", "#ff9999", "#70ad47", "#7f6000", "#bf9000", "#ffd966", "#ffe699","#fff2cc", "#f4b183", "#c55a11", "#8497b0", "#a9d18e", "#767171", "#ffffff")

tissueb <- read.table("BW_final2.meth", header = F, sep = "\t")
colnames(tissueb) <- c("State", "Me")
tissueb <- merge(emision, tissueb)

tiff(file = "Fig. 2f.tif", res = 300, width = 1100, height = 1200,compression = "lzw")
par(mar=c(4,10,2,2))
boxplot(Me ~ New, tissueb, col = Mycolors, boxwex = 0.8, pch = 19, xaxt = "n", yaxt = "n", cex = 0.5, xlab = "", ylab = "", las = 2, outline = F, varwidth = F,
        cex.axis = 0.6, cex.lab = 0.8,horizontal = T, notch = F, medlwd = 1, range=1)
axis(3,seq(0, 1, 0.2), las = 1, cex.axis = 0.6)
dev.off()


# Fig. 2g
Overlaps_ATAC1 <- read.table("BW_15_overlap.txt",header = T,stringsAsFactors = F,sep="\t")
Overlaps_ATAC_order1 <- Overlaps_ATAC1[c(8,10,9,7,5,13,6,4,12,14,11,15,3,2,1),]
colnames(Overlaps_ATAC_order1) <- gsub(".bed.gz","",colnames(Overlaps_ATAC_order1))
Overlaps_ATAC2 <- read.table("AW_15_overlap.txt",header = T,stringsAsFactors = F,sep="\t")
Overlaps_ATAC_order2 <- Overlaps_ATAC2[c(8,10,9,7,5,13,6,4,12,14,11,15,3,2,1),]
colnames(Overlaps_ATAC_order2) <- gsub(".bed.gz","",colnames(Overlaps_ATAC_order2))
Overlaps_ATAC3 <- read.table("CO_15_overlap.txt",header = T,stringsAsFactors = F,sep="\t")
Overlaps_ATAC_order3 <- Overlaps_ATAC3[c(8,10,9,7,5,13,6,4,12,14,11,15,3,2,1),]
colnames(Overlaps_ATAC_order3) <- gsub(".bed.gz","",colnames(Overlaps_ATAC_order3))
Overlaps_ATAC4 <- read.table("BT_15_overlap.txt",header = T,stringsAsFactors = F,sep="\t")
Overlaps_ATAC_order4 <- Overlaps_ATAC4[c(8,10,9,7,5,13,6,4,12,14,11,15,3,2,1),]
colnames(Overlaps_ATAC_order4) <- gsub(".bed.gz","",colnames(Overlaps_ATAC_order4))

Convserved_regions1 <- Overlaps_ATAC_order1[,"C_Conserved"]
Convserved_regions2 <- Overlaps_ATAC_order2[,"C_Conserved"]
Convserved_regions3 <- Overlaps_ATAC_order3[,"C_Conserved"]
Convserved_regions4 <- Overlaps_ATAC_order4[,"C_Conserved"]

Convserved_regions_all <- cbind(Convserved_regions1, Convserved_regions2, Convserved_regions3, Convserved_regions4)
Convserved_regions_all <- as.data.frame(Convserved_regions_all)
Convserved_regions_all$mean <- apply(Convserved_regions_all, 1, mean)
Convserved_regions_all$sd <- apply(Convserved_regions_all[,1:4], 1, sd)
Convserved_regions_all$state <- 1:15

Mycolors <- c("white","#767171","#a9d18e","#8497b0","#c54611","#f4b183","#fff2cc","#ffe699","#ffd966","#bf9000","#7f6000","#70ad47","#ff9999","#ff5050","#ff0000")

library(dplyr)
library(forcats)
library(ggplot2)
Convserved_regions_all$state <- as.factor(Convserved_regions_all$state)

tiff(file = "Fig. 2g.tif", res = 300, width = 1000, height = 1850,compression = "lzw")
par(mar=c(4,4,2,2))
Convserved_regions_all %>%
  mutate(state = fct_reorder(state, desc(state))) %>%
  ggplot(aes(x=state, y=mean)) +
  geom_bar(stat="identity", fill = rev(Mycolors), alpha=.8, width=1, color = "black") +
  geom_errorbar(aes(x=state, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="grey", alpha=0.9, size=1) + 
  coord_flip() + ylim(0,5) +
  theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("") + ylab(paste("Enrichment for non-coding\n", "mammalian-conserved elements\n (GERP)",sep=""))
dev.off()


# Fig. 2h & i
TSS <- read.table("BW_15_bw_exp_TSS_neighborhood.txt",header = T,stringsAsFactors = F,sep="\t")
TSS <- TSS[,-1]
TSS <- TSS[c(8,10,9,7,5,13,6,4,12,14,11,15,3,2,1),]
State_names <- c("15 Quies", "14 ReprPC","13 BivFlnk","12 ReprWkCTCF",
                 "11 ATAC","10 EnhWkCTCFATAC","9 EnhPoisATAC","8 EnhPois","7 EnhWk",
                 "6 EnhAATAC","5 EnhA", "4 TxFlnk","3 TssAFlnk", "2 TssAATACCTCF","1 TssA")
row.names(TSS) <- rev(State_names)
my_palette <- colorRampPalette(c("white","red"))(n = 199)
col_breaks = c(seq(0,30,length=100), seq(30.1,60,length=100))
tiff(file = "Fig. 2h1.tiff", res = 300, width = 1800, height = 1200,compression = "lzw")
heatmap.2(as.matrix(TSS), notecol="black", density.info="none", trace="none", margins =c(4,12), key.par=list(mar=c(2,0.5,0,0)), col=my_palette, breaks=col_breaks, dendrogram="none",    
          labRow = rev(State_names), labCol = seq(-2000,2000,by = 200), cexRow=1.2,cexCol=1.2, key=T, keysize=0.6,key.title=NA,Colv=F,Rowv=F, symm=T,symkey=F,symbreaks=T, 
          scale="none",srtCol=90, adjCol = c(0.95,0), key.xlab=expression(paste(italic(EnriFold),sep="")))          
dev.off()

TES <- read.table("BW_15_bw_exp_TES_neighborhood.txt",header = T,stringsAsFactors = F,sep="\t")
TES <- TES[,-1]
TES <- TES[c(8,10,9,7,5,13,6,4,12,14,11,15,3,2,1),]
row.names(TES) <- rev(State_names)
my_palette <- colorRampPalette(c("white","red"))(n = 199)
col_breaks = c(seq(0,30,length=100), seq(30.1,60,length=100))
tiff(file = "Fig. 2h2.tiff",res = 300, width = 1800, height = 1200,compression = "lzw")
heatmap.2(as.matrix(TES), notecol="black", density.info="none", trace="none", margins =c(4,12), key.par=list(mar=c(2,0.5,0,0)), col=my_palette, breaks=col_breaks, dendrogram="none",
          labRow = rev(State_names), labCol = seq(-2000,2000,by = 200), cexRow=1.2,cexCol=1.2, key=T, keysize=0.6,key.title=NA,Colv=F,Rowv=F, symm=T,symkey=F,symbreaks=T, 
          scale="none",srtCol=90, adjCol = c(0.95,0), key.xlab=expression(paste(italic(EnriFold),sep="")))
dev.off()


TSS <- read.table("BW_15_bw_repr_TSS_neighborhood.txt",header = T,stringsAsFactors = F,sep="\t")
TSS <- TSS[,-1]
TES <- TES[c(8,10,9,7,5,13,6,4,12,14,11,15,3,2,1),]
row.names(TSS) <- rev(State_names)
my_palette <- colorRampPalette(c("white","red"))(n = 199)
col_breaks = c(seq(0,30,length=100),seq(30.1,60,length=100))
tiff(file = "Fig. 2i1.tiff",res = 300, width = 1800, height = 1200,compression = "lzw")
heatmap.2(as.matrix(TSS), notecol="black", density.info="none", trace="none", margins =c(4,12), key.par=list(mar=c(2,0.5,0,0)), col=my_palette, breaks=col_breaks, dendrogram="none",
          labRow = rev(State_names), labCol = seq(-2000,2000,by = 200), cexRow=1.2,cexCol=1.2, key=T, keysize=0.6,key.title=NA,Colv=F,Rowv=F, symm=T,symkey=F,symbreaks=T, 
          scale="none",srtCol=90, adjCol = c(0.95,0), key.xlab=expression(paste(italic(EnriFold),sep="")))
dev.off()


TES <- read.table("BW_b_15_bw_repr_TES_neighborhood.txt",header = T,stringsAsFactors = F,sep="\t")
TES <- TES[,-1]
TES <- TES[c(8,10,9,7,5,13,6,4,12,14,11,15,3,2,1),]
row.names(TES) <- rev(State_names)
my_palette <- colorRampPalette(c("white","red"))(n = 199)
col_breaks = c(seq(0,30,length=100), seq(30.1,60,length=100))
tiff(file = "Fig. 2i2.tiff",res = 300, width = 1800, height = 1200,compression = "lzw")
heatmap.2(as.matrix(TES), notecol="black", density.info="none", trace="none", margins =c(4,12), key.par=list(mar=c(2,0.5,0,0)), col=my_palette, breaks=col_breaks, dendrogram="none",
          labRow = rev(State_names), labCol = seq(-2000,2000,by = 200), cexRow=1.2,cexCol=1.2, key=T, keysize=0.6,key.title=NA,Colv=F,Rowv=F, symm=T,symkey=F,symbreaks=T, 
          scale="none",srtCol=90, adjCol = c(0.95,0), key.xlab=expression(paste(italic(EnriFold),sep="")))
dev.off()



#========#
# Fig. 3 #
#========#
# Fig. 3a
Mycolors <- c("white","#767171","#a9d18e","#8497b0","#c54611","#f4b183","#fff2cc","#ffe699","#ffd966","#bf9000","#7f6000","#70ad47","#ff9999","#ff5050","#ff0000")
QTLdb_regions1 <- Overlaps_ATAC_order1[,"C_QTLdb_bed1"]
QTLdb_regions2 <- Overlaps_ATAC_order2[,"C_QTLdb_bed1"]
QTLdb_regions3 <- Overlaps_ATAC_order3[,"C_QTLdb_bed1"]
QTLdb_regions4 <- Overlaps_ATAC_order4[,"C_QTLdb_bed1"]

QTLdb_regions_all <- cbind(QTLdb_regions1, QTLdb_regions2, QTLdb_regions3, QTLdb_regions4)
QTLdb_regions_all <- as.data.frame(QTLdb_regions_all)
QTLdb_regions_all$mean <- apply(QTLdb_regions_all, 1, mean)
QTLdb_regions_all$sd <- apply(QTLdb_regions_all[,1:4], 1, sd)
QTLdb_regions_all$state <- 1:15

library(dplyr)
library(forcats)
library(ggplot2)
QTLdb_regions_all$state <- as.factor(QTLdb_regions_all$state)

tiff(file = "Fig. 3a.tif", res = 300, width = 1000, height = 1850,compression = "lzw")
par(mar=c(4,4,2,2))
QTLdb_regions_all %>%
  mutate(state = fct_reorder(state, desc(state))) %>%
  ggplot(aes(x=state, y=mean)) +
  geom_bar(stat="identity", fill = rev(Mycolors), alpha=.8, width=1, color = "black") +
  geom_errorbar(aes(x=state, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="grey", alpha=0.9, size=1) + 
  coord_flip() + ylim(0,5) +
  theme(panel.grid.major = element_blank(),  panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("") + ylab(paste("Enrichment for 108274 cattle QTLs\n in 317 complex traits",sep=""))
dev.off()


# Fig. 3b
eqtl <- colnames(Overlaps_ATAC_order1)
cis_eqtl <- grep('cis_eqtl', eqtl, value = T)
for (i in cis_eqtl) {
  cis_eqtl_regions1 <- Overlaps_ATAC_order1[, i]
  cis_eqtl_regions2 <- Overlaps_ATAC_order2[, i]
  cis_eqtl_regions3 <- Overlaps_ATAC_order3[, i]
  cis_eqtl_regions4 <- Overlaps_ATAC_order4[, i]
  cis_eqtl_regions_all <- cbind(cis_eqtl_regions1, cis_eqtl_regions2, cis_eqtl_regions3, cis_eqtl_regions4)
  cis_eqtl_regions_all <- as.data.frame(cis_eqtl_regions_all)
  cis_eqtl_regions_all$mean <- apply(cis_eqtl_regions_all, 1, mean)
  cis_eqtl_regions_all$sd <- apply(cis_eqtl_regions_all[,1:4], 1, sd)
  cis_eqtl_regions_all$state <- 1:15
  cis_eqtl_regions_all$state <- as.factor(cis_eqtl_regions_all$state)
  
  tiff(file = paste(i,"Fig. 3b.tif", sep = ""),res = 300, width = 1000, height = 1850,compression = "lzw")
  par(mar=c(4,4,2,2))
  print(cis_eqtl_regions_all %>%
          mutate(state = fct_reorder(state, desc(state))) %>%
          ggplot(aes(x=state, y=mean)) +
          geom_bar(stat="identity", fill = rev(Mycolors), alpha=.8, width=1, color = "black") +
          geom_errorbar(aes(x=state, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="grey", alpha=0.9, size=1) + 
          coord_flip() + ylim(0,5) +
          theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
          xlab("") + ylab(paste("Enrichment for ",i,"\n in cattle",sep="")))
  dev.off()
}


# Fig. 3c
sqtl <- grep('sqtl', eqtl, value = T)
for (i in sqtl) {
  sqtl_regions1 <- Overlaps_ATAC_order1[, i]
  sqtl_regions2 <- Overlaps_ATAC_order2[, i]
  sqtl_regions3 <- Overlaps_ATAC_order3[, i]
  sqtl_regions4 <- Overlaps_ATAC_order4[, i]
  sqtl_regions_all <- cbind(sqtl_regions1, sqtl_regions2, sqtl_regions3, sqtl_regions4)
  sqtl_regions_all <- as.data.frame(sqtl_regions_all)
  sqtl_regions_all$mean <- apply(sqtl_regions_all, 1, mean)
  sqtl_regions_all$sd <- apply(sqtl_regions_all[,1:4], 1, sd)
  sqtl_regions_all$state <- 1:15
  sqtl_regions_all$state <- as.factor(sqtl_regions_all$state)
  
  tiff(file = paste(i,"Fig. 3c.tif", sep = ""),res = 300, width = 1000, height = 1850,compression = "lzw")
  par(mar=c(4,4,2,2))
  print(sqtl_regions_all %>%
          mutate(state = fct_reorder(state, desc(state))) %>%
          ggplot(aes(x=state, y=mean)) +
          geom_bar(stat="identity", fill = rev(Mycolors), alpha=.8, width=1, color = "black") +
          geom_errorbar(aes(x=state, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="grey", alpha=0.9, size=1) + 
          coord_flip() + ylim(0,5) +
          theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
          xlab("") + ylab(paste("Enrichment for ",i,"\n in cattle",sep="")))
  dev.off()
}


# Fig. 3d
ss <- colnames(Overlaps_ATAC_order1)
ss <- ss[3:9]
for (i in ss) {
  ss_regions1 <- Overlaps_ATAC_order1[, i]
  ss_regions2 <- Overlaps_ATAC_order2[, i]
  ss_regions3 <- Overlaps_ATAC_order3[, i]
  ss_regions4 <- Overlaps_ATAC_order4[, i]
  ss_regions_all <- cbind(ss_regions1, ss_regions2, ss_regions3, ss_regions4)
  ss_regions_all <- as.data.frame(ss_regions_all)
  ss_regions_all$mean <- apply(ss_regions_all, 1, mean)
  ss_regions_all$sd <- apply(ss_regions_all[,1:4], 1, sd)
  ss_regions_all$state <- 1:15
  ss_regions_all$state <- as.factor(ss_regions_all$state)
  
  tiff(file = paste(i,"Fig. 3d.tif", sep = ""),res = 300, width = 1000, height = 1850,compression = "lzw")
  par(mar=c(4,4,2,2))
  print(ss_regions_all %>%
          mutate(state = fct_reorder(state, desc(state))) %>%
          ggplot(aes(x=state, y=mean)) +
          geom_bar(stat="identity", fill = rev(Mycolors), alpha=.8, width=1, color = "black") +
          geom_errorbar(aes(x=state, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="grey", 
                        alpha=0.9, size=1) + 
          coord_flip() + ylim(0,5) +
          theme(panel.grid.major = element_blank(), 
                panel.background = element_blank(),
                axis.line = element_line(colour = "black")) +
          xlab("") + ylab(paste("Enrichment for ",i,"\n in cattle",sep="")))
  dev.off()
}


# Fig. 3e
library(gplots)

cis_eqtl1 <- gsub("_cis_eqtl","",cis_eqtl)
cis_eqtl_regions_all1 <- list()
for (i in cis_eqtl) {
  cis_eqtl_regions1 <- Overlaps_ATAC_order1[, i]
  cis_eqtl_regions2 <- Overlaps_ATAC_order2[, i]
  cis_eqtl_regions3 <- Overlaps_ATAC_order3[, i]
  cis_eqtl_regions4 <- Overlaps_ATAC_order4[, i]
  cis_eqtl_regions_all <- cbind(cis_eqtl_regions1, cis_eqtl_regions2, cis_eqtl_regions3, cis_eqtl_regions4)
  cis_eqtl_regions_all <- as.data.frame(cis_eqtl_regions_all)
  cis_eqtl_regions_all$mean <- apply(cis_eqtl_regions_all, 1, mean)
  cis_eqtl_regions_all1 <- cbind(cis_eqtl_regions_all1, cis_eqtl_regions_all$mean)
  cis_eqtl_regions_all1 <- data.frame(matrix(unlist(cis_eqtl_regions_all1), nrow=15, byrow=F),stringsAsFactors=FALSE)
}

cis_eqtl2 <- c("Ileum","Jejunum","Rumen","Adipose","Blood","Leukocyte","Lymph_node","Macrophage",
               "Milk_cell","Monocytes","Hypothalamus","Pituitary","Salivary_gland","Embryo",
               "Ovary","Oviduct","Uterus","Liver","Lung","Testis","Mammary","Muscle",
               "Intramuscular_fat","Skin_fibroblast")
cis_eqtl_regions_all2 <- cis_eqtl_regions_all1[,c(5,7,20,1,2,8,11,12,14,15,4,19,21,3,17,18,24,9,10,23,13,16,6,22)]

my_palette <- colorRampPalette(c("white","pink", "red"))(n = 299)
col_breaks = c(seq(0,2,length=100), seq(2.1,3,length=100), seq(3.1,6,length=100))

tiff(file="tissue_cis_eqtl.tiff", width = 24, height = 16.7, units = "cm", res = 600, pointsize = 8,compression= "lzw")
heatmap.2(as.matrix(cis_eqtl_regions_all2),cellnote = round(as.matrix(cis_eqtl_regions_all2),digits = 2), notecol="black", density.info="none", trace="none", margins =c(10,2),
          key.par=list(mar=c(3.5,1,5,0)), col=my_palette, breaks=col_breaks, dendrogram="none", labRow = NA, labCol = cis_eqtl2, cexRow=1.5,cexCol=1.5, key=F, keysize=2,key.title=NA,Colv=F,Rowv=F,
          symm=T,symkey=F,symbreaks=T, scale="none",srtCol=90, adjCol = c(0.95,0), key.xlab=expression(paste(italic(Overlap),sep=""))) 
dev.off()


sqtl1 <- gsub("_sqtl","",sqtl)
sqtl_regions_all1 <- list()
for (i in sqtl) {
  sqtl_regions1 <- Overlaps_ATAC_order1[, i]
  sqtl_regions2 <- Overlaps_ATAC_order2[, i]
  sqtl_regions3 <- Overlaps_ATAC_order3[, i]
  sqtl_regions4 <- Overlaps_ATAC_order4[, i]
  sqtl_regions_all <- cbind(sqtl_regions1, sqtl_regions2, sqtl_regions3, sqtl_regions4)
  sqtl_regions_all <- as.data.frame(sqtl_regions_all)
  sqtl_regions_all$mean <- apply(sqtl_regions_all, 1, mean)
  sqtl_regions_all1 <- cbind(sqtl_regions_all1, sqtl_regions_all$mean)
  sqtl_regions_all1 <- data.frame(matrix(unlist(sqtl_regions_all1), nrow=15, byrow=F),stringsAsFactors=FALSE)
}

sqtl2 <- c("Ileum","Jejunum","Rumen","Adipose","Blood","Leukocyte","Lymph_node","Macrophage",
           "Milk_cell","Monocytes","Hypothalamus","Pituitary","Salivary_gland","Embryo",
           "Ovary","Oviduct","Uterus","Liver","Lung","Testis","Mammary","Muscle_cross",
           "Muscle_indicus","Muscle","Muscle_taurus","Intramuscular_fat","Skin_fibroblast")
sqtl_regions_all2 <- sqtl_regions_all1[,c(5,7,23,1,2,8,11,12,14,15,4,22,24,3,20,21,27,9,10,26,13,16,17,18,19,6,25)]
my_palette <- colorRampPalette(c("white","pink", "red"))(n = 299)
col_breaks = c(seq(0,2,length=100), seq(2.1,3,length=100), seq(3.1,6,length=100))

tiff(file="Fig. 3e.tiff", width = 27, height = 16.7, units = "cm", res = 600, pointsize = 8,compression= "lzw")
heatmap.2(as.matrix(sqtl_regions_all2), cellnote = round(as.matrix(sqtl_regions_all2),digits = 2), notecol="black", density.info="none", trace="none", margins =c(10,2),
          key.par=list(mar=c(3.5,1,5,0)),col=my_palette,breaks=col_breaks,dendrogram="none",labRow = NA,labCol = sqtl2, cexRow=1.5,cexCol=1.5, key=F, keysize=2,key.title=NA,Colv=F,Rowv=F,
          symm=T,symkey=F,symbreaks=T, scale="none",srtCol=90, adjCol = c(0.95,0), key.xlab=expression(paste(italic(Overlap),sep=""))) 
dev.off()


# Fig. 3g
e_sig <- read.table("tissue_b_eall.txt", header = T, sep = "\t")
colnames(e_sig) <- c("Motif", "State", "Log10P")
e_sig$Log10P <- -e_sig$Log10P

tiff(file = "Fig. 3g.tiff", res = 300, width = 1400, height = 1850,compression = "lzw")
ggplot(e_sig, aes(x = State, y=Motif)) +
  geom_point(aes(size = Log10P),shape = 21, colour = "black",fill = "red", stroke = 1) +
  scale_size_continuous(range = c(1, 15),name = expression(paste(-log[10],italic(P),sep="")))+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 14),axis.text.y = element_text(size = 14))
dev.off()


# Fig. 3h
Overlaps_ATAC1 <- read.table("BW_15_overlap_rumen_ca_tpm.txt",header = T,stringsAsFactors = F,sep="\t")
Overlaps_ATAC_order1 <- Overlaps_ATAC1[c(8,10,9,7,5,13,6,4,12,14,11,15,3,2,1),]
colnames(Overlaps_ATAC_order1) <- gsub(".bed.gz","",colnames(Overlaps_ATAC_order1))
Overlaps_ATAC2 <- read.table("AW_15_overlap_rumen_ca_tpm.txt",header = T,stringsAsFactors = F,sep="\t")
Overlaps_ATAC_order2 <- Overlaps_ATAC2[c(8,10,9,7,5,13,6,4,12,14,11,15,3,2,1),]
colnames(Overlaps_ATAC_order2) <- gsub(".bed.gz","",colnames(Overlaps_ATAC_order2))
Overlaps_ATAC3 <- read.table("CO_overlap_rumen_ca_tpm.txt",header = T,stringsAsFactors = F,sep="\t")
Overlaps_ATAC_order3 <- Overlaps_ATAC3[c(8,10,9,7,5,13,6,4,12,14,11,15,3,2,1),]
colnames(Overlaps_ATAC_order3) <- gsub(".bed.gz","",colnames(Overlaps_ATAC_order3))
Overlaps_ATAC4 <- read.table("BT_15_overlap_rumen_ca_tpm.txt",header = T,stringsAsFactors = F,sep="\t")
Overlaps_ATAC_order4 <- Overlaps_ATAC4[c(8,10,9,7,5,13,6,4,12,14,11,15,3,2,1),]
colnames(Overlaps_ATAC_order4) <- gsub(".bed.gz","",colnames(Overlaps_ATAC_order4))

cate_regions_all1 <- list()
for (i in 3:7) {
  cate_regions1 <- Overlaps_ATAC_order1[, i]
  cate_regions2 <- Overlaps_ATAC_order2[, i]
  cate_regions3 <- Overlaps_ATAC_order3[, i]
  cate_regions4 <- Overlaps_ATAC_order4[, i]
  cate_regions_all <- cbind(cate_regions1, cate_regions2, cate_regions3, cate_regions4)
  cate_regions_all <- as.data.frame(cate_regions_all)
  cate_regions_all$mean <- apply(cate_regions_all, 1, mean)
  cate_regions_all$sd  <- apply(cate_regions_all[,1:4], 1, sd)
  cate_regions_all$state <- 1:15
  cate_regions_all$state <- as.factor(cate_regions_all$state)
  cate_regions_all$cate <- paste("cate", i-1)
  cate_regions_all1 <- rbind(cate_regions_all1, cate_regions_all)
}

Mycolors <- c("#80ca7f","#beafd4","#f6c085","#ea397e","#356fb2")

p <- ggplot(cate_regions_all1,aes(x=state, y=mean, color = cate, shape = cate, fill = cate)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) +
  scale_color_manual(values = rev(Mycolors)) +
  scale_shape_manual(values = c(19,19,19,19,19,19,19,19,19,19,19,19,19,19,19))+
  geom_line() + geom_point() + theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),legend.position="none",
        axis.line = element_line(colour = "black"),        
        axis.title.y = element_blank())

tiff("Fig. 3h1.tif", res = 300, width = 1000, height = 700, compression = "lzw")
p + coord_cartesian(ylim = c(0, 7.5)) + theme(legend.position='none')
dev.off()

tiff("Fig. 3h2.tif", res = 300, width = 1000, height = 300,compression = "lzw")
p + coord_cartesian(ylim = c(9.5, 50)) +  
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(),legend.position = 'none', 
        axis.line.x=element_line(colour="white"))
dev.off()


# Fig. 3i
library("ggpubr")
library("Rmisc")
library("gridExtra")
library(ggplot2)

rumen <- read.table("Rumen_gene1.txt", header = T, sep = " ", row.names = 1)
rumen <- t(rumen)
rumen <- as.data.frame(rumen)
rumen <- rumen[rowSums(rumen)!=0,]
rumen <- rumen[rowSums(rumen>=10)>5,colSums(rumen>=10)>1000]
rumen <- log10(rumen+1)
rumen_mad_gene <- apply(rumen, 1, mad, na.rm = TRUE)
rumen_mad_gene <- as.data.frame(rumen_mad_gene)
rumen_mad_gene$gene <- rownames(rumen_mad_gene)
colnames(rumen_mad_gene) <- c("mad", "gene")

rumen_mad_gene$category[rumen_mad_gene$mad < 0.1] <- 1
rumen_mad_gene$category[rumen_mad_gene$mad >=0.1 & rumen_mad_gene$mad <0.2] <- 2
rumen_mad_gene$category[rumen_mad_gene$mad >=0.2 & rumen_mad_gene$mad <0.3] <- 3
rumen_mad_gene$category[rumen_mad_gene$mad >=0.3 & rumen_mad_gene$mad <0.4] <- 4
rumen_mad_gene$category[rumen_mad_gene$mad >=0.4 & rumen_mad_gene$mad <0.5] <- 5
rumen_mad_gene$category[rumen_mad_gene$mad >=0.5 & rumen_mad_gene$mad <0.6] <- 6
rumen_mad_gene$category[rumen_mad_gene$mad >=0.6 & rumen_mad_gene$mad <0.7] <- 7
rumen_mad_gene$category[rumen_mad_gene$mad >=0.7 & rumen_mad_gene$mad <0.8] <- 8
rumen_mad_gene$category[rumen_mad_gene$mad >=0.9 & rumen_mad_gene$mad <1.0] <- 9
rumen_mad_gene$category[rumen_mad_gene$mad >=1.0] <- 10

for (i in 1:9) {
  cate <- rumen_mad_gene[rumen_mad_gene$category == i,]
  write.table(cate$gene, paste("mad",i,"_gene.txt", sep = ""), col.names = F, row.names = F, quote = F)
}

Overlaps_ATAC1 <- read.table("BW_15_overlap_rumen_ca_mad.txt",header = T,stringsAsFactors = F,sep="\t")
Overlaps_ATAC_order1 <- Overlaps_ATAC1[c(8,10,9,7,5,13,6,4,12,14,11,15,3,2,1),]
colnames(Overlaps_ATAC_order1) <- gsub(".bed.gz","",colnames(Overlaps_ATAC_order1))
Overlaps_ATAC2 <- read.table("AW_15_overlap_rumen_ca_mad.txt",header = T,stringsAsFactors = F,sep="\t")
Overlaps_ATAC_order2 <- Overlaps_ATAC2[c(8,10,9,7,5,13,6,4,12,14,11,15,3,2,1),]
colnames(Overlaps_ATAC_order2) <- gsub(".bed.gz","",colnames(Overlaps_ATAC_order2))
Overlaps_ATAC3 <- read.table("CO_overlap_rumen_ca_mad.txt",header = T,stringsAsFactors = F,sep="\t")
Overlaps_ATAC_order3 <- Overlaps_ATAC3[c(8,10,9,7,5,13,6,4,12,14,11,15,3,2,1),]
colnames(Overlaps_ATAC_order3) <- gsub(".bed.gz","",colnames(Overlaps_ATAC_order3))
Overlaps_ATAC4 <- read.table("BT_15_overlap_rumen_ca_mad.txt",header = T,stringsAsFactors = F,sep="\t")
Overlaps_ATAC_order4 <- Overlaps_ATAC4[c(8,10,9,7,5,13,6,4,12,14,11,15,3,2,1),]
colnames(Overlaps_ATAC_order4) <- gsub(".bed.gz","",colnames(Overlaps_ATAC_order4))

cate_regions_all1 <- list()
for (i in 3:11) {
  cate_regions1 <- Overlaps_ATAC_order1[, i]
  cate_regions2 <- Overlaps_ATAC_order2[, i]
  cate_regions3 <- Overlaps_ATAC_order3[, i]
  cate_regions4 <- Overlaps_ATAC_order4[, i]
  cate_regions_all <- cbind(cate_regions1, cate_regions2, cate_regions3, cate_regions4)
  cate_regions_all <- as.data.frame(cate_regions_all)
  cate_regions_all$mean <- apply(cate_regions_all, 1, mean)
  cate_regions_all$sd  <- apply(cate_regions_all[,1:4], 1, sd)
  cate_regions_all$state <- 1:15
  cate_regions_all$state <- as.factor(cate_regions_all$state)
  cate_regions_all$cate <- paste("cate", i-2)
  cate_regions_all1 <- rbind(cate_regions_all1, cate_regions_all)
}

Mycolors <- c("#fff2cc","#ffe699","#ffd966","#bf9000","#7f6000","#70ad47","#ff9999","#ff5050","#ff0000")
Mycolors <- c("black","red","green3","blue","cyan","magenta","yellow","gray","#ffe699","#ff5050")
Mycolors <- brewer.pal(10,"BrBG")
Mycolors <- colorRampPalette(c("blue","white","red"))(n = 10)
cate <- c("<0.1","[0.1,0.2)","[0.2,0.3)","[0.3,0.4)","[0.5,0.6)","[0.6,0.7)","[0.7,0.8)","[0.8,0.9)","[0.9,1.0)",">=1.0")

p <- ggplot(cate_regions_all1,aes(x=state, y=mean, color = cate, shape = cate, fill = cate)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) +
  scale_color_manual(values = Mycolors) +
  scale_shape_manual(values = c(19,19,19,19,19,19,19,19,19,19,19,19,19,19,19))+
  geom_line() + geom_point() + theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),legend.position="none",
        axis.line = element_line(colour = "black"),        
        axis.title.y = element_blank())

tiff("Fig. 3i.tif", res = 300, width = 1000, height = 700, compression = "lzw")
p + annotate("point", x=12,y=c(11,12,13,14,15,16,17,18,19,20),shape=15, color=Mycolors, size=1) +
  annotate("text", x=12.2, y=c(11,12,13,14,15,16,17,18,19,20), label=cate, hjust=0, size=1.9)
dev.off()



#========#
# Fig. 4 #
#========#
# Fig. 4a
library(data.table)
library(GenomicRanges)

Original_Order <- paste("E",rev(c(8,10,9,7,5,13,6,4,12,14,11,15,3,2,1)),sep = "")
names(Original_Order) <- c("15 Quies", "14 ReprPC","13 BivFlnk","12 ReprWkCTCF",
                           "11 ATAC","10 EnhWkCTCFATAC","9 EnhPoisATAC","8 EnhPois","7 EnhWk",
                           "6 EnhAATAC","5 EnhA", "4 TxFlnk","3 TssAFlnk", "2 TssAATACCTCF","1 TssA")

bw_segments <- read.table("BW_15_segments.bed", header = F, sep = "\t", stringsAsFactors = F)
bw_segments$length <- bw_segments$V3 - bw_segments$V2
bw_states_length <- tapply(bw_segments$length, bw_segments$V4, sum)
sum(bw_states_length)

aw_segments <- read.table("AW_15_segments.bed", header = F, sep = "\t", stringsAsFactors = F)
aw_segments$length <- aw_segments$V3 - aw_segments$V2
aw_states_length <- tapply(aw_segments$length, aw_segments$V4, sum)
sum(aw_states_length)

aw_induced <- ((aw_states_length - bw_states_length)/bw_states_length)
aw_induced_order <- aw_induced[c(1,8:15,2:7)]
aw_induced_order_1 <- aw_induced_order[c(8,10,9,7,5,13,6,4,12,14,11,15,3,2,1)]
names(aw_induced_order_1) <- rev(c("15 Quies", "14 ReprPC","13 BivFlnk","12 ReprWkCTCF",
                                   "11 ATAC","10 EnhWkCTCFATAC","9 EnhPoisATAC","8 EnhPois","7 EnhWk",
                                   "6 EnhAATAC","5 EnhA", "4 TxFlnk","3 TssAFlnk", "2 TssAATACCTCF","1 TssA"))

Mycolors <- c("white","#767171","#a9d18e","#8497b0","#c54611","#f4b183","#fff2cc","#ffe699","#ffd966","#bf9000","#7f6000","#70ad47","#ff9999","#ff5050","#ff0000")

tiff(file = "Fig. 4a.tif", res = 300, width = 1400, height = 1800,compression = "lzw")
par(mar=c(6,9,2,2))
barplot(rev(aw_induced_order_1),las=2, horiz = T, xlim = c(-1,6), cex.names = 0.95, xlab=paste("Relative proportion of altered regions\n between BW and AW", sep = ""), col = Mycolors)
abline(v = 0)
dev.off()


# Fig. 4b
TSS_up <- read.table("tissue_tss_up.txt",header = T,stringsAsFactors = F,sep="\t")
TSS_up <- TSS_up[,-1]
TSS_up <- TSS_up[c(8,10,9,7,5,13,6,4,12,14,11,15,3,2,1),]
State_names <- c("15 Quies", "14 ReprPC","13 BivFlnk","12 ReprWkCTCF",
                 "11 ATAC","10 EnhWkCTCFATAC","9 EnhPoisATAC","8 EnhPois","7 EnhWk",
                 "6 EnhAATAC","5 EnhA", "4 TxFlnk","3 TssAFlnk", "2 TssAATACCTCF","1 TssA")
row.names(TSS_up) <- rev(State_names)
my_palette <- colorRampPalette(c("blue","white","red"))(n = 299)
col_breaks = c(seq(-20,-5,length=100), seq(-4.9,5,length=100), seq(5.1,20,length=100))
tiff(file = "Fig. 4b1.tiff", res = 300, width = 1800, height = 1200,compression = "lzw")
heatmap.2(as.matrix(TSS_up), notecol="black", density.info="none", trace="none", margins =c(4,12), key.par=list(mar=c(2,0.5,0,0)), col=my_palette, breaks=col_breaks,
          dendrogram="none", labRow = rev(State_names), labCol = seq(-2000,2000,by = 200), cexRow=1.2,cexCol=1.2, key=T, keysize=0.6,key.title=NA,Colv=F,Rowv=F,
          symm=T,symkey=F,symbreaks=T, scale="none",srtCol=90, adjCol = c(0.95,0), key.xlab=expression(paste(italic(EnriFold),sep="")))
dev.off()


TSS_down<- read.table("tissue_tss_down.txt",header = T,stringsAsFactors = F,sep="\t")
TSS_down<- TSS_down[,-1]
TSS_down<- TSS_down[c(8,10,9,7,5,13,6,4,12,14,11,15,3,2,1),]
State_names <- c("15 Quies", "14 ReprPC","13 BivFlnk","12 ReprWkCTCF",
                 "11 ATAC","10 EnhWkCTCFATAC","9 EnhPoisATAC","8 EnhPois","7 EnhWk",
                 "6 EnhAATAC","5 EnhA", "4 TxFlnk","3 TssAFlnk", "2 TssAATACCTCF","1 TssA")
row.names(TSS_down) <- rev(State_names)
my_palette <- colorRampPalette(c("blue","white","red"))(n = 299)
col_breaks = c(seq(-20,-5,length=100), seq(-4.9,5,length=100), seq(5.1,20,length=100))
tiff(file = "Fig. 4b2.tiff", res = 300, width = 1800, height = 1200,compression = "lzw")
heatmap.2(as.matrix(TSS_down), notecol="black", density.info="none", trace="none", margins =c(4,12), key.par=list(mar=c(2,0.5,0,0)), col=my_palette, breaks=col_breaks,
          dendrogram="none", labRow = rev(State_names), labCol = seq(-2000,2000,by = 200), cexRow=1.2,cexCol=1.2, key=T, keysize=0.6,key.title=NA,Colv=F,Rowv=F,
          symm=T,symkey=F,symbreaks=T, scale="none",srtCol=90, adjCol = c(0.95,0), key.xlab=expression(paste(italic(EnriFold),sep="")))
dev.off()


# Fig. 4c
library(ggplot2)
module_gene <- read.table("pathway.txt", sep = "\t", header = T, row.names = 1)
color <- c("blue", "red")
tiff(file = "Fig. 4c.tiff", width = 11, height = 7, units = "cm", res = 600, pointsize = 8, compression= "lzw")
p <- ggplot(module_gene, aes(x = Module, y = Description))+
  geom_point(aes(color = p.adjust, size = Count))+
  scale_color_gradient(low = "red", high = "blue")+
  xlab("") + ylab("Pathways") + theme_bw()
p + theme(axis.text.x = element_text(size=14, color = color, angle = 0), axis.text.y = element_text(size=10))
dev.off()



#========#
# Fig. 5 #
#========#
# Fig. 5a
co_segments <- read.table("CO_15_segments.bed", header = F, sep = "\t", stringsAsFactors = F)
co_segments$length <- co_segments$V3 - co_segments$V2
co_states_length <- tapply(co_segments$length, co_segments$V4, sum)

bt_segments <- read.table("BT_15_segments.bed", header = F, sep = "\t", stringsAsFactors = F)
bt_segments$length <- bt_segments$V3 - bt_segments$V2
bt_states_length <- tapply(bt_segments$length, bt_segments$V4, sum)
bt_induced <- ((bt_states_length - co_states_length)/co_states_length)
bt_induced_order <- bt_induced[c(1,8:15,2:7)]
bt_induced_order_1 <- bt_induced_order[c(8,10,9,7,5,13,6,4,12,14,11,15,3,2,1)]
names(bt_induced_order_1) <- rev(c("15 Quies", "14 ReprPC","13 BivFlnk","12 ReprWkCTCF",
                                   "11 ATAC","10 EnhWkCTCFATAC","9 EnhPoisATAC","8 EnhPois","7 EnhWk",
                                   "6 EnhAATAC","5 EnhA", "4 TxFlnk","3 TssAFlnk", "2 TssAATACCTCF","1 TssA"))

Mycolors <- c("white","#767171","#a9d18e","#8497b0","#c54611","#f4b183","#fff2cc","#ffe699","#ffd966","#bf9000","#7f6000","#70ad47","#ff9999","#ff5050","#ff0000")
tiff(file = "Fig. 5a.tif", res = 300, width = 1400, height = 1800,compression = "lzw")
par(mar=c(6,9,2,2))
barplot(rev(bt_induced_order_1),las=2, horiz = T, xlim = c(-1,6),cex.names = 0.95,xlab=paste("Relative proportion of altered regions\n between co and bt", sep = ""), col = Mycolors)
abline(v = 0)
dev.off()
change <- data.frame(aw_induced_order_1, bt_induced_order_1)
change$trend <- change$aw_induced_order_1 * change$bt_induced_order_1



#========#
# Fig. 6 #
#========#
# Fig. 6c
emissions <- read.table("emissions_15.txt",header = T,stringsAsFactors = F,sep="\t")
emissions <- emissions[,c("H3K4me3","H3K4me1","H3K27ac","ATAC","CTCF","H3K27me3")]
emissions_order <- emissions[c(4,5,11,2,7,3,1,8,6,9,15,10,12:14),]

library(gplots)
library(RColorBrewer)
library(data.table)

my_palette <- colorRampPalette(c("white","darkblue"))(n = 199)
col_breaks = c(seq(0,0.5,length=100), seq(0.51,1,length=100))

tiff(file="Fig. 6c.tiff", width = 6, height = 16.7, units = "cm", res = 600, pointsize = 8,compression= "lzw")
heatmap.2(as.matrix(emissions_order), notecol="black", density.info="none", trace="none", margins =c(7,2), key.par=list(mar=c(3.5,1,5,0)), col=my_palette, breaks=col_breaks,
          dendrogram="none", labRow = NA, labCol = colnames(emissions), cexRow=1.5,cexCol=1.5, key=T, keysize=1.2,key.title=NA,Colv=F,Rowv=F,
          symm=T,symkey=F,symbreaks=T, scale="none",srtCol=90, adjCol = c(0.95,0), key.xlab=expression(paste(italic(Emissions),sep=)))
dev.off()


# Fig. 6d
library(dplyr)
library(forcats)
library(ggplot2)

gc <- read.table("genome_coverage.txt", header = T, sep = "\t")
gc <- as.data.frame(gc)
gc <- gc[,-1]
gc$mean <- apply(gc, 1, mean)
gc$sd <- apply(gc[,1:8], 1, sd)
gc$state <- 1:15
gc$state <- as.factor(gc$state)

Mycolors <- c("white","#9f9b9b","#767171","#a9d18e","#8497b0","#f4b183","#c55a11","#ffe699","#ffd966","#cca533","#bf9000","#7f6000","#ff9999","#ff5050","#ff0000")

gc1 <- gc[,c(9,9)]
gc1 <- as.matrix(gc1)
my_palette <- colorRampPalette(c("white","darkblue"))(n = 199)
col_breaks = c(seq(0,2,length=100), seq(2.1,4,length=100))

tiff("Fig. 6d.tiff", width = 9, height = 16.7, units = "cm", res = 600, pointsize = 8,compression= "lzw")
heatmap.2(gc1,
          cellnote = gc[,c(12,12)],
          notecol=c("white","white","white","white","black","black","black","black",
                    "black","black","black","black","black","black","black",
                    "black","black","black","black","black","black","black","black",
                    "black","black","black","black","black","black","black"),
          notecex=1.5, density.info="none", trace="none", margins =c(7,2), key.par=list(mar=c(3.5,1,5,0)), col=my_palette, breaks=col_breaks, dendrogram="none", 
          labRow = NA, labCol = NA, cexRow=1.5,cexCol=1.5, key=F, keysize=1.2,key.title=NA,Colv=F,Rowv=F, symm=T,symkey=F,symbreaks=T, scale="none",srtCol=90, adjCol = c(0.95,0), 
          key.xlab=expression(paste(italic(Emissions),sep=)))  
dev.off()


# Fig. 6e
Overlaps_ATAC1 <- read.table("Overlap_ATAC_Annotation1.txt",header = T,row.names = 1,sep="\t")
Overlaps_ATAC2 <- read.table("Overlap_ATAC_Annotation2.txt",header = T,row.names = 1,sep="\t")
Overlaps_ATAC3 <- read.table("Overlap_ATAC_Annotation3.txt",header = T,row.names = 1,sep="\t")
Overlaps_ATAC4 <- read.table("Overlap_ATAC_Annotation4.txt",header = T,row.names = 1,sep="\t")
Overlaps_ATAC5 <- read.table("Overlap_ATAC_Annotation5.txt",header = T,row.names = 1,sep="\t")
Overlaps_ATAC6 <- read.table("Overlap_ATAC_Annotation6.txt",header = T,row.names = 1,sep="\t")
Overlaps_ATAC <- rbind(Overlaps_ATAC1,Overlaps_ATAC2,Overlaps_ATAC3,Overlaps_ATAC4,Overlaps_ATAC5,Overlaps_ATAC6)
rownames(Overlaps_ATAC) <- gsub("1","Colon",rownames(Overlaps_ATAC))
rownames(Overlaps_ATAC) <- gsub("2","Foreskin",rownames(Overlaps_ATAC))
rownames(Overlaps_ATAC) <- gsub("3","Heart",rownames(Overlaps_ATAC))
rownames(Overlaps_ATAC) <- gsub("4","Pancreas",rownames(Overlaps_ATAC))
rownames(Overlaps_ATAC) <- gsub("5","Stomach",rownames(Overlaps_ATAC))
rownames(Overlaps_ATAC) <- gsub("6","Testis",rownames(Overlaps_ATAC))
Overlaps_ATAC <- Overlaps_ATAC[c(1,6,11,16,21,26,2,7,12,17,22,27,3,8,13,18,23,28,4,9,14,19,24,29,5,10,15,20,25,30),]

Overlaps_ATAC_Pr_cattle <- Overlaps_ATAC[,8]+Overlaps_ATAC[,9]+Overlaps_ATAC[,10]+Overlaps_ATAC[,11]
Overlaps_ATAC_En_cattle <- Overlaps_ATAC[,12]+Overlaps_ATAC[,13]+Overlaps_ATAC[,14]+Overlaps_ATAC[,15]
Overlaps_ATAC_ATAC_cattle <- Overlaps_ATAC[,18]+Overlaps_ATAC[,16]
Overlaps_ATAC_CTCF_cattle <- Overlaps_ATAC[,17]+Overlaps_ATAC[,19]
Overlaps_ATAC_Bi_cattle <- Overlaps_ATAC[,20]
Overlaps_ATAC_state_cattle <- cbind(Overlaps_ATAC[,1:7], Overlaps_ATAC_state_cattle)

Mynames <- c("CpG_island","Exon","Gene", "TES","TSS","TSS2kb", "Pr", "En", "ATAC", "CTCF", "Bi")
my_palette <- colorRampPalette(c("white","red"))(n = 399)
col_breaks = seq(0,241,length=400)

tiff(file = "Fig. 6e.tiff", res = 300, width = 2700, height = 4300,compression = "lzw")
heatmap.2(as.matrix(Overlaps_ATAC_state_cattle[,-c(1)]),
          cellnote = round(as.matrix(Overlaps_ATAC_state_cattle[,-c(1)]),digits = 2), notecol="black", density.info="none", trace="none", margins =c(14,8), key.par=list(mar=c(3.5,1,5,3.5)),
          col=my_palette, breaks=col_breaks, dendrogram="none", labRow = rownames(Overlaps_ATAC_state_cattle), labCol = Mynames, cexRow=1,cexCol=1.5, key=F, 
          keysize=1.4,key.title=NA,Colv=F,Rowv=F, symm=T,symkey=F,symbreaks=T, scale="none",srtCol=90, adjCol = c(0.95,0), key.xlab=expression(paste(italic(Overlap),sep=""))) 
dev.off()
