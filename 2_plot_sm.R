#=========#
# Fig. S1 #
#=========#

# Fig. S1b
All_peaks <- read.table("All_peaks.bed", header = T, sep = "\t")

unique(All_peaks$Sample)
tiff(file = "Fig. S1b1.tif",res = 300, width = 1800, height = 2400,compression = "lzw")
par(mar = c(4,8,1,2))
XX <- barplot(rev(table(All_peaks$Sample)),las=1, horiz = T, xlab="Number of peaks",col = "deepskyblue1",xlim = c(0,150000))
text(y = XX, x = rev(table(All_peaks$Sample)), label = rev(table(All_peaks$Sample)), pos = 4, cex = 0.8, col = "red")
dev.off()

Genome_Length <- 2649685036
All_peaks_coverage <- tapply(All_peaks$Length,All_peaks$Sample,sum)
All_peaks_coverage_perc <- (All_peaks_coverage/Genome_Length)*100
tiff(file = "Fig. S1b2.tif",res = 300, width = 1800, height = 2400,compression = "lzw")
par(mar = c(4,8,1,2))
XX <- barplot(rev(All_peaks_coverage_perc),las=1, horiz = T,xlim=c(0,15),  xlab="Genome coverage(%)",col = "black")
text(y = XX, x = rev(All_peaks_coverage_perc), label = round(rev(All_peaks_coverage_perc),digits = 2), pos = 4, cex = 0.8, col = "red")
dev.off()


# Fig. S1f. Refer to Fig. 1b


#=========#
# Fig. S2 #
#=========#
Transiton <- read.table("transitions_15.txt",header = T,stringsAsFactors = F,sep="\t")
Transiton <- Transiton[,-1]
Transiton <- Transiton[c(8,10,9,7,5,13,6,4,12,14,11,15,3,2,1),c(8,10,9,7,5,13,6,4,12,14,11,15,3,2,1)]
State_names <- c("15 Quies", "14 ReprPC","13 BivFlnk","12 ReprWkCTCF",
                 "11 ATAC","10 EnhWkCTCFATAC","9 EnhPoisATAC","8 EnhPois","7 EnhWk",
                 "6 EnhAATAC","5 EnhA", "4 TxFlnk","3 TssAFlnk", "2 TssAATACCTCF","1 TssA")
colnames(Transiton) <- rev(State_names)
row.names(Transiton) <- rev(State_names)
my_palette <- colorRampPalette(c("white","red"))(n = 199)
col_breaks = c(seq(0,0.1,length=100), seq(0.11,1,length=100))
tiff(file = "Fig. S2.tiff", res = 300, width = 2450, height = 2450,compression = "lzw")
heatmap.2(as.matrix(Transiton),
          cellnote = round(Transiton,digits = 2), notecol="black", density.info="none", trace="none", margins =c(14,14), key.par=list(mar=c(2,1.5,0,0)), col=my_palette,      
          breaks=col_breaks,  dendrogram="none", cexRow=1.5,cexCol=1.5, key=T, keysize=0.6,key.title=NA,Colv=F,Rowv=F, symm=T,symkey=F,symbreaks=T, 
          scale="none",srtCol=90, adjCol = c(0.95,0), key.xlab=expression(paste(italic(Transiton),sep="")))     


#=========#
# Fig. S3 #
#=========#
# Fig. S3a. Refer to Fig. 3d

# Fig. S3b
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

rumen_mean_gene <- apply(rumen, 1, mean, na.rm = TRUE)
rumen_mean_gene <- as.data.frame(rumen_mean_gene)
rumen_mean_gene$gene <- rownames(rumen_mean_gene)
colnames(rumen_mean_gene) <- c("mean", "gene")

rumen_mean_mad <- merge(rumen_mean_gene, rumen_mad_gene)
rumen_mean_mad <- rumen_mean_mad[order(rumen_mean_mad$mad, decreasing = T),]
rumen_mean_mad <- rumen_mean_mad[-1,]

p1 <- ggplot(rumen_mean_mad, aes(mad, mean)) +
  geom_point() +
  geom_smooth(method = lm, level = 0.95, fullrange = TRUE) +
  stat_cor(method = "pearson", label.x = 0.75, label.y = 5, size=4.5) +
  labs(x="MAD",y="Mean") +
  theme_bw()+ 
  theme(axis.text=element_text(size=12), axis.title = element_blank(), axis.title.x = element_text(), axis.title.y = element_text())
p2 <- ggplot(rumen_mean_mad, aes(x=mad)) + geom_density(fill='#bdd7ee') + theme_void()
p3 <- ggplot(rumen_mean_mad, aes(y=mean)) + geom_density(fill='#bdd7ee') + theme_void()

tiff(file = "Fig. S3b.tiff",res = 300, width = 1200, height = 1200, compression = "lzw")
p1 %>% insert_top(p2, height=.2) %>% insert_right(p3, width=.2)
dev.off()


# Fig. S3d
rumen_mean_mad1 <- rumen_mean_mad[rumen_mean_mad$mad > 1.1,]
rumen_mean_mad1 <- rumen_mean_mad1[rumen_mean_mad1$mean < 2,]
rumen_mean_mad2 <- rumen_mean_mad[rumen_mean_mad$mean > 3,]
rumen_mean_mad2 <- rumen_mean_mad2[rumen_mean_mad2$mad < 1,]
write.table(rumen_mean_mad1$gene, "top_mad_gene.txt", quote = F, col.names = F, row.names = F)
write.table(rumen_mean_mad2$gene, "top_mean_gene.txt", quote = F, col.names = F, row.names = F)

rumen_mean_mad <- rumen_mean_mad[order(rumen_mean_mad$mad, decreasing = F),]
rumen_mean_mad3 <- head(rumen_mean_mad, n = 130)
rumen_mean_mad3 <- rumen_mean_mad3[rumen_mean_mad3$mad < 0.5,]
rumen_mean_mad4 <- rumen_mean_mad[rumen_mean_mad$mad > 0.13,]
rumen_mean_mad4 <- rumen_mean_mad4[rumen_mean_mad4$mad < 0.4,]
rumen_mean_mad4 <- rumen_mean_mad4[rumen_mean_mad4$mean > 1,]
rumen_mean_mad4 <- rumen_mean_mad4[rumen_mean_mad4$mean > 1,]
rumen_mean_mad4 <- rumen_mean_mad4[rumen_mean_mad4$mean < 1.01,]
write.table(rumen_mean_mad3$gene, "bot_mad_gene.txt", quote = F, col.names = F, row.names = F)
write.table(rumen_mean_mad4$gene, "bot_mean_gene.txt", quote = F, col.names = F, row.names = F)

library(clusterProfiler)
library(org.Bt.eg.db)

mad_gene1 <- rumen_mean_mad1$gene
mad_gene1 <- as.character(mad_gene1)
mad_gene1_id <- bitr(mad_gene1, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Bt.eg.db")
mad_gene1_list <- mad_gene1_id$ENTREZID
mad_gene1_list[duplicated(mad_gene1_list)]
mad_gene1_go <- enrichGO(mad_gene1_list, OrgDb = org.Bt.eg.db, ont='ALL', pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.2, keyType = 'ENTREZID')
tiff(file = "Fig. S3d1.tiff", width = 35, height = 30, units = "cm", res = 600, pointsize = 8,compression= "lzw")
dotplot(mad_gene1_go, showCategory=50, font.size =20) + theme(legend.text = element_text(size = 15), legend.title = element_text(size = 18))
dev.off()

mad_gene1_kegg <- enrichKEGG(mad_gene1_list, organism = 'bta', keyType = 'kegg', pvalueCutoff = 0.05, pAdjustMethod = 'BH', minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
tiff(file = "Fig. S3d2.tiff", width = 25, height = 20, units = "cm", res = 600, pointsize = 8, compression= "lzw")
dotplot(mad_gene1_kegg, showCategory = 30, font.size = 20) + theme(legend.text = element_text(size = 15),legend.title = element_text(size = 18))
dev.off()

mean_gene1 <- rumen_mean_mad2$gene
mean_gene1 <- as.character(mean_gene1)
mean_gene1_id <- bitr(mean_gene1, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Bt.eg.db")
mean_gene1_list <- mean_gene1_id$ENTREZID
mean_gene1_list[duplicated(mean_gene1_list)]
mean_gene1_go <- enrichGO(mean_gene1_list, OrgDb = org.Bt.eg.db, ont='ALL', pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.2, keyType = 'ENTREZID')
tiff(file = "Fig. S3d3.tiff", width = 35, height = 30, units = "cm", res = 600, pointsize = 8,compression= "lzw")
dotplot(mean_gene1_go, showCategory=50, font.size =20) + theme(legend.text = element_text(size = 15), legend.title = element_text(size = 18))
dev.off()

mean_gene1_kegg <- enrichKEGG(mean_gene1_list, organism = 'bta', keyType = 'kegg', pvalueCutoff = 0.05, pAdjustMethod = 'BH', minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
tiff(file = "Fig. S3d4.tiff", width = 25, height = 20, units = "cm", res = 600, pointsize = 8, compression= "lzw")
dotplot(mean_gene1_kegg, showCategory = 30, font.size = 20) + theme(legend.text = element_text(size = 15),legend.title = element_text(size = 18))
dev.off()


# Fig. S3c
library(dplyr)
library(forcats)
library(ggplot2)

Overlaps_ATAC1 <- read.table("BW_15_overlap_rumen.txt",header = T,stringsAsFactors = F,sep="\t")
Overlaps_ATAC_order1 <- Overlaps_ATAC1[c(8,10,9,7,5,13,6,4,12,14,11,15,3,2,1),]
colnames(Overlaps_ATAC_order1) <- gsub(".bed.gz","",colnames(Overlaps_ATAC_order1))
Overlaps_ATAC2 <- read.table("AW_15_overlap_rumen.txt",header = T,stringsAsFactors = F,sep="\t")
Overlaps_ATAC_order2 <- Overlaps_ATAC2[c(8,10,9,7,5,13,6,4,12,14,11,15,3,2,1),]
colnames(Overlaps_ATAC_order2) <- gsub(".bed.gz","",colnames(Overlaps_ATAC_order2))
Overlaps_ATAC3 <- read.table("CO_overlap_rumen.txt",header = T,stringsAsFactors = F,sep="\t")
Overlaps_ATAC_order3 <- Overlaps_ATAC3[c(8,10,9,7,5,13,6,4,12,14,11,15,3,2,1),]
colnames(Overlaps_ATAC_order3) <- gsub(".bed.gz","",colnames(Overlaps_ATAC_order3))
Overlaps_ATAC4 <- read.table("BT_15_overlap_rumen.txt",header = T,stringsAsFactors = F,sep="\t")
Overlaps_ATAC_order4 <- Overlaps_ATAC4[c(8,10,9,7,5,13,6,4,12,14,11,15,3,2,1),]
colnames(Overlaps_ATAC_order4) <- gsub(".bed.gz","",colnames(Overlaps_ATAC_order4))

Mycolors <- c("white","#767171","#a9d18e","#8497b0","#c54611","#f4b183","#fff2cc","#ffe699","#ffd966","#bf9000","#7f6000","#70ad47","#ff9999","#ff5050","#ff0000")

top_mad1 <- Overlaps_ATAC_order1[,3]
top_mad2 <- Overlaps_ATAC_order2[,3]
top_mad3 <- Overlaps_ATAC_order3[,3]
top_mad4 <- Overlaps_ATAC_order4[,3]
top_mad_all <- cbind(top_mad1, top_mad2, top_mad3, top_mad4)
top_mad_all <- as.data.frame(top_mad_all)
top_mad_all$mean <- apply(top_mad_all, 1, mean)
top_mad_all$sd <- apply(top_mad_all[,1:4], 1, sd)
top_mad_all$state <- 1:15
top_mad_all$state <- as.factor(top_mad_all$state)
tiff(file = "Fig. S3c1.tif",res = 300, width = 1000, height = 1850,compression = "lzw")
par(mar=c(4,4,2,2))
print(top_mad_all %>%
        mutate(state = fct_reorder(state, desc(state))) %>%
        ggplot(aes(x=state, y=mean)) +
        geom_bar(stat="identity", fill = rev(Mycolors), alpha=.8, width=1, color = "black") +
        geom_errorbar(aes(x=state, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="grey", alpha=0.9, size=1) + 
        coord_flip() + ylim(0, 80) +
        theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        xlab("") + ylab(paste("Enrichment for top MAD genes\n of the rumen in cattle",sep="")))
dev.off()

top_mean1 <- Overlaps_ATAC_order1[,4]
top_mean2 <- Overlaps_ATAC_order2[,4]
top_mean3 <- Overlaps_ATAC_order3[,4]
top_mean4 <- Overlaps_ATAC_order4[,4]
top_mean_all <- cbind(top_mean1, top_mean2, top_mean3, top_mean4)
top_mean_all <- as.data.frame(top_mean_all)
top_mean_all$mean <- apply(top_mean_all, 1, mean)
top_mean_all$sd <- apply(top_mean_all[,1:4], 1, sd)
top_mean_all$state <- 1:15
top_mean_all$state <- as.factor(top_mean_all$state)
tiff(file = "Fig. S3c2.tif",res = 300, width = 1000, height = 1850,compression = "lzw")
par(mar=c(4,4,2,2))
print(top_mean_all %>%
        mutate(state = fct_reorder(state, desc(state))) %>%
        ggplot(aes(x=state, y=mean)) +
        geom_bar(stat="identity", fill = rev(Mycolors), alpha=.8, width=1, color = "black") +
        geom_errorbar(aes(x=state, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="grey", alpha=0.9, size=1) + 
        coord_flip() + ylim(0, 80) +
        theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        xlab("") + ylab(paste("Enrichment for top TPM genes\n of the rumen in cattle",sep="")))
dev.off()

bot_mad1 <- Overlaps_ATAC_order1[,5]
bot_mad2 <- Overlaps_ATAC_order2[,5]
bot_mad3 <- Overlaps_ATAC_order3[,5]
bot_mad4 <- Overlaps_ATAC_order4[,5]
bot_mad_all <- cbind(bot_mad1, bot_mad2, bot_mad3, bot_mad4)
bot_mad_all <- as.data.frame(bot_mad_all)
bot_mad_all$mean <- apply(bot_mad_all, 1, mean)
bot_mad_all$sd <- apply(bot_mad_all[,1:4], 1, sd)
bot_mad_all$state <- 1:15
bot_mad_all$state <- as.factor(bot_mad_all$state)
tiff(file = "Fig. S3c3.tif",res = 300, width = 1000, height = 1850,compression = "lzw")
par(mar=c(4,4,2,2))
print(bot_mad_all %>%
        mutate(state = fct_reorder(state, desc(state))) %>%
        ggplot(aes(x=state, y=mean)) +
        geom_bar(stat="identity", fill = rev(Mycolors), alpha=.8, width=1, color = "black") +
        geom_errorbar(aes(x=state, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="grey", alpha=0.9, size=1) + 
        coord_flip() + ylim(0, 80) +
        theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        xlab("") + ylab(paste("Enrichment for bot MAD genes\n of the rumen in cattle",sep="")))
dev.off()

bot_mean1 <- Overlaps_ATAC_order1[,6]
bot_mean2 <- Overlaps_ATAC_order2[,6]
bot_mean3 <- Overlaps_ATAC_order3[,6]
bot_mean4 <- Overlaps_ATAC_order4[,6]
bot_mean_all <- cbind(bot_mean1, bot_mean2, bot_mean3, bot_mean4)
bot_mean_all <- as.data.frame(bot_mean_all)
bot_mean_all$mean <- apply(bot_mean_all, 1, mean)
bot_mean_all$sd <- apply(bot_mean_all[,1:4], 1, sd)
bot_mean_all$state <- 1:15
bot_mean_all$state <- as.factor(bot_mean_all$state)
tiff(file = "Fig. S3c4.tif",res = 300, width = 1000, height = 1850,compression = "lzw")
par(mar=c(4,4,2,2))
print(bot_mean_all %>%
        mutate(state = fct_reorder(state, desc(state))) %>%
        ggplot(aes(x=state, y=mean)) +
        geom_bar(stat="identity", fill = rev(Mycolors), alpha=.8, width=1, color = "black") +
        geom_errorbar(aes(x=state, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="grey", alpha=0.9, size=1) + 
        coord_flip() + ylim(0, 80) +
        theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        xlab("") + ylab(paste("Enrichment for bot TPM genes\n of the rumen in cattle",sep="")))
dev.off()



#=========#
# Fig. S5 #
#=========#
tissue <- list.files(getwd(), pattern = '*overlap.txt$')
for (i in 1:6) {
  Overlaps_ATAC <- read.table(tissue[i],header = T,stringsAsFactors = F,sep="\t")
  Overlaps_ATAC <- Overlaps_ATAC[,c(1:8,15:23,9:14)]
  Overlaps_ATAC <- Overlaps_ATAC[c(4,5,11,2,7,3,1,8,6,9,15,10,12:14),]
  Overlaps_ATAC <- Overlaps_ATAC[,c(1:8,16,18,17,15,13,21,14,12,20,22,19,23,11,10,9)]
  State_names <- c("15 Quies", "14 ReprPC","13 BivFlnk","12 ReprWkCTCF",
                   "11 ATAC","10 EnhWkCTCFATAC","9 EnhPoisATAC","8 EnhPois","7 EnhWk",
                   "6 EnhAATAC","5 EnhA", "4 TxFlnk","3 TssAFlnk", "2 TssAATACCTCF","1 TssA")
  Mynames <- c("CpG_island","Exon","Gene", "TES","TSS","TSS2kb", rev(State_names))
  my_palette <- colorRampPalette(c("white","red"))(n = 199)
  col_breaks = seq(0,130,length=200)
  
  tiff(file=paste("Fig. S5",i,".tiff",sep = ""), width = 20, height = 16.7, units = "cm", res = 600, pointsize = 8,compression= "lzw")
  heatmap.2(as.matrix(Overlaps_ATAC[,-c(1,2)]),
            cellnote = round(as.matrix(Overlaps_ATAC[,-c(1,2)]),digits = 2), notecol="black", density.info="none", trace="none", margins =c(14,2), key.par=list(mar=c(3.5,1,5,3.5)),
            col=my_palette, breaks=col_breaks, dendrogram="none", labRow = NA, labCol = Mynames, cexRow=1.5,cexCol=1.5, key=F, keysize=1.4,key.title=NA,Colv=F,Rowv=F,
            symm=T,symkey=F,symbreaks=T, scale="none",srtCol=90, adjCol = c(0.95,0), key.xlab=expression(paste(italic(Overlap),sep=""))) 
  dev.off()
  
  Overlaps_ATAC_new <- Overlaps_ATAC[,-1]
  Overlaps_ATAC_promoter <- Overlaps_ATAC_new[1,]+Overlaps_ATAC_new[2,]+Overlaps_ATAC_new[3,]
  Overlaps_ATAC_enhancer <- Overlaps_ATAC_new[4,]+Overlaps_ATAC_new[5,]+Overlaps_ATAC_new[6,]+Overlaps_ATAC_new[7,]+Overlaps_ATAC_new[8,]+Overlaps_ATAC_new[9,]
  Overlaps_ATAC_ATAC <- Overlaps_ATAC_new[10,]
  Overlaps_ATAC_CTCF <- Overlaps_ATAC_new[11,]
  Overlaps_ATAC_BivFlnk <- Overlaps_ATAC_new[12,]
  Overlaps_ATAC_new <- rbind(Overlaps_ATAC_promoter,Overlaps_ATAC_enhancer,Overlaps_ATAC_ATAC,Overlaps_ATAC_CTCF,Overlaps_ATAC_BivFlnk)
  rownames(Overlaps_ATAC_new) <- c(paste("Pr_",i,sep=""),paste("En_",i,sep=""),paste("ATAC_",i,sep=""),paste("CTCF_",i,sep=""),paste("Bi_",i,sep=""))
  write.table(Overlaps_ATAC_new, paste("Overlap_ATAC_Annotation",i,".txt",sep=""), quote = F, sep = "\t")
}
