

#### loading dataset #### 
# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"; path <- paste(urld, "acc=GSE198256", "file=GSE198256_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
GSE198256_count <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)

#### retreive metadata #### 
library(GEOquery)
gds <- getGEO("GSE198256")
Meta_GSE198256 <- pData(gds$GSE198256_series_matrix.txt.gz@phenoData)
Meta_GSE198256 <- Meta_GSE198256[,c("title","source_name_ch1","characteristics_ch1","characteristics_ch1.1","description","cell type:ch1","disease state:ch1")]
Factors_GSE198256 <- Meta_GSE198256[,c("disease state:ch1")]
Factors_GSE198256 <- data.frame(Meta_GSE198256 [colnames(GSE198256_count),c("disease state:ch1")])
colnames(Factors_GSE198256)[1]<- "Group"
# (nrow(Factors_GSE198256) == ncol(GSE198256_count))
Factors_GSE198256[Factors_GSE198256=="Covid19: Acute infection"] <- "Covid19AI"
Factors_GSE198256[Factors_GSE198256=="Covid19: Recovery 3Mo"] <- "Covid193Mo"
Factors_GSE198256[Factors_GSE198256=="Covid19: Recovery 6Mo"] <- "Covid196Mo"
rownames(Factors_GSE198256) <- rownames(Meta_GSE198256)
Factors_GSE198256$Group <- factor(Factors_GSE198256$Group)


#### Metadata of genes #### 
# 1) Biological information: GC, Gene length...etc. (using Biomart: https://asia.ensembl.org/info/data/biomart/index.html)
genes_annotations <- read.table(file = "~/Library/Mobile Documents/com~apple~CloudDocs/BESE394/Week_1/mart_export.txt", sep = "\t", header = T)
# Filter the information
genes_annotations <- genes_annotations[genes_annotations$Chromosome %in% c(as.character(1:22) ,"X","Y"),]
# sum(rownames(GSE198256_count) %in% genes_annotations_filt$Entrezgene)
genes_annotations_filt <- genes_annotations[!duplicated(genes_annotations$Entrezgene),] # removing duplicates 
# sum(rownames(GSE198256_count) %in% genes_annotations_filt$Entrezgene)
# sum(genes_annotations_filt$Entrezgene %in% rownames(GSE198256_count))
# Overlap between annotation and genes
rownames(genes_annotations_filt) <- as.character(genes_annotations_filt$Entrezgene)
# sum(as.character(rownames(genes_annotations_filt)) %in% rownames(GSE198256_count))

##  Work with the annotated genes!
GSE198256_count_filt <- GSE198256_count[rownames(GSE198256_count) %in% rownames(genes_annotations_filt),]
GSE198256_count_exc <-GSE198256_count[!(rownames(GSE198256_count) %in% rownames(genes_annotations_filt)),]
annotgene_ord <- genes_annotations_filt[rownames(GSE198256_count_filt ),]
# sum(rownames(annotgene_ord)==rownames(GSE198256_count_filt))

#### QC: Exploratory  #### 

# Additional biological annotation - NOISeq 
library(NOISeq)
data_NOISEQ <- readData(data = GSE198256_count_filt,
                        length=abs(annotgene_ord$end-annotgene_ord$start),
                        gc=annotgene_ord$GC,
                        biotype= annotgene_ord$type ,
                        chromosome = annotgene_ord[,c("Chromosome","start","end")],
                        factors = Factors_GSE198256)


lengthuse <- abs(annotgene_ord$end-annotgene_ord$start)
names(lengthuse) <- rownames(annotgene_ord)
gc <- annotgene_ord$GC
names(gc) <- rownames(annotgene_ord)
biotype <-annotgene_ord$type
names(biotype) <- rownames(annotgene_ord)
chromosome <- annotgene_ord[,c("Chromosome","start","end")]

data_NOISEQ <- readData(data = GSE198256_count_filt,
                        length=lengthuse,
                        gc=gc,
                        biotype= biotype ,
                        chromosome = annotgene_ord[,c("Chromosome","start","end")],
                        factors = Factors_GSE198256)

# Length bias
mylengthbias = dat(data_NOISEQ, factor = "Group", type = "lengthbias")
explo.plot(mylengthbias, samples = NULL, toplot = "global")
# GC content bias 
myGCbias = dat(data_NOISEQ, factor = "Group", type = "GCbias")
explo.plot(myGCbias, samples = NULL, toplot = "global")
# RNA composition 
mycd = dat(data_NOISEQ, type = "cd", norm = FALSE, refColumn = 1)
explo.plot(mycd, samples = "GSM5942351")



# Creating a DESeq object and import the cpint data and sample information 









#### Normalization -- DESeq #### https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html #### 
library(DESeq2)

# creating DESeq object 
dds <- DESeqDataSetFromMatrix(countData = GSE198256_count_filt, colData = Factors_GSE198256, design = ~ Group)
View(counts(dds))

colData(dds)
sum(assay(dds)[,1])
barplot(log(colSums(assay(dds))))

is_expressed <- assay(dds) >= 10
df <- data.frame(Expressed = rowSums(is_expressed))
ggplot(df, aes(x=Expressed)) + geom_bar() + labs(x = "Total # samples Expressed (count >= 10)", y= "# genes")
keep <- rowSums(assay(dds) >= 10) >= 2
dds_f <- dds[keep,]
table(keep)

# inspection of Normalized data 
librarySizes <- colSums(assay(dds))
barplot(librarySizes, 
        names=names(librarySizes), 
        las=2, 
        main="Barplot of library sizes")
abline(h=20e6, lty=2)

# Get log2 counts per million
logcounts <- log2(assay(dds_f) + 1)
boxplot(logcounts, 
        xlab="", 
        ylab="Log2(Unormalizaed_Counts)",
        las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(as.matrix(logcounts)), col="blue")



#### the mean versus variance in read count data
mean_readCounts <- apply(assay(dds_f)[,1:3], 1, mean)
var_readCounts <- apply(assay(dds_f)[,1:3], 1, var)
df <- data.frame(mean_readCounts, var_readCounts)
ggplot(df) +
  geom_point(aes(x=mean_readCounts, y= var_readCounts)) +
  scale_y_log10() +
  scale_x_log10() +
  xlab("Mean counts per gene") +
  ylab("Variance per gene") +
  labs(title = "DESeq2 model - Dispersion")


mean_readCounts <- apply(GSE198256_count_filt[,1:3], 1, mean)
var_readCounts <- apply(GSE198256_count_filt[,1:3], 1, var)
df <- data.frame(mean_readCounts, var_readCounts)
ggplot(df) +
  geom_point(aes(x=mean_readCounts, y= var_readCounts)) +
  scale_y_log10() +
  scale_x_log10() +
  xlab("Mean counts per gene") +
  ylab("Variance per gene") +
  labs(title = "DESeq2 model - Dispersion")



#### PCA 
vsd <- vst(dds_f,blind=TRUE)
plotPCA(vsd,intgroup=c("Group"))

#### DE 
GSE198256_DESeq <- DESeq(dds_f)
GSE198256_DESeq_unfilt <-  DESeq(dds)
GSE198256_res <- results(GSE198256_DESeq, contrast = c("Group", "Covid19AI", "Covid193Mo"))
GSE198256_res_unfilt <- results(GSE198256_DESeq_unfilt, contrast = c("Group", "Covid19AI", "Covid193Mo"))
resultsNames(GSE198256_DESeq)


ggplot(as.data.frame(GSE198256_res_unfilt), aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.05, "red", "black")), alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  labs(title = "Volcano Plot for DEG of unfiltered counts",
       x = "Log2 Fold Change",
       y = "-log10(Adjusted p-value)")

EnhancedVolcano(as.data.frame(GSE198256_res),
                lab = rownames(as.data.frame(GSE198256_res)),
                x = 'log2FoldChange',
                y = 'pvalue')

#### dispersion
plotDispEsts(GSE198256_DESeq, main="Dispersion plot")

#### Shrinkage  
# before 
plotMA(GSE198256_res, ylim=c(-2,2))
# after
plotMA(lfcShrink(GSE198256_DESeq,coef=c("Group_Covid196Mo_vs_Covid193Mo")), ylim=c(-2,2))

resLFC <- lfcShrink(GSE198256_DESeq, res =GSE198256_res, coef="Group_Healthy_vs_Covid193Mo", type="apeglm")

#### Test 
wald_dss <- DESeq(dds_f, test = "Wald")
wald_dss_res <- results(wald_dss, contrast = c("Group", "Covid19AI", "Healthy"))
LRT_dss <- DESeq(dds_f, test="LRT", reduced=~1)
LRT_dss_res <- results(LRT_dss, contrast = c("Group", "Covid19AI", "Healthy"))

wald_dss <- DESeq(dds_f, test = "Wald")
wald_dssRes <- results(wald_dss, alpha = 0.05, lfcThreshold = 1)
summary(wald_dssRes)
significant_wald_dss <- subset(wald_dssRes, padj < 0.05)

#### comparisons of test 
summary(LRT_dss_res, alpha = 0.05, lfcThreshold = 2)
summary(wald_dss_res, alpha = 0.05, lfcThreshold = 2)
# additional
summary(results(wald_dss, alpha = 0.05, lfcThreshold = 0.58))
wald_dss <- DESeq(dds_f, test = "Wald")

#### acute vs early  
dds_a <- DESeq(dds_f)
contrast <- c("Group", "Covid19AI", "Covid193Mo")
result_acute_vs_early <- results(dds_a, contrast = contrast)
resLFC_acute_vs_early <- lfcShrink(dds_a, res =result_acute_vs_early, coef="Group_Covid19AI_vs_Covid193Mo", type="apeglm")
resLFC_acute_vs_early_f <- subset(result_acute_vs_early, padj <0.05 & log2FoldChange >2 | log2FoldChange < -2)
acute_genes <- rownames(result_acute_vs_early_f)
  
  
#### acute vs healthy   
dds_a <- DESeq(dds_f)
contrast <- c("Group", "Covid19AI", "Healthy")
result_acute_vs_healthy <- results(dds_a, contrast = contrast)
result_acute_vs_healthy_f_up <- subset(result_acute_vs_healthy, padj <  0.05 & log2FoldChange > 2)
result_acute_vs_healthy_f_down <- subset(result_acute_vs_healthy, padj <  0.05 & log2FoldChange < -2)
acute_genes <- rownames(result_acute_vs_healthy_f_up)

library("ggVennDiagram")

x <- list(
  up = rownames(result_acute_vs_healthy_f_up), 
  down = rownames(result_acute_vs_healthy_f_down)
)
ggVennDiagram(x[1:2], label_alpha = 0) + coord_flip()



#### DE Boxplots 
plotCounts(dds_f, gene="3248", intgroup = "Group", returnData = TRUE) %>%
  ggplot() + aes(Group, count) + geom_boxplot(aes(fill=Group)) + scale_y_log10() + theme_bw() + ggtitle("HPGD")
plotCounts(dds_f, gene="3553", intgroup = "Group", returnData = TRUE) %>%
  ggplot() + aes(Group, count) + geom_boxplot(aes(fill=Group)) + scale_y_log10() + theme_bw() + ggtitle("IL1B")

plotCounts(dds_f, gene="2180", intgroup="Group")
plotCounts(dds_f, gene="3248", intgroup="Group")
plotCounts(dds_f, gene="3553", intgroup="Group")


#### DE heatmaps 
#1
library("pheatmap")
acute_genes
vsd <- assay(vst(dds_f))[acute_genes[1:200],]
Z <- t(scale(t(vsd)))
pheatmap(Z, cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)



#2
library(ComplexHeatmap)
result_acute_vs_healthy_sig <- subset(result_acute_vs_healthy, padj <  0.05 & log2FoldChange > 2 |  log2FoldChange < -2 )
acute_vs_healthy_df <- as.data.frame(result_acute_vs_healthy_sig)
mat <- counts(dds_a, normalized = T)[rownames(acute_vs_healthy_df), ]
mat.z <- t(apply(mat, 1, scale))
colnames(mat.z) <- rownames(Factors_GSE198256)
Heatmap(mat.z, cluster_rows = T, cluster_columns = F, column_labels = colnames(mat.z), name = "z-score")


##### GSEA
library(clusterProfiler)
library(enrichplot)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)


res_acuteVsearly_genes <- result_acute_vs_early$log2FoldChange
names(res_acuteVsearly_genes) <- rownames(result_acute_vs_early)
gene_list<-na.omit(res_acuteVsearly_genes)
gene_list = sort(gene_list, decreasing = TRUE)
gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "ENTREZID", 
             minGSSize = 200, 
             maxGSSize = 1000, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "BH")
# dotplot 
require(DOSE)
View(gse@result)
which(gse@result$Description == "killing of cells of another organism")
dotplot(gse, showCategory=15, split=".sign") + facet_grid(.~.sign)

gseaplot2(gse, gse@result$ID[159], title = gse@result$Description[159])

