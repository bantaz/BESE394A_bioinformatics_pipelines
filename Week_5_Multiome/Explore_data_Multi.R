# module load gsl

options(width = 210)
library(Seurat)
library(TFBSTools)
library(motifmatchr)
library(Signac)
library(ggplot2)
library(Azimuth)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(JASPAR2020)
# Multiome https://www.10xgenomics.com/datasets/mouse-brain-nuclei-isolated-with-chromium-nuclei-isolation-kit-saltyez-protocol-and-10x-complex-tissue-dp-ct-sorted-and-ct-unsorted-1-standard
# wget https://cf.10xgenomics.com/samples/cell-arc/2.0.2/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_filtered_feature_bc_matrix.tar.gz
# wget https://cf.10xgenomics.com/samples/cell-arc/2.0.2/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_atac_fragments.tsv.gz
# wget https://cf.10xgenomics.com/samples/cell-arc/2.0.2/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_atac_fragments.tsv.gz.tbi

plot_QC_features <- function(seurat_object, pdf_path){
  pdf(pdf_path, width =10)
      # shows the distribution of the transcripts per cells
      print(VlnPlot(
          object = seurat_object,
          features = c("nCount_RNA", "nCount_ATAC", 'nFeature_RNA',  'mitoPct', 
                      'FRiP', 'log10GenesPerUMI', "TSS.enrichment", 
                      "nucleosome_signal"),
          pt.size = 0, ncol =4))
        DefaultAssay(seurat_object) <- "ATAC"
        print(FragmentHistogram(object = seurat_object, region = 'chr1-1-10000000', group.by = 'nucleosome_group'))
        print(TSSPlot(seurat_object, group.by = 'high.tss') + NoLegend())
        print(DensityScatter(seurat_object, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE))
      print(ggplot(seurat_object@meta.data, aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
          geom_density(alpha = 0.2) + 
          theme_classic() +
          scale_x_log10() + 
          geom_vline(xintercept = 300) + ggtitle('GENESpercell'))
      # Visualize the distribution of genes detected per cell via boxplot
      print(ggplot(seurat_object@meta.data, aes(x=orig.ident, y=log10(nFeature_RNA), fill=orig.ident)) + 
          geom_boxplot() + 
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
          theme(plot.title = element_text(hjust=0.5, face="bold")) +
          ggtitle("NCells vs NGenes"))
      # correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
      print(ggplot(seurat_object@meta.data, aes(x=nCount_RNA, y=nFeature_RNA, colour=mitoPct, group=orig.ident)) + 
          geom_point() + 
          scale_colour_gradient(low = "gray90", high = "black", limits=c(0,100)) +
          stat_smooth(method=lm) +
          scale_x_log10() + 
          scale_y_log10() + 
          theme_classic() +
          geom_vline(xintercept = 500) +
      geom_hline(yintercept = 6000) +
          geom_hline(yintercept = 250) + ggtitle(paste0('UMIvsGENESpercell  Ncell:: ', ncol(seurat_object))))
  dev.off()

}

### Loading and storing scRNA-seq data from a directory. The data is filtered feature-barcode matrices from a 10X Genomics experiment
counts <- Read10X("/Users/bantanai/Library/Mobile Documents/com~apple~CloudDocs/BESE394/Week_5/Data/filtered_feature_bc_matrix")
### A directory path to the ATAC-seq fragments file which contains information on accessible chromatin regions 
fragpath <- '/Users/bantanai/Library/Mobile Documents/com~apple~CloudDocs/BESE394/Week_5/Data/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_atac_fragments.tsv.gz'
# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

# create a Seurat object containing the RNA adata
mouse_brain <- CreateSeuratObject(
  counts = counts[['Gene Expression']],
  assay = "RNA"
)

# create ATAC assay and add it to the object
# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
atac_counts <- counts$Peaks
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]

mouse_brain[["ATAC"]] <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

mouse_brain[["ATAC"]]@counts


################################################################################
# Group assignment: Group 5 
################################################################################

# MACS2 has been downloaded locally and the path to tool has been defined
# $ pip install MACS2
MASC2_path <- "/Users/bantanai/opt/anaconda3/bin/MACS2"

# calling peaks using MACS2 
peaks_macs2 <- CallPeaks(
  object = mouse_brain[["ATAC"]],
  macs2.path = MASC2_path)
 # saveRDS(peaks, file ="/Users/bantanai/Library/Mobile Documents/com~apple~CloudDocs/BESE394/Week_5/Data/peaks.rds")

# Visualizing MACS2 peak calls alongside the 10x Cellranger peak calls of two chosen regions. Cellranger peaks are shown in grey and the MACS2 peaks in red. 
DefaultAssay(mouse_brain) <- "ATAC"
CoveragePlot(
  object = mouse_brain,
  region = "chr9-102498973-102499900",
  ranges = peaks_macs2,
  ranges.title = "MACS2"
)
CoveragePlot(
  object = mouse_brain,
  region = "chr1-3229745-3230583",
  ranges = peaks_macs2,
  ranges.title = "MACS2"
)


mouse_brain_2 <- mouse_brain
# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(mouse_brain_2),
  features = peaks,
  cells = colnames(mouse_brain_2)
)
## add precomputed feature matrix to the seurat object:
mouse_brain_2 <- mouse_brain_2[,colnames(mouse_brain_2) %in% colnames(macs2_counts)]
macs2_counts <- macs2_counts[,colnames(mouse_brain_2)]
# remove previous ATAC assay 
DefaultAssay(mouse_brain_2) <- "RNA"
mouse_brain_2[["ATAC"]] <- NULL
# create a new assay using the MACS2 peak set and add it to the Seurat object
mouse_brain_2[["ATAC"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments  = fragpath,
  annotation = annotation
)
DefaultAssay(mouse_brain_2) <- "ATAC"
mouse_brain_2 <- NucleosomeSignal(mouse_brain_2)
mouse_brain_2 <- TSSEnrichment(mouse_brain_2, fast=FALSE)
total_fragments <- CountFragments(fragments = fragpath)
rownames(total_fragments) <- total_fragments$CB
mouse_brain_2$fragments <- total_fragments[colnames(mouse_brain_2), "frequency_count"]
mouse_brain_2 <- FRiP(object = mouse_brain_2, assay = 'ATAC', total.fragments = 'fragments')
mouse_brain_2$nucleosome_group <- ifelse(mouse_brain_2$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
mouse_brain_2$high.tss <- ifelse(mouse_brain_2$TSS.enrichment > 3, 'High', 'Low')

DefaultAssay(mouse_brain_2) <- "RNA"
mouse_brain_2$mitoPct <- PercentageFeatureSet(mouse_brain_2, pattern = "^mt-")
mouse_brain_2$RPSPct  <- PercentageFeatureSet(object = mouse_brain_2, pattern = "^Rp[sl]")
mouse_brain_2$log10GenesPerUMI <- log10(mouse_brain_2$nFeature_RNA) / log10(mouse_brain_2$nCount_RNA)

plot_QC_features(mouse_brain_2, '/Users/bantanai/Library/Mobile Documents/com~apple~CloudDocs/BESE394/Week_5/Plots/mouse_brain_2_QC.pdf')

## Comparing the number of peaks  
DefaultAssay(mouse_brain) <- "ATAC"
DefaultAssay(mouse_brain_2) <- "ATAC"
N_peaks_cellrangerArc <- nrow(mouse_brain)
N_peaks_macs2 <- nrow(mouse_brain_2)
df <- data.frame(
  Method = c("Cell Ranger ARC", "MACS2"),
  Number_of_Peaks = c(N_peaks_cellrangerArc, N_peaks_macs2)
)
ggplot(df, aes(x = Method, y = Number_of_Peaks, fill = Method)) +
  geom_bar(stat = "identity", color = "black") +
  theme_minimal() +
  labs(title = "Comparison of Peak Counts",
       x = "Method",
       y = "Number of Peaks") +
  scale_fill_manual(values = c("Cell Ranger ARC" = "skyblue", "MACS2" = "salmon")) +
  geom_text(aes(label = Number_of_Peaks), vjust = -0.3, size = 6)


## filtering 
mouse_brain_2_filt <- subset(
  x = mouse_brain_2,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1 & 
    FRiP > 0.2
)
plot_QC_features(mouse_brain_2_filt, '/Users/bantanai/Library/Mobile Documents/com~apple~CloudDocs/BESE394/Week_5/Plots/Mouse_brain_2_QC_Post.pdf')


################################################################################
# Calculate metrics for QC 
################################################################################
# here, the QC metrices are computed separately for each assay (RNA and ATAC). QC filtering will also applied separately. 

DefaultAssay(mouse_brain) <- "ATAC" # to compute the QC matrices for the ATAC assay 
mouse_brain <- NucleosomeSignal(mouse_brain) # to identify nucleosome-free regions versus mononucleosome or dinucleosome.
mouse_brain <- TSSEnrichment(mouse_brain, fast=FALSE) # to assess how well the ATAC-seq signal is enriched at TSS
total_fragments <- CountFragments(fragments = fragpath) # calculating the total number of ATAC-seq fragments, which will be used for calculating the fraction of reads in peaks (FRiP) score. 
rownames(total_fragments) <- total_fragments$CB
mouse_brain$fragments <- total_fragments[colnames(mouse_brain), "frequency_count"]
mouse_brain <- FRiP(object = mouse_brain, assay = 'ATAC', total.fragments = 'fragments') # the fraction of sequenced reads that fall into the called peaks
# categorizinh cells into two groups based on their nucleosome signal. If a cell's nucleosome signal is greater than 4, it is assigned to the 'NS > 4' group; otherwise, it is assigned to the 'NS < 4' group
mouse_brain$nucleosome_group <- ifelse(mouse_brain$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
mouse_brain$high.tss <- ifelse(mouse_brain$TSS.enrichment > 3, 'High', 'Low')
# categorizing cells based on their TSS enrichment score. If a cell's score is greater than 3, it is labeled as having 'High' TSS enrichment, otherwise 'Low' -- will be used in the QC filtering step 


DefaultAssay(mouse_brain) <- "RNA" # to compute the QC matrices for the RNA assay 
mouse_brain$mitoPct <- PercentageFeatureSet(mouse_brain, pattern = "^mt-") # compute % of mitochondrial genes 
mouse_brain$RPSPct  <- PercentageFeatureSet(object = mouse_brain, pattern = "^Rp[sl]") # compute % of ribosomal genes  
mouse_brain$log10GenesPerUMI <- log10(mouse_brain$nFeature_RNA) / log10(mouse_brain$nCount_RNA) # computing log10GenesPerUMI reflect on the complexity 

# plotting QC and saving the plots in a pdf file 
plot_QC_features(mouse_brain, '/Users/bantanai/Library/Mobile Documents/com~apple~CloudDocs/BESE394/Week_5/Plots/5k_mouse_brain_GEX_QC_Post.pdf')




################################################################################
# Filter cells based on QC metrics
################################################################################
# filtering cells (considering both RNA and ATAC assays)
mouse_brain <- subset(
  x = mouse_brain,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1 & 
    FRiP > 0.2
)

saveRDS(mouse_brain, '/Users/bantanai/Library/Mobile Documents/com~apple~CloudDocs/BESE394/Week_5/Data/mouse_brain_multiome_Post_QC.rds')
# mouse_brain <- readRDS('/ibex/user/serrang/Projects_data/MultiomeCourse/Data/mouse_brain_multiome_Post_QC.rds')
plot_QC_features(mouse_brain, '/Users/bantanai/Library/Mobile Documents/com~apple~CloudDocs/BESE394/Week_5/Plots/5k_mouse_brain_GEX_QC_Post.pdf')



################################################################################
# Process the data (normalization, scaling, PCA, clustering, UMAP, ....)
################################################################################

DefaultAssay(mouse_brain) <- "RNA"
# normalization step to adjusts raw gene expression counts for making fair comparisons between cells
mouse_brain <- NormalizeData(mouse_brain) 
#  Identifies the most variably expressed genes across cells using a variance-stabilizing transformation (VST). This focuses subsequent analyses on genes that are most informative about cellular heterogeneity.
mouse_brain <- FindVariableFeatures(mouse_brain, 
            selection.method = "vst", 
            nfeatures = 2000)
# Scales and centers the expression data for each gene which is necessary for  dimensionality reduction
mouse_brain <- ScaleData(mouse_brain)
# Calculates a K-nearest neighbors graph based on the cells' expression profiles in reduced dimensions, which lays the foundation for clustering and UMAP
mouse_brain <- FindNeighbors(mouse_brain, dims = 1:30)
mouse_brain <- FindClusters(mouse_brain, 
            resolution = 0.4, 
            algorithm = 3, 
            cluster.name="RNA_clusters_03")
mouse_brain <- RunPCA(mouse_brain)
mouse_brain <- RunUMAP(mouse_brain, dims=1:30, 
            reduction = "pca", 
            reduction.name = "rna_umap")
mouse_brain <- Azimuth::RunAzimuth(mouse_brain, 
            reference = "mousecortexref")
# plotting 
pdf('/Users/bantanai/Library/Mobile Documents/com~apple~CloudDocs/BESE394/Week_5/Plots/UMAP_RNA.pdf')
DimPlot(mouse_brain, reduction = "rna_umap", group.by='RNA_clusters_03', 
        label = TRUE, repel = TRUE) + NoLegend()
DimPlot(mouse_brain, reduction = "rna_umap", group.by = "predicted.subclass", 
        label = TRUE, label.size = 3) + NoLegend()
dev.off()


DefaultAssay(mouse_brain) <- "ATAC"
# Identifies the most accessible regions (peaks) across the genome by selecting features (peaks) that appear in at least min.cutoff number of cells.
mouse_brain <- FindTopFeatures(mouse_brain, min.cutoff = 5)
#  Transforms raw accessibility counts into a normalized measure 
mouse_brain <- RunTFIDF(mouse_brain)
# Performs a dimensionality reduction on the TF-IDF-transformed data
mouse_brain <- RunSVD(mouse_brain)

pdf('/Users/bantanai/Library/Mobile Documents/com~apple~CloudDocs/BESE394/Week_5/Plots/Correlation.pdf')
DepthCor(mouse_brain)
dev.off()

mouse_brain <- RunUMAP(mouse_brain, dims=1:30, 
                reduction = "lsi", 
                reduction.name = "atac_umap")

pdf('/Users/bantanai/Library/Mobile Documents/com~apple~CloudDocs/BESE394/Week_5/Plots/UMAP_ATAC.pdf')

DimPlot(mouse_brain, reduction = "atac_umap", label = TRUE, repel = TRUE) + NoLegend()
DimPlot(mouse_brain, reduction = "atac_umap", group.by = "predicted.subclass",
label = TRUE, label.size = 3) + NoLegend() # labeling the clusters with predicted subclass annotation. Such umap helps to identify the hetergenity and/or similar clusters, and identify patterns in high-dimensional data
dev.off()

saveRDS(mouse_brain, '/Users/bantanai/Library/Mobile Documents/com~apple~CloudDocs/BESE394/Week_5/Data/mouse_brain_multiome.rds')

################################################################################
# Integrate the RNA and ATAC data
################################################################################

# integration is a crucial step, it helps linking the genomic landscape and transcriptomic activity, providing a multidimensional view of the cellular machinery
# a weighted nearest neighbor (WNN) graph is built, integrating the information from both the PCA of RNA-seq data and the LSI of ATAC-seq data to find neighbors
mouse_brain <- FindMultiModalNeighbors(
  object = mouse_brain,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 1:40),
  verbose = TRUE
)

# build a joint (WNN) UMAP visualization
mouse_brain <- RunUMAP(
  object = mouse_brain, 
  reduction.name = "wnn.umap", dims = 1:30,
  assay = "RNA",
  verbose = TRUE
)

saveRDS(mouse_brain, '/Users/bantanai/Library/Mobile Documents/com~apple~CloudDocs/BESE394/Week_5/Data/mouse_brain_multiome_Integrated.rds')

pdf('/Users/bantanai/Library/Mobile Documents/com~apple~CloudDocs/BESE394/Week_5/Plots/Umap_integrated.pdf')
DimPlot(mouse_brain, reduction = "rna_umap", label = TRUE, group.by = "predicted.subclass",
        repel = TRUE) + NoLegend()
DimPlot(mouse_brain, reduction = "atac_umap", label = TRUE, group.by = "predicted.subclass",
        repel = TRUE) + NoLegend()
DimPlot(mouse_brain, reduction = "wnn.umap", label = TRUE, group.by = "predicted.subclass",
        repel = TRUE) + NoLegend()
dev.off()

# Here, there's a balance between RNA and ATAC contributions, but higher contribution might be slightly coming from the RNA-seq data 
pdf('/Users/bantanai/Library/Mobile Documents/com~apple~CloudDocs/BESE394/Week_5/Plots/Umap_integrated_weights.pdf')
Idents(mouse_brain) <- "predicted.subclass"
VlnPlot(mouse_brain, features = c("RNA.weight", "ATAC.weight"),  pt.size = 0, ncol = 1) 
dev.off()

################################################################################
# Perform a Differential accessibility analysis
################################################################################

# finding differential accessible peaks between two groups of cells, Oligo and Astro using likelihood ratio (LR) test  
DefaultAssay(mouse_brain) <- 'ATAC'
Idents(mouse_brain) <- "predicted.subclass"
da_peaks <- FindMarkers(
  object = mouse_brain,
  ident.1 = c("Oligo"), 
  ident.2 = c("Astro"),
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)
saveRDS(da_peaks, '/Users/bantanai/Library/Mobile Documents/com~apple~CloudDocs/BESE394/Week_5/Data/DA_peaks.rds')

# applying sorting from most positive (more accessible differentially expressed peaks) to most negative values of average log2 fold change 
da_peaks <- da_peaks[order(da_peaks$avg_log2FC, decreasing = TRUE), ]

pdf('/Users/bantanai/Library/Mobile Documents/com~apple~CloudDocs/BESE394/Week_5/Plots/Peaks_DA.pdf', width=12, height=20)
cowplot::plot_grid(
  VlnPlot(
    object = mouse_brain,
    assay = 'ATAC',
    features = rownames(da_peaks)[1:3],
    pt.size = 0.1,
    group.by='predicted.subclass', 
    ncol=3
  ),
  VlnPlot(
    object = mouse_brain,
    assay = 'RNA',
    features = ClosestFeature(mouse_brain, rownames(da_peaks)[1:3])$gene_name,
    pt.size = 0.1,
    group.by='predicted.subclass', 
    ncol=3
  ),
  FeaturePlot(
    object = mouse_brain,
    reduction = 'wnn.umap', 
    order=TRUE,
    features = rownames(da_peaks)[1:3],
    pt.size = 0.1,
    max.cutoff = 'q95',
    ncol=3
  ) & NoLegend(),
  FeaturePlot(
    object = mouse_brain,
    reduction = 'wnn.umap', 
    order=TRUE,
    features = ClosestFeature(mouse_brain, rownames(da_peaks)[1:3])$gene_name,
    pt.size = 0.1,
    max.cutoff = 'q95',
    ncol=3
  ) & NoLegend(),
nrow=4)
dev.off()


# Annotate peaks with his closest feature
open_Oligo <- rownames(da_peaks[da_peaks$avg_log2FC > 3, ])
open_Astro <- rownames(da_peaks[da_peaks$avg_log2FC < 3, ])
closest_Oligo <- ClosestFeature(mouse_brain, open_Oligo)
closest_Astro <- ClosestFeature(mouse_brain, open_Astro)

# https://www.cellsignal.com/pathways/neuronal-and-glial-cell-markers
# Visualize the coverage of the peaks
pdf('/Users/bantanai/Library/Mobile Documents/com~apple~CloudDocs/BESE394/Week_5/Plots/Coverage_selected.pdf', height=12)
CoveragePlot(
  object = mouse_brain,
  region = c('Olig1', 
             'Gfap'),
  extend.upstream = 1000,
  extend.downstream = 1000,
  ncol = 1
)
dev.off()



################################################################################
# Motif analysis
################################################################################

# Motif analysis with the DA peaks using JASPAR database, specifying for the CORE collection which contains a curated set of high-quality, non-redundant transcription factor binding profiles
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
mouse_brain <- AddMotifs(
  object = mouse_brain,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)


top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])

DefaultAssay(mouse_brain) <- "ATAC"
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 'Mus musculus', all_versions = FALSE))
motif.matrix <- CreateMotifMatrix(features = granges(mouse_brain), pwm = pwm_set, genome = 'mm10', use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
mouse_brain <- SetAssayData(mouse_brain, assay = 'ATAC', layer = 'motifs', new.data = motif.object)

# Note that this step can take 30-60 minutes 
mouse_brain <- RunChromVAR(
  object = mouse_brain,
  genome = BSgenome.Mmusculus.UCSC.mm10
)


enriched.motifs <- FindMotifs(
  object = mouse_brain,
  features = top.da.peak
)

pdf('/Users/bantanai/Library/Mobile Documents/com~apple~CloudDocs/BESE394/Week_5/Plots/Motif_enrichment.pdf')
MotifPlot(
  object = mouse_brain,
  motifs = head(rownames(enriched.motifs))
)
dev.off()

################################################################################
# Generate a RNA activity matrix based on the ATAC-seq data
################################################################################

# We can create a proxy of the gene expression from the ATAC-seq data using  GeneActivity function. 
# We can also create a proxy of the gene expression from the ATAC-seq data using the chromVAR package. 
# This package uses the motif accessibility to infer the gene expression. 
# We can use the motif models from JASPAR2020 to perform this analysis.

gene.activities <- GeneActivity(mouse_brain)
mouse_brain[['RNA_ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
mouse_brain <- NormalizeData(
  object = mouse_brain,
  assay = 'RNA_ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(mouse_brain$nCount_RNA)
)

DefaultAssay(mouse_brain) <- 'RNA'
RNA_plot <- FeaturePlot(
    object = mouse_brain,
    order=TRUE,
    features =c("Sema5a","Dennd4a","Nkain1"),
    pt.size = 0.1,
    max.cutoff = 'q95',
    ncol = 3
  )& NoLegend()

DefaultAssay(mouse_brain) <- 'RNA_ACTIVITY'
RNA_Activity_plot <- FeaturePlot(
    object = mouse_brain,
    order=TRUE,
    features =c("Sema5a","Dennd4a","Nkain1"),
    pt.size = 0.1,
    max.cutoff = 'q95',
    ncol = 3
  )& NoLegend()

pdf('/Users/bantanai/Library/Mobile Documents/com~apple~CloudDocs/BESE394/Week_5/Plots/RNA_comparison.pdf', height=12)
cowplot::plot_grid(
  RNA_plot,
  RNA_Activity_plot,
nrow=2)
dev.off()
