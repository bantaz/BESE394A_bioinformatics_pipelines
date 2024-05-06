
# module load gsl
# BiocManager::install("TFBSTools")
# remotes::install_github('satijalab/azimuth', ref = 'master')

library(Seurat)
library(ggplot2)
library(Azimuth)

# GEX https://www.10xgenomics.com/datasets/5k-adult-mouse-brain-nuclei-isolated-with-chromium-nuclei-isolation-kit-3-1-standard
# ATAC https://www.10xgenomics.com/datasets/8k-adult-mouse-cortex-cells-atac-v2-chromium-controller-2-standard
# Multiome https://www.10xgenomics.com/datasets/mouse-brain-nuclei-isolated-with-chromium-nuclei-isolation-kit-saltyez-protocol-and-10x-complex-tissue-dp-ct-sorted-and-ct-unsorted-1-standard

plot_QC_features <- function(seurat_object, pdf_path){
    pdf(pdf_path)
        # shows the distribution of the transcripts per cells
        print(ggplot(seurat_object@meta.data,aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
            geom_density(alpha = 0.2) + 
            scale_x_log10() + 
            theme_classic() +
            ylab("Cell density") +
            geom_vline(xintercept = 500) + ggtitle('UMISpercell'))
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
            geom_hline(yintercept = 250) + ggtitle(paste0('UMIvsGENESpercell  Ncell:: ', ncol(mouse_brain))))
        print(VlnPlot(mouse_brain, features = c("nFeature_RNA", "nCount_RNA", "mitoPct", "RPSPct", 'log10GenesPerUMI'), ncol = 3, pt.size=0, layer='counts'))
    dev.off()
}

how_many_PCs <- function(obj, pdf_name){
  # determine the correct PC number
  # Determine percent of variation associated with each PC
  pct <- obj[["pca"]]@stdev / sum(obj[["pca"]]@stdev) * 100
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  # last point where change of % of variation is more than 0.1%.
  # Minimum of the two calculation
  pcs <- min(co1, co2)
  plot_df <- data.frame(pct = pct, 
           cumu = cumu, 
           rank = 1:length(pct))
  # Elbow plot to visualize 
  pdf(pdf_name)
  print(ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw())
  print(ElbowPlot(obj))
  dev.off()
  return(pcs)
}



mouse_brain <- Read10X(data.dir = '/ibex/user/serrang/Projects_data/MultiomeCourse/Data/5k_mouse_brain_GEX/filtered_feature_bc_matrix')
mouse_brain <- CreateSeuratObject(
  counts = mouse_brain,
  project = 'mouse_brain',
  assay = "RNA",
  min.cells = 3, 
  min.features = 200)

mouse_brain$mitoPct <- PercentageFeatureSet(mouse_brain, pattern = "^mt-")
mouse_brain$RPSPct  <- PercentageFeatureSet(object = mouse_brain, pattern = "^Rp[sl]")
mouse_brain$log10GenesPerUMI <- log10(mouse_brain$nFeature_RNA) / log10(mouse_brain$nCount_RNA)

pdf('/ibex/user/serrang/Projects_data/MultiomeCourse/Plots/5k_mouse_brain_GEX_QC.pdf')
VlnPlot(mouse_brain, features = c("nFeature_RNA", "nCount_RNA", "mitoPct"), ncol = 3)
plot1 <- FeatureScatter(mouse_brain, feature1 = "nCount_RNA", feature2 = "mitoPct")
plot2 <- FeatureScatter(mouse_brain, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(plot1 + plot2)
dev.off()


plot_QC_features(mouse_brain, '/ibex/user/serrang/Projects_data/MultiomeCourse/Plots/5k_mouse_brain_GEX_QC_Pre.pdf')

# three to five times of standard deviation or median absolute deviation from the median
mouse_brain <- subset(mouse_brain, 
                      subset =  nFeature_RNA > 400 & 
                                nFeature_RNA < 6000 & 
                                RPSPct < 10 &
                                mitoPct < sd(mouse_brain$mitoPct) *3 & 
                                log10GenesPerUMI > 0.8)

plot_QC_features(mouse_brain, '/ibex/user/serrang/Projects_data/MultiomeCourse/Plots/5k_mouse_brain_GEX_QC_Post.pdf')

# preprocessing
mouse_brain <- NormalizeData(mouse_brain)
mouse_brain <- FindVariableFeatures(mouse_brain, nfeatures = 3000)
mouse_brain <- ScaleData(mouse_brain)
mouse_brain <- RunPCA(mouse_brain, npcs = 100)
mouse_brain <- RunTSNE(mouse_brain, dims = 1:30)
PCA_dims <- how_many_PCs(mouse_brain, '/ibex/user/serrang/Projects_data/MultiomeCourse/Plots/Num_of_PCs.pdf')
mouse_brain <- FindNeighbors(mouse_brain, dims = 1:30)
mouse_brain <- FindClusters(mouse_brain, resolution = 0.4, algorithm = 3)
mouse_brain <- RunUMAP(mouse_brain, dims = 1:30)

pdf('/ibex/user/serrang/Projects_data/MultiomeCourse/Plots/5k_mouse_brain_GEX_QC_DimRed.pdf')
DimPlot(mouse_brain, reduction = "umap", pt.size = 0.9, label = TRUE, label.size = 3) + NoLegend()
DimPlot(mouse_brain, reduction = "tsne", pt.size = 0.9, label = TRUE, label.size = 3) + NoLegend()
dev.off()


mouse_brain <- Azimuth::RunAzimuth(mouse_brain, reference = "mousecortexref")
pdf('/ibex/user/serrang/Projects_data/MultiomeCourse/Plots/5k_mouse_brain_GEX_Clusters.pdf')
DimPlot(mouse_brain, reduction = "umap", pt.size = 0.9, label = TRUE, label.size = 3) + NoLegend()
DimPlot(mouse_brain, reduction = "umap", group.by = "predicted.subclass", pt.size = 0.9, label = TRUE, label.size = 3) + NoLegend()
dev.off()



saveRDS(mouse_brain, '/ibex/user/serrang/Projects_data/MultiomeCourse/Data/5k_mouse_brain_GEX.rds')