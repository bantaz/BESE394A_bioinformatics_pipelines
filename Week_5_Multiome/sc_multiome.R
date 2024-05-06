# libraries 
library(Signac)
library(Seurat)
library(SeuratObject)
library(Azimuth)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
library(chromVAR)
library(motifmatchr)



install.packages('devtools')
devtools::install_github('stuart-lab/signac', ref = 'develop')
install.packages('Signac')
BiocManager::install("TFBSTools")
BiocManager::install("Azimuth")
BiocManager::install('BSgenome.Mmusculus.UCSC.mm10')
BiocManager::install('EnsDb.Mmusculus.v79')
BiocManager::install('chromVAR')
BiocManager::install('motifmatchr')
remotes::install_github('satijalab/azimuth', ref = 'master')

