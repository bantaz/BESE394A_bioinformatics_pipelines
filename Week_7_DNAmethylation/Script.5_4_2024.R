##################################
#Runing shinyepico user interface
##################################
library("shinyepico")
run_shinyepico()

#############
#Runing ChAMP
#############
library(ChAMP)
#load the data & process raw signals...
myload <- champ.load(directory = "/home/sobahytm/class_7/GSE188573/idat",
    method="ChAMP",
    methValue="B",
    autoimpute=TRUE,
    filterDetP=TRUE,
    ProbeCutoff=0,
    SampleCutoff=0.1,
    detPcut=0.01,
    filterBeads=TRUE,
    beadCutoff=0.05,
    filterNoCG=TRUE,
    filterSNPs=TRUE,
    population=NULL,
    filterMultiHit=TRUE,
    filterXY=TRUE,
    force=FALSE,
    arraytype="EPIC")
# runc QC...
champ.QC(beta = myload$beta,
    pheno=myload$pd$Sample_Group,
    mdsPlot=TRUE,
    densityPlot=TRUE,
    dendrogram=TRUE,
    PDFplot=TRUE,
    Rplot=TRUE,
    Feature.sel="None",
    resultsDir="/home/sobahytm/class_7/GSE188573/CHAMP_QCimages/")
# normalize...
mynorm <- champ.norm(beta=myload$beta,
    rgSet=myload$rgSet,
    mset=myload$mset,
    resultsDir="/home/sobahytm/class_7/GSE188573/CHAMP_Normalization/",
    method="BMIQ",
    plotBMIQ=FALSE,
    arraytype="EPIC",
    cores=3)
# detect DMP...
myDMP <- champ.DMP(beta = mynorm,
    pheno = myload$pd$Sample_Group,
    compare.group = NULL,
    adjPVal = 0.05,
    adjust.method = "BH",
    arraytype = "EPIC")
# save ...
capture.output(myDMP, file = "/home/sobahytm/class_7/GSE188573/CHAMP_DMP.txt") 
# detect DMR...
myDMR <- champ.DMR(beta=mynorm,
    pheno=myload$pd$Sample_Group,
    compare.group=NULL,
    arraytype="EPIC",
    method = "Bumphunter",
    minProbes=7,
    adjPvalDmr=0.05,
    cores=3,
    maxGap=300,
    cutoff=NULL,
    pickCutoff=TRUE,
    smooth=TRUE,
    smoothFunction=loessByCluster,
    useWeights=FALSE,
    permutations=NULL,
    B=250,
    nullMethod="bootstrap")
# save ...
capture.output(myDMR, file = "/home/sobahytm/class_7/GSE188573/CHAMP_DMR.txt")
# run GSEA ...

myGSEA_DMR <- champ.GSEA(beta=mynorm,
    DMR=myDMR,
    CpGlist=NULL,
    Genelist=NULL,
    pheno=myload$pd$Sample_Group,
    method="fisher",
    arraytype="EPIC",
    Rplot=TRUE,
    adjPval=0.05,
    cores=1)

# save ...
capture.output(myGSEA_DMR, file = "/home/sobahytm/class_7/GSE188573/CHAMP_GSEA_DMR.txt")
