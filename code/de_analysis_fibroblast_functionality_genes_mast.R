# Load the relevant packages
library(dplyr)
library(patchwork)
library(ggplot2)
library(SummarizedExperiment)
library(zellkonverter)
library(SingleCellExperiment)
library(MAST)
library(EnhancedVolcano)

fibroblast_sce <- readH5AD("../data/integratedfibroblastsdata.h5ad")

# Load the different gene lists 
fib_signal_genes <- read.csv("../data/fibroblast_signal_ligands.csv", header=T)$X0
fib_ecm_structural_genes <- read.csv("../datafibroblast_ecm_structural_genes.csv", header=T)$X0


# Subset for the various different gene lists
fib_signal_sce <- fib_sce[fib_signal_genes, ]
fib_ecm_structural_sce <- fib_sce[fib_ecm_structural_genes, ]

# Convert to SCA objects
fib_signal_sca <- SceToSingleCellAssay(fib_signal_sce, class="SingleCellAssay")
fib_ecm_structural_sca <-  SceToSingleCellAssay(fib_ecm_structural_sce, class="SingleCellAssay")

colData(fib_signal_sca)$n_genes <- scale(colData(fib_signal_sca)$n_genes)
colData(fib_ecm_structural_sca)$n_genes <- scale(colData(fib_ecm_structural_sca)$n_genes)

# Filter out genes with zero expression
fib_signal_sca_filt <- fib_signal_sca[rowSums(assay(fib_signal_sca)) != 0, ]
fib_ecm_structural_sca_filt <- fib_ecm_structural_sca[rowSums(assay(fib_ecm_structural_sca)) != 0, ]

# Relevel factors so that the unwounded condition and mainly unwounded cluster, FIB-III, and the lower cells are the reference levels
levels(colData(fib_signal_sca_filt)$sample) <- c("LW_FIB_PWD18","LW_PWD12", "LW_PWD14","LW_REG_PWD18", "SW_PWD4", "SW_PWD7", "UW_P49", "UW_P21")
cond <- factor(colData(fib_signal_sca_filt)$sample)
cond <- relevel(cond,"UW_P49")
colData(fib_signal_sca_filt)$sample <- cond
clust <- factor(colData(fib_signal_sca_filt)$leiden_sub)
clust <- relevel(clust, "FIB-III")
colData(fib_signal_sca_filt)$leiden_sub <- clust
pos <- factor(colData(fib_signal_sca_filt)$position)
pos <- relevel(pos,"Lower")
colData(fib_signal_sca_filt)$position <- pos

levels(colData(fib_ecm_structural_sca_filt)$sample) <- c("LW_FIB_PWD18","LW_PWD12", "LW_PWD14","LW_REG_PWD18", "SW_PWD4", "SW_PWD7", "UW_P49", "UW_P21")
cond <- factor(colData(fib_ecm_structural_sca_filt)$sample)
cond <- relevel(cond,"UW_P49")
colData(fib_ecm_structural_sca_filt)$sample <- cond
clust <- factor(colData(fib_ecm_structural_sca_filt)$leiden_sub)
clust <- relevel(clust, "FIB-III")
colData(fib_ecm_structural_sca_filt)$leiden_sub <- clust
pos <- factor(colData(fib_ecm_structural_sca_filt)$position)
pos <- relevel(pos,"Lower")
colData(fib_ecm_structural_sca_filt)$position <- pos


# Build the Hurdle models, from previous analyses, you don't really get batch effects with mouse data, so we won't include sub_sample as a random effect
zlmModelSignalPosition <- zlm(formula = ~ position + n_genes, sca = fib_signal_sca_filt)
zlmModelEcmStructurePosition <- zlm(formula = ~ position + n_genes, sca = fib_ecm_structural_sca_filt)

zlmModelSignalCluster <- zlm(formula = ~ leiden_sub + n_genes, sca = fib_signal_sca_filt)
zlmModelEcmStructureCluster <- zlm(formula = ~ leiden_sub + n_genes, sca = fib_ecm_structural_sca_filt)

zlmModelSignal <- zlm(formula = ~ sample + n_genes, sca = fib_signal_sca_filt)
zlmModelEcmStructure <- zlm(formula = ~ sample + n_genes, sca = fib_ecm_structural_sca_filt)

### Test for DEGS between spatial position
contrastPositions <- c("positionUpper")
lrTestSignalPositions <- summary(zlmModelSignalPosition, doLRT=contrastPositions)
lrTestEcmStructurePositions <- summary(zlmModelEcmStructurePosition, doLRT=contrastPositions)

### Test for DEGs across clusters
contrastClusters <- c("leiden_subFIB-VIII", "leiden_subFIB-IX", "leiden_subFIB-IV", "leiden_subFIB-X", "leiden_subFIB-V", "leiden_subFIB-I", "leiden_subFIB-II", "leiden_subFIB-VI")
lrTestSignalClusters <- summary(zlmModelSignalCluster, doLRT=contrastClusters)
lrTestEcmStructureClusters <- summary(zlmModelEcmStructureCluster, doLRT=contrastClusters)

# Get the differentially expressed signalling genes per contrast
lrTestSignalSummaryDt <- lrTestSignalPositions$datatable

lrTestSignalResults <- merge(lrTestSignalSummaryDt[contrast=='positionUpper' & component=='H',.(primerid, `Pr(>Chisq)`)], #P-vals
                             lrTestSignalSummaryDt[contrast=='positionUpper' & component=='logFC', .(primerid, coef)],
                             by='primerid') #logFC coefficients

lrTestSignalHurdle <- merge(lrTestSignalSummaryDt[contrast=='positionUpper' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                            lrTestSignalSummaryDt[contrast=='positionUpper' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
FCTHRESHOLD <- log2(1.5)

lrTestSignalHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
lrTestSignalHurdleSig <- merge(lrTestSignalHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(fib_signal_sca_filt)), by='primerid')
setorder(lrTestSignalHurdleSig, fdr)

#Correct for multiple testing (FDR correction) and filtering
lrTestSignalResults[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')]
lrTestSignalDe <- lrTestSignalResults[lrTestSignalResults$FDR<0.01,, drop=F]
lrTestSignalDe <- lrTestSignalDe[order(lrTestSignalDe$FDR),]

lrTestSignalDeVolcano <-lrTestSignalDe[!is.na(lrTestSignalDe$coef)]
lrTestSignalVolcanoPlot <- EnhancedVolcano(lrTestSignalDeVolcano,
                                           lab = lrTestSignalDeVolcano$primerid,
                                           x = 'coef',
                                           y = 'Pr(>Chisq)', 
                                           title = 'Lower vs Upper',
                                           FCcutoff = 0.25,
                                           pCutoff = 10e-100,
                                           legendPosition = 'right',
                                           legendLabSize = 12,
                                           legendIconSize = 4.0,
                                           drawConnectors = TRUE,
                                           widthConnectors = 0.75)

lrTestSignalVolcanoPlot

lrTestSignalSummaryDt <- lrTestSignalPositions$datatable

lrTestSignalResults <- merge(lrTestSignalSummaryDt[contrast=='positionUpper' & component=='H',.(primerid, `Pr(>Chisq)`)], #P-vals
                             lrTestSignalSummaryDt[contrast=='positionUpper' & component=='logFC', .(primerid, coef)],
                             by='primerid') #logFC coefficients

lrTestSignalHurdle <- merge(lrTestSignalSummaryDt[contrast=='positionUpper' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                            lrTestSignalSummaryDt[contrast=='positionUpper' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
FCTHRESHOLD <- log2(1.5)

lrTestSignalHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
lrTestSignalHurdleSig <- merge(lrTestSignalHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(fib_signal_sca_filt)), by='primerid')
setorder(lrTestSignalHurdleSig, fdr)

#Correct for multiple testing (FDR correction) and filtering
lrTestSignalResults[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')]
lrTestSignalDe <- lrTestSignalResults[lrTestSignalResults$FDR<0.01,, drop=F]
lrTestSignalDe <- lrTestSignalDe[order(lrTestSignalDe$FDR),]

lrTestSignalDeVolcano <-lrTestSignalDe[!is.na(lrTestSignalDe$coef)]
lrTestSignalVolcanoPlot <- EnhancedVolcano(lrTestSignalDeVolcano,
                                           lab = lrTestSignalDeVolcano$primerid,
                                           x = 'coef',
                                           y = 'Pr(>Chisq)', 
                                           title = 'Lower vs Upper',
                                           FCcutoff = 0.5,
                                           pCutoff = 10e-100,
                                           legendPosition = 'right',
                                           legendLabSize = 12,
                                           legendIconSize = 4.0,
                                           drawConnectors = TRUE,
                                           widthConnectors = 0.75)

lrTestSignalVolcanoPlot

# Get the differentially expressed ECM structure genes per contrast
lrTestEcmStructureSummaryDt <- lrTestEcmStructurePositions$datatable

lrTestEcmStructureResults <- merge(lrTestEcmStructureSummaryDt[contrast=='positionUpper' & component=='H',.(primerid, `Pr(>Chisq)`)], #P-vals
                                   lrTestEcmStructureSummaryDt[contrast=='positionUpper' & component=='logFC', .(primerid, coef)],
                                   by='primerid') #logFC coefficients

lrTestEcmStructureHurdle <- merge(lrTestEcmStructureSummaryDt[contrast=='positionUpper' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                                  lrTestEcmStructureSummaryDt[contrast=='positionUpper' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
FCTHRESHOLD <- log2(1.5)

lrTestEcmStructureHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
lrTestEcmStructureHurdleSig <- merge(lrTestEcmStructureHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(fib_ecm_structural_sca_filt)), by='primerid')
setorder(lrTestEcmStructureHurdleSig, fdr)

#Correct for multiple testing (FDR correction) and filtering
lrTestEcmStructureResults[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')]
lrTestEcmStructureDe <- lrTestEcmStructureResults[lrTestEcmStructureResults$FDR<0.01,, drop=F]
lrTestEcmStructureDe <- lrTestEcmStructureDe[order(lrTestEcmStructureDe$FDR),]

lrTestEcmStructureDeVolcano <-lrTestEcmStructureDe[!is.na(lrTestEcmStructureDe$coef)]
lrTestEcmStructureVolcanoPlot <- EnhancedVolcano(lrTestEcmStructureDeVolcano,
                                                 lab = lrTestEcmStructureDeVolcano$primerid,
                                                 x = 'coef',
                                                 y = 'Pr(>Chisq)', 
                                                 title = 'Lower vs Upper',
                                                 FCcutoff = 0.25,
                                                 pCutoff = 10e-100,
                                                 legendPosition = 'right',
                                                 legendLabSize = 12,
                                                 legendIconSize = 4.0,
                                                 drawConnectors = TRUE,
                                                 widthConnectors = 0.75)

lrTestEcmStructureVolcanoPlot

### Test for DEGs across clusters
contrastClusters <- c("leiden_subFIB-VIII", "leiden_subFIB-IX", "leiden_subFIB-IV", "leiden_subFIB-X", "leiden_subFIB-V", "leiden_subFIB-I", "leiden_subFIB-II", "leiden_subFIB-VI")
lrTestSignalClusters <- summary(zlmModelSignalCluster, doLRT=contrastClusters)
lrTestEcmStructureClusters <- summary(zlmModelEcmStructureCluster, doLRT=contrastClusters)

# Get the differentially expressed signalling genes per contrast
lrTestSignalSummaryDt <- lrTestSignalClusters$datatable

lrTestSignalResults <- merge(lrTestSignalSummaryDt[contrast=='leiden_subFIB-II' & component=='H',.(primerid, `Pr(>Chisq)`)], #P-vals
                             lrTestSignalSummaryDt[contrast=='leiden_subFIB-II' & component=='logFC', .(primerid, coef)],
                             by='primerid') #logFC coefficients

lrTestSignalHurdle <- merge(lrTestSignalSummaryDt[contrast=='leiden_subFIB-II' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                            lrTestSignalSummaryDt[contrast=='leiden_subFIB-II' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
FCTHRESHOLD <- log2(1.5)

lrTestSignalHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
lrTestSignalHurdleSig <- merge(lrTestSignalHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(fib_signal_sca_filt)), by='primerid')
setorder(lrTestSignalHurdleSig, fdr)

#Correct for multiple testing (FDR correction) and filtering
lrTestSignalResults[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')]
lrTestSignalDe <- lrTestSignalResults[lrTestSignalResults$FDR<0.01,, drop=F]
lrTestSignalDe <- lrTestSignalDe[order(lrTestSignalDe$FDR),]

lrTestSignalDeVolcano <-lrTestSignalDe[!is.na(lrTestSignalDe$coef)]
lrTestSignalVolcanoPlot <- EnhancedVolcano(lrTestSignalDeVolcano,
                                           lab = lrTestSignalDeVolcano$primerid,
                                           x = 'coef',
                                           y = 'Pr(>Chisq)', 
                                           title = 'FIB-III vs FIB-II',
                                           FCcutoff = 0.5,
                                           pCutoff = 10e-100,
                                           legendPosition = 'right',
                                           legendLabSize = 12,
                                           legendIconSize = 4.0,
                                           drawConnectors = TRUE,
                                           widthConnectors = 0.75)

lrTestSignalVolcanoPlot

# Get the differentially expressed ECM structure genes per contrast
lrTestEcmStructureSummaryDt <- lrTestEcmStructureClusters$datatable

lrTestEcmStructureResults <- merge(lrTestEcmStructureSummaryDt[contrast=='leiden_subFIB-II' & component=='H',.(primerid, `Pr(>Chisq)`)], #P-vals
                                   lrTestEcmStructureSummaryDt[contrast=='leiden_subFIB-II' & component=='logFC', .(primerid, coef)],
                                   by='primerid') #logFC coefficients

lrTestEcmStructureHurdle <- merge(lrTestEcmStructureSummaryDt[contrast=='leiden_subFIB-II' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                                  lrTestEcmStructureSummaryDt[contrast=='leiden_subFIB-II' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
FCTHRESHOLD <- log2(1.5)

lrTestEcmStructureHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
lrTestEcmStructureHurdleSig <- merge(lrTestEcmStructureHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(fib_ecm_structural_sca_filt)), by='primerid')
setorder(lrTestEcmStructureHurdleSig, fdr)

#Correct for multiple testing (FDR correction) and filtering
lrTestEcmStructureResults[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')]
lrTestEcmStructureDe <- lrTestEcmStructureResults[lrTestEcmStructureResults$FDR<0.01,, drop=F]
lrTestEcmStructureDe <- lrTestEcmStructureDe[order(lrTestEcmStructureDe$FDR),]

lrTestEcmStructureDeVolcano <-lrTestEcmStructureDe[!is.na(lrTestEcmStructureDe$coef)]
lrTestEcmStructureVolcanoPlot <- EnhancedVolcano(lrTestEcmStructureDeVolcano,
                                                 lab = lrTestEcmStructureDeVolcano$primerid,
                                                 x = 'coef',
                                                 y = 'Pr(>Chisq)', 
                                                 title = 'FIB-III vs FIB-II',
                                                 FCcutoff = 0.5,
                                                 pCutoff = 10e-100,
                                                 legendPosition = 'right',
                                                 legendLabSize = 12,
                                                 legendIconSize = 4.0,
                                                 drawConnectors = TRUE,
                                                 widthConnectors = 0.75, 
                                                 xlim=c(0, 3))

lrTestEcmStructureVolcanoPlot

########################
### Test to see any genes differentially expressed across time
contrasts <- c("sampleUW_P21", "sampleSW_PWD7", "sampleLW_PWD12", "sampleLW_FIB_PWD18", "sampleLW_REG_PWD18")

lrTestFunctionalitySamples <- summary(zlmModelFunctionality, doLRT=contrasts)
lrTestSignalSamples <- summary(zlmModelSignal, doLRT=contrasts)
lrTestEcmStructureSamples <- summary(zlmModelEcmStructure, doLRT=contrasts)
lrTestEcmModifyingSamples <- summary(zlmModelEcmModifying, doLRT=contrasts)

# Get the differentially expressed signalling genes per contrast
lrTestFunctionalitySummaryDt <-lrTestFunctionalitySamples$datatable

lrTestFunctionalityResults <- merge(lrTestFunctionalitySummaryDt[contrast=='sampleLW_PWD12' & component=='H',.(primerid, `Pr(>Chisq)`)], #P-vals
                                    lrTestFunctionalitySummaryDt[contrast=='sampleLW_PWD12' & component=='logFC', .(primerid, coef)],
                                    by='primerid') #logFC coefficients

lrTestFunctionalityHurdle <- merge(lrTestFunctionalitySummaryDt[contrast=='sampleLW_PWD12' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                                   lrTestFunctionalitySummaryDt[contrast=='sampleLW_PWD12' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
FCTHRESHOLD <- log2(1.5)

lrTestFunctionalityHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
lrTestFunctionalityHurdleSig <- merge(lrTestFunctionalityHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(fib_signal_sca_filt)), by='primerid')
setorder(lrTestSignalHurdleSig, fdr)

#Correct for multiple testing (FDR correction) and filtering
lrTestFunctionalityResults[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')]
lrTestFunctionalityDe <- lrTestFunctionalityResults[lrTestFunctionalityResults$FDR<0.01,, drop=F]
lrTestFunctionalityDe <- lrTestFunctionalityDe[order(lrTestFunctionalityDe$FDR),]

lrTestFunctionalityDeVolcano <-lrTestFunctionalityDe[!is.na(lrTestFunctionalityDe$coef)]
lrTestFunctionalityVolcanoPlot <- EnhancedVolcano(lrTestFunctionalityDeVolcano,
                                                  lab = lrTestFunctionalityDeVolcano$primerid,
                                                  x = 'coef',
                                                  y = 'Pr(>Chisq)', 
                                                  title = 'UW P49 vs LW PWD12',
                                                  FCcutoff = 0.5,
                                                  pCutoff = 10e-100)
# legendPosition = 'right',
# legendLabSize = 12,
# legendIconSize = 4.0,
# drawConnectors = TRUE,
# widthConnectors = 0.75)

lrTestFunctionalityVolcanoPlot

lrTestSignalSummaryDt <- lrTestSignalSamples$datatable

lrTestSignalResults <- merge(lrTestSignalSummaryDt[contrast=='sampleLW_PWD12' & component=='H',.(primerid, `Pr(>Chisq)`)], #P-vals
                             lrTestSignalSummaryDt[contrast=='sampleLW_PWD12' & component=='logFC', .(primerid, coef)],
                             by='primerid') #logFC coefficients

lrTestSignalHurdle <- merge(lrTestSignalSummaryDt[contrast=='sampleLW_PWD12' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                            lrTestSignalSummaryDt[contrast=='sampleLW_PWD12' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
FCTHRESHOLD <- log2(1.5)

lrTestSignalHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
lrTestSignalHurdleSig <- merge(lrTestSignalHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(fib_signal_sca_filt)), by='primerid')
setorder(lrTestSignalHurdleSig, fdr)

#Correct for multiple testing (FDR correction) and filtering
lrTestSignalResults[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')]
lrTestSignalDe <- lrTestSignalResults[lrTestSignalResults$FDR<0.01,, drop=F]
lrTestSignalDe <- lrTestSignalDe[order(lrTestSignalDe$FDR),]

lrTestSignalDeVolcano <-lrTestSignalDe[!is.na(lrTestSignalDe$coef)]
lrTestSignalVolcanoPlot <- EnhancedVolcano(lrTestSignalDeVolcano,
                                           lab = lrTestSignalDeVolcano$primerid,
                                           x = 'coef',
                                           y = 'Pr(>Chisq)', 
                                           title = 'UW P49 vs LW PWD12',
                                           FCcutoff = 0.5,
                                           pCutoff = 10e-100,
                                           legendPosition = 'right',
                                           legendLabSize = 12,
                                           legendIconSize = 4.0,
                                           drawConnectors = TRUE,
                                           widthConnectors = 0.75)

lrTestSignalVolcanoPlot

# Get the differentially expressed ECM structure genes per contrast
lrTestEcmStructureSummaryDt <- lrTestEcmStructureSamples$datatable

lrTestEcmStructureResults <- merge(lrTestEcmStructureSummaryDt[contrast=='sampleLW_PWD12' & component=='H',.(primerid, `Pr(>Chisq)`)], #P-vals
                                   lrTestEcmStructureSummaryDt[contrast=='sampleLW_PWD12' & component=='logFC', .(primerid, coef)],
                                   by='primerid') #logFC coefficients

lrTestEcmStructureHurdle <- merge(lrTestEcmStructureSummaryDt[contrast=='sampleLW_PWD12' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                                  lrTestEcmStructureSummaryDt[contrast=='sampleLW_PWD12' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
FCTHRESHOLD <- log2(1.5)

lrTestEcmStructureHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
lrTestEcmStructureHurdleSig <- merge(lrTestEcmStructureHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(fib_ecm_structural_sca_filt)), by='primerid')
setorder(lrTestEcmStructureHurdleSig, fdr)

#Correct for multiple testing (FDR correction) and filtering
lrTestEcmStructureResults[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')]
lrTestEcmStructureDe <- lrTestEcmStructureResults[lrTestEcmStructureResults$FDR<0.01,, drop=F]
lrTestEcmStructureDe <- lrTestEcmStructureDe[order(lrTestEcmStructureDe$FDR),]

lrTestEcmStructureDeVolcano <-lrTestEcmStructureDe[!is.na(lrTestEcmStructureDe$coef)]
lrTestEcmStructureVolcanoPlot <- EnhancedVolcano(lrTestEcmStructureDeVolcano,
                                                 lab = lrTestEcmStructureDeVolcano$primerid,
                                                 x = 'coef',
                                                 y = 'Pr(>Chisq)', 
                                                 title = 'UW P49 vs LW PWD12',
                                                 FCcutoff = 0.5,
                                                 pCutoff = 10e-100,
                                                 legendPosition = 'right',
                                                 legendLabSize = 12,
                                                 legendIconSize = 4.0,
                                                 drawConnectors = TRUE,
                                                 widthConnectors = 0.75)

lrTestEcmStructureVolcanoPlot

# Get the differentially expressed ECM modifying genes per contrast
lrTestEcmModifyingSummaryDt <- lrTestEcmModifyingSamples$datatable

lrTestEcmModifyingResults <- merge(lrTestEcmModifyingSummaryDt[contrast=='sampleLW_REG_PWD18' & component=='H',.(primerid, `Pr(>Chisq)`)], #P-vals
                                   lrTestEcmModifyingSummaryDt[contrast=='sampleLW_REG_PWD18' & component=='logFC', .(primerid, coef)],
                                   by='primerid') #logFC coefficients

lrTestEcmModifyingHurdle <- merge(lrTestEcmModifyingSummaryDt[contrast=='sampleLW_REG_PWD18' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                                  lrTestEcmModifyingSummaryDt[contrast=='sampleLW_REG_PWD18' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
FCTHRESHOLD <- log2(1.5)

lrTestEcmModifyingHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
lrTestEcmModifyingHurdleSig <- merge(lrTestEcmModifyingHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(fib_ecm_modifying_sca_filt)), by='primerid')
setorder(lrTestEcmModifyingHurdleSig, fdr)

#Correct for multiple testing (FDR correction) and filtering
lrTestEcmModifyingResults[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')]
lrTestEcmModifyingDe <- lrTestEcmModifyingResults[lrTestEcmModifyingResults$FDR<0.01,, drop=F]
lrTestEcmModifyingDe <- lrTestEcmModifyingDe[order(lrTestEcmModifyingDe$FDR),]

lrTestEcmModifyingDeVolcano <-lrTestEcmModifyingDe[!is.na(lrTestEcmModifyingDe$coef)]
lrTestEcmModifyingVolcanoPlot <- EnhancedVolcano(lrTestEcmModifyingDeVolcano,
                                                 lab = lrTestEcmModifyingDeVolcano$primerid,
                                                 x = 'coef',
                                                 y = 'Pr(>Chisq)', 
                                                 title = 'UW P49 vs LW REG PWD18',
                                                 FCcutoff = 0.5,
                                                 pCutoff = 10e-100)
# legendPosition = 'right',
# legendLabSize = 12,
# legendIconSize = 4.0,
# drawConnectors = TRUE,
# widthConnectors = 0.75)

lrTestEcmModifyingVolcanoPlot

