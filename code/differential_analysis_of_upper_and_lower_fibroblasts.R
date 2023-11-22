library(mgcv)
library(SingleCellExperiment)
library(zellkonverter)
library(tradeSeq)
library(dplyr)
library(pheatmap)
library(clusterExperiment)
library(ComplexHeatmap)
library(circlize)


fibroblasts.sce <- readH5AD("/Users/axelalmet/Documents/CellCellCommunicationModelling/Data/IntegratedData/integratedfibroblastsdata.h5ad")

# Only consider hvgs
fibroblasts.sce <- fibroblasts.sce[rowData(fibroblasts.sce)$highly_variable == TRUE]

# Get the upper and lower scores as pseudotimes
upper.scores <- fibroblasts.sce$upper
lower.scores <- fibroblasts.sce$lower

# Fit for upper scores
fibroblast.counts <- counts(fibroblasts.sce)

pseudotime <- data.frame(curve1=upper.scores, curve2=lower.scores)
rownames(pseudotime) <- colnames(fibroblast.counts)

# In this case, I think all fibroblasts contribute to the upper/lower scores
cellweights <- data.frame(curve1=array(1, dim=dim(pseudotime)[1]), curve1=array(1, dim=dim(pseudotime)[1]))
rownames(cellweights) <- colnames(fibroblast.counts)

sce <- fitGAM(counts = fibroblast.counts,
              pseudotime = pseudotime,
              cellWeights = cellweights,
              nknots = 6, verbose = FALSE)

assocRes <- associationTest(sce)
assocRes$pvalue_adj <- p.adjust(assocRes$pvalue, method="bonferroni")

assocResSig <- assocRes[assocRes$pvalue_adj < 0.05, ]
assocResSig <- assocResSig[order(assocResSig$pvalue_adj), ]
assocResSig <- assocResSig[order(assocResSig$pvalue_adj, -assocResSig$meanLogFC), ]

assocGenes <- rownames(assocResSig)

tradeSeqOutput <- rowData(sce)$tradeSeq[rowData(sce)$tradeSeq$name %in% assocGenes, ]

# Cluster based on patterns (we really just want the y-hat scaled)
nPointsClus <- 20
upperGenes <- c('Crabp1', 'Dpp4', 'Fabp5', 'Lef1', 'Prdm1', 'Prss35', 'Runx1')
lowerGenes <- c('Cnn1', 'Dlk1', 'Fmo1', 'Ly6a', 'Mest')

clusPat <- clusterExpressionPatterns(sce, nPoints = nPointsClus,
                                     genes = assocGenes)

yHatUpper <- clusPat$yhatScaled[rownames(clusPat$yhatScaled) %in% assocGenes, 1:nPointsClus]
yHatLower <- clusPat$yhatScaled[rownames(clusPat$yhatScaled) %in% assocGenes, seq((nPointsClus + 1), 2*nPointsClus)]

pseudotime.vals.upper <- seq(min(pseudotime$curve1), max(pseudotime$curve1), length.out=20)
pseudotime.vals.lower <- seq(min(pseudotime$curve2), max(pseudotime$curve2), length.out=20)

upper.correlations <- array(0, dim=length(assocGenes))
lower.correlations <- array(0, dim=length(assocGenes))

for (i in seq(1, length(assocGenes)))
{
  upper.correlations[i] <- cor.test(yHatUpper[i, ], pseudotime.vals.upper, method="spearm")$estimate
  lower.correlations[i] <- cor.test(yHatLower[i, ], pseudotime.vals.lower, method="spearm")$estimate
}

yHatLowerDf <- data.frame(yHatLower)
yHatUpperDf <- data.frame(yHatUpper)

yHatLowerDf <- yHatLowerDf[order(lower.correlations), ]
yHatUpperDf <- yHatUpperDf[order(upper.correlations), ]

# Sort it manually because this is stupid
col_fun = colorRamp2(c(-4, 0, 4), hcl_palette = "Spectral", reverse=TRUE)
p_lower <- Heatmap(yHatLower[order(lower.correlations), ], col = col_fun, cluster_rows = FALSE, cluster_columns=FALSE, show_row_names=FALSE, show_column_names=FALSE)

p_upper <- Heatmap(yHatUpper[order(upper.correlations), ], col = col_fun, cluster_rows = FALSE, cluster_columns=FALSE, show_row_names=FALSE, show_column_names=FALSE)

write.csv(yHatLowerDf, "/Users/axelalmet/Documents/CellCellCommunicationModelling/Data/IntegratedData/integratedfibroblasts_lower_score_gam_fits.csv")
write.csv(yHatUpperDf, "/Users/axelalmet/Documents/CellCellCommunicationModelling/Data/IntegratedData/integratedfibroblasts_upper_score_gam_fits.csv")
