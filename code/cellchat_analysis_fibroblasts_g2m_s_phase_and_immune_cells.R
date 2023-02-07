### Application of Suoqin Jin's CellChat package to various datas on unwounded and wounded murine skin tissue data.
### We're following the tutorial found here:https://github.com/sqjin/CellChat/blob/master/vignettes/walkthrough_wound.html

# Load the relevant packages
library(dplyr)
library(CellChat)
library(SummarizedExperiment)
library(zellkonverter)

### Set the ligand-receptor database. Here we will use the "Secreted signalling" database for cell-cell communication (let's look at ECM-receptor in the future )
CellChatDB <- CellChatDB.mouse # The othe roption is CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # Other options include ECM-Receptor and Cell-Cell Contact

# To facilitate cell-cell communication, we will set a minimum percentage of the total population
min_percentage <- 0.1 # This should remove most of the bad cells (we do this in Scanpy to fix visualisation issues)

integrated_sce <- readH5AD("../data/integratedskindata.h5ad")

integrated.data.input <- assay(integrated_sce, "X")

integrated.meta.data <- data.frame(colData(integrated_sce))

##################################################################################################################################################################
##################                    Analyse CCC between G2/M + S-phase fibroblast and immune cells                             #################################
##################################################################################################################################################################

### SW PWD4
swpwd4.meta.data <- integrated.meta.data[(integrated.meta.data$sample == "SW PWD4"), ]
swpwd4.meta.data$leiden_sub <- droplevels(swpwd4.meta.data$leiden_sub)
levels(swpwd4.meta.data$leiden_sub)

# Get the relevant immuen adn fibroblast clusters of interest
swpwd4_immune_clusters_of_interest <- c("DC-1", "DC-2", "DC-3", "MAC-1", "MAC-2", "MAC-3", "MAC-5", "BASO", "MAST", "NEU", "NK", "TCELL", "pDC")
swpwd4_fib_clusters_of_interest <- c("FIB-IV")

swpwd4.meta.data.subset <- swpwd4.meta.data[(swpwd4.meta.data$leiden_sub %in% swpwd4_immune_clusters_of_interest)|
                                              ((swpwd4.meta.data$leiden_sub %in% swpwd4_fib_clusters_of_interest)&(swpwd4.meta.data$phase %in% c("G2M", "S"))),] 

swpwd4.meta.data.subset$leiden_sub <- droplevels(swpwd4.meta.data.subset$leiden_sub)

swpwd4.use <- rownames(swpwd4.meta.data.subset)
swpwd4.data.input <- integrated.data.input[, swpwd4.use]

swpwd4.identity <- data.frame(group = swpwd4.meta.data.subset$leiden_sub, row.names = row.names(swpwd4.meta.data.subset))

# Create the cellchat object
swpwd4.cc <-createCellChat(object = swpwd4.data.input, do.sparse = T, meta = swpwd4.identity, group.by = "group")
levels(swpwd4.cc@idents) # show factor levels of the cell labels

swpwd4GroupSize <- as.numeric(table(swpwd4.cc@idents)) # Get the number of cells in each group
swpwd4MinCells <- floor(min_percentage / 100.0 * sum(swpwd4GroupSize))

swpwd4.cc@DB <- CellChatDB.use # Set the database for the unwounded data

# We now identify over-expressed ligands/receptors in a cell group and then project gene expression data onto the protein-protein interaction network
swpwd4.cc <- subsetData(swpwd4.cc) # We subset the expression data of signalling genes to save on computational cost
swpwd4.cc <- identifyOverExpressedGenes(swpwd4.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
swpwd4.cc <- identifyOverExpressedInteractions(swpwd4.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
swpwd4.cc <- projectData(swpwd4.cc, PPI.mouse) # Other option includes PPI.human We're told that we may have to comment these out.

# We now infer the cell-cell communication network by calculating the communication probabilities
swpwd4.cc <- computeCommunProb(swpwd4.cc, raw.use = FALSE, population.size = FALSE)
swpwd4.cc <- filterCommunication(swpwd4.cc, min.cells = swpwd4MinCells) # Filter out the clusters with small numbers
swpwd4.cc <- computeCommunProbPathway(swpwd4.cc) # Calculate the probabilities at the signalling level
swpwd4.cc <- aggregateNet(swpwd4.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

swpwd4.net <- subsetCommunication(swpwd4.cc)

### SW PWD7
swpwd7.meta.data <- integrated.meta.data[(integrated.meta.data$sample == "SW PWD7"), ]
swpwd7.meta.data$leiden_sub <- droplevels(swpwd7.meta.data$leiden_sub)
levels(swpwd7.meta.data$leiden_sub)

# Get the relevant immuen adn fibroblast clusters of interest
swpwd7_immune_clusters_of_interest <- c("DC-1", "DC-2", "DC-3", "MAC-1", "MAC-2", "MAC-3", "MAC-4", "BASO", "NK", "TCELL")
swpwd7_fib_clusters_of_interest <- c("FIB-II", "FIB-III", "FIB-IX")

swpwd7.meta.data.subset <- swpwd7.meta.data[(swpwd7.meta.data$leiden_sub %in% swpwd7_immune_clusters_of_interest)|
                                              ((swpwd7.meta.data$leiden_sub %in% swpwd7_fib_clusters_of_interest)&(swpwd7.meta.data$phase %in% c("G2M", "S"))),] 
swpwd7.meta.data.subset$leiden_sub <- droplevels(swpwd7.meta.data.subset$leiden_sub)

swpwd7.use <- rownames(swpwd7.meta.data.subset)
swpwd7.data.input <- integrated.data.input[, swpwd7.use]

swpwd7.identity <- data.frame(group = swpwd7.meta.data.subset$leiden_sub, row.names = row.names(swpwd7.meta.data.subset))

# Create the cellchat object
swpwd7.cc <-createCellChat(object = swpwd7.data.input, do.sparse = T, meta = swpwd7.identity, group.by = "group")
levels(swpwd7.cc@idents) # show factor levels of the cell labels

swpwd7GroupSize <- as.numeric(table(swpwd7.cc@idents)) # Get the number of cells in each group
swpwd7MinCells <- floor(min_percentage / 100.0 * sum(swpwd7GroupSize))

swpwd7.cc@DB <- CellChatDB.use # Set the database for the unwounded data

# We now identify over-expressed ligands/receptors in a cell group and then project gene expression data onto the protein-protein interaction network
swpwd7.cc <- subsetData(swpwd7.cc) # We subset the expression data of signalling genes to save on computational cost
swpwd7.cc <- identifyOverExpressedGenes(swpwd7.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
swpwd7.cc <- identifyOverExpressedInteractions(swpwd7.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
swpwd7.cc <- projectData(swpwd7.cc, PPI.mouse) # Other option includes PPI.human We're told that we may have to comment these out.

# We now infer the cell-cell communication network by calculating the communication probabilities
swpwd7.cc <- computeCommunProb(swpwd7.cc, raw.use = FALSE, population.size = FALSE)
swpwd7.cc <- filterCommunication(swpwd7.cc, min.cells = swpwd7MinCells) # Filter out the clusters with small numbers
swpwd7.cc <- computeCommunProbPathway(swpwd7.cc) # Calculate the probabilities at the signalling level
swpwd7.cc <- aggregateNet(swpwd7.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

swpwd7.net <- subsetCommunication(swpwd7.cc)

### LW PWD12
lwpwd12.meta.data <- integrated.meta.data[(integrated.meta.data$sample == "LW PWD12"), ]
lwpwd12.meta.data$leiden_sub <- droplevels(lwpwd12.meta.data$leiden_sub)
levels(lwpwd12.meta.data$leiden_sub)

# Get the relevant immuen adn fibroblast clusters of interest
lwpwd12_immune_clusters_of_interest <- c("DC-1", "DC-2", "DC-3", "MAC-1", "MAC-2","MAC-3", "MAC-4", "MAC-5", "BASO", "MAST", "NEU", "pDC", "TCELL")
lwpwd12_fib_clusters_of_interest <- c("FIB-I", "FIB-V")

lwpwd12.meta.data.subset <- lwpwd12.meta.data[(lwpwd12.meta.data$leiden_sub %in% lwpwd12_immune_clusters_of_interest)|
                                                ((lwpwd12.meta.data$leiden_sub %in% lwpwd12_fib_clusters_of_interest)&(lwpwd12.meta.data$phase %in% c("G2M", "S"))),] 
lwpwd12.meta.data.subset$leiden_sub <- droplevels(lwpwd12.meta.data.subset$leiden_sub)

lwpwd12.use <- rownames(lwpwd12.meta.data.subset)
lwpwd12.data.input <- integrated.data.input[, lwpwd12.use]
lwpwd12.identity <- data.frame(group = lwpwd12.meta.data.subset$leiden_sub, row.names = row.names(lwpwd12.meta.data.subset))

# Create the cellchat object
lwpwd12.cc <-createCellChat(object = lwpwd12.data.input, do.sparse = T, meta = lwpwd12.identity, group.by = "group")
levels(lwpwd12.cc@idents) # show factor levels of the cell labels

lwpwd12GroupSize <- as.numeric(table(lwpwd12.cc@idents)) # Get the number of cells in each group
lwpwd12MinCells <- floor(min_percentage / 100.0 * sum(lwpwd12GroupSize))

lwpwd12.cc@DB <- CellChatDB.use # Set the database for the unwounded data

# We now identify over-expressed ligands/receptors in a cell group and then project gene expression data onto the protein-protein interaction network
lwpwd12.cc <- subsetData(lwpwd12.cc) # We subset the expression data of signalling genes to save on computational cost
lwpwd12.cc <- identifyOverExpressedGenes(lwpwd12.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
lwpwd12.cc <- identifyOverExpressedInteractions(lwpwd12.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
lwpwd12.cc <- projectData(lwpwd12.cc, PPI.mouse) # Other option includes PPI.human We're told that we may have to comment these out.

# We now infer the cell-cell communication network by calculating the communication probabilities
lwpwd12.cc <- computeCommunProb(lwpwd12.cc, raw.use = FALSE, population.size = FALSE)
lwpwd12.cc <- filterCommunication(lwpwd12.cc, min.cells = lwpwd12MinCells) # Filter out the clusters with small numbers
lwpwd12.cc <- computeCommunProbPathway(lwpwd12.cc) # Calculate the probabilities at the signalling level
lwpwd12.cc <- aggregateNet(lwpwd12.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

lwpwd12.net <- subsetCommunication(lwpwd12.cc)

### LW PWD18
lwpwd18.meta.data <- integrated.meta.data[(integrated.meta.data$sample == "LW REG PWD18"), ]
lwpwd18.meta.data$leiden_sub <- droplevels(lwpwd18.meta.data$leiden_sub)
levels(lwpwd18.meta.data$leiden_sub)

# Get the relevant immuen adn fibroblast clusters of interest
lwpwd18_immune_clusters_of_interest <- c("DC-1", "DC-2", "DC-3", "MAC-1", "MAC-2","MAC-3", "BASO", "MAST", "NEU", "NK", "TCELL")
lwpwd18_fib_clusters_of_interest <- c("FIB-II", "FIB-VI")

lwpwd18.meta.data.subset <- lwpwd18.meta.data[(lwpwd18.meta.data$leiden_sub %in% lwpwd18_immune_clusters_of_interest)|
                                                ((lwpwd18.meta.data$leiden_sub %in% lwpwd18_fib_clusters_of_interest)&(lwpwd18.meta.data$phase %in% c("G2M", "S"))),] 
lwpwd18.meta.data.subset$leiden_sub <- droplevels(lwpwd18.meta.data.subset$leiden_sub)

lwpwd18.use <- rownames(lwpwd18.meta.data.subset)
lwpwd18.data.input <- integrated.data.input[, lwpwd18.use]
lwpwd18.identity <- data.frame(group = lwpwd18.meta.data.subset$leiden_sub, row.names = row.names(lwpwd18.meta.data.subset))

# Create the cellchat object
lwpwd18.cc <-createCellChat(object = lwpwd18.data.input, do.sparse = T, meta = lwpwd18.identity, group.by = "group")
levels(lwpwd18.cc@idents) # show factor levels of the cell labels

lwpwd18GroupSize <- as.numeric(table(lwpwd18.cc@idents)) # Get the number of cells in each group
lwpwd18MinCells <- floor(min_percentage / 100.0 * sum(lwpwd18GroupSize))

lwpwd18.cc@DB <- CellChatDB.use # Set the database for the unwounded data

# We now identify over-expressed ligands/receptors in a cell group and then project gene expression data onto the protein-protein interaction network
lwpwd18.cc <- subsetData(lwpwd18.cc) # We subset the expression data of signalling genes to save on computational cost
lwpwd18.cc <- identifyOverExpressedGenes(lwpwd18.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
lwpwd18.cc <- identifyOverExpressedInteractions(lwpwd18.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
lwpwd18.cc <- projectData(lwpwd18.cc, PPI.mouse) # Other option includes PPI.human We're told that we may have to comment these out.

# We now infer the cell-cell communication network by calculating the communication probabilities
lwpwd18.cc <- computeCommunProb(lwpwd18.cc, raw.use = FALSE, population.size = FALSE)
lwpwd18.cc <- filterCommunication(lwpwd18.cc, min.cells = lwpwd18MinCells) # Filter out the clusters with small numbers
lwpwd18.cc <- computeCommunProbPathway(lwpwd18.cc) # Calculate the probabilities at the signalling level
lwpwd18.cc <- aggregateNet(lwpwd18.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

lwpwd18.net <- subsetCommunication(lwpwd18.cc)

# Save the communication activities for later
write.csv(swpwd4.net, file = "../output/integratedfibroblasts_communications_dividing_swpwd4.csv")
write.csv(swpwd7.net, file = "../output/integratedfibroblasts_communications_dividing_swpwd7.csv")
write.csv(lwpwd12.net, file = "../output/integratedfibroblasts_communications_dividing_lwpwd12.csv")
write.csv(lwpwd18.net, file = "../output/integratedfibroblasts_communications_dividing_lwpwd18.csv")
