getwd()
setwd("/Users/priyarai/Documents/Researchproject/data/hydrogenperoxide")

#--------------------------------------------------------------------
#For WGCNA

library(WGCNA)
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)

#i've got the GSE42035 FPKM tracking data from Phytozome

gse34826_fpkm_data <- read_tsv("GSE34826_genes.fpkm_tracking")

head(gse34826_fpkm_data)
dim(gse34826_fpkm_data)

#in gse34826_fpkm_data I have 17,741 genes

#in gse34826_fpkm_data I want columns: 
# 4 (gene_id)
# 10 SRR394058_FPKM Hydrogen peroxide, t=0h after treatment (1)
# 14 SRR394059_FPKM Hydrogen peroxide, t=0h after treatment (2)
# 18 SRR394060_FPKM Hydrogen peroxide, t=0.5h after treatment (1)
# 22 SRR394061_FPKM Hydrogen peroxide, t=0.5h after treatment (2)
# 26 SRR394062_FPKM Hydrogen peroxide, t=1h after treatment (1)
# 30 SRR394063_FPKM Hydrogen peroxide, t=1h after treatment (2)

analysis_gse34826_fpkm_data <- gse34826_fpkm_data %>%
  select(4, 10, 14, 18, 22, 26, 30)

head(analysis_gse34826_fpkm_data)
dim(analysis_gse34826_fpkm_data)
#17741 7

write.csv(analysis_gse34826_fpkm_data, "analysis_gse34826_fpkm_data.csv", row.names = FALSE)

#---------------------------------------------------------------------
# Log-transform the FPKM values (adding a small constant to avoid log(0))
analysis_gse34826_log_fpkm_data <- analysis_gse34826_fpkm_data %>%
  mutate(across(starts_with("SRR"), ~ log2(. + 1)))

head(analysis_gse34826_log_fpkm_data)
dim(analysis_gse34826_log_fpkm_data)
#17741 7
#---------------------------------------------------------------------

#Just a table of genes of interest

genes_of_interest <- c("Cre07.g317250", "Cre06.g270500", "Cre06.g273100")

filtered_genes <- analysis_gse34826_log_fpkm_data %>%
  filter(gene_id %in% genes_of_interest)

head(filtered_genes)
dim(filtered_genes)

#---------------------------------------------------------------------
#now trying to make WGCNA

library(WGCNA)
library(dplyr)
library(readr)
library(tibble)

#---------------------------------------------------------------------
#when running the pickSoftThreshold the 0s are causing the NA so
#need to remove them

# Display column names to verify
print(names(analysis_gse34826_log_fpkm_data))

# Define the columns of interest based on the actual column names
columns_of_interest <- c("SRR394058_FPKM", "SRR394059_FPKM", "SRR394060_FPKM", "SRR394061_FPKM", "SRR394062_FPKM", "SRR394063_FPKM")

# Verify that the columns exist in the dataframe
missing_columns <- setdiff(columns_of_interest, names(analysis_gse34826_log_fpkm_data))
if (length(missing_columns) > 0) {
  stop("The following columns are missing in the dataframe: ", paste(missing_columns, collapse = ", "))
}

# Convert the specified columns to numeric
analysis_gse34826_log_fpkm_data[columns_of_interest] <- lapply(analysis_gse34826_log_fpkm_data[columns_of_interest], as.numeric)

# Create a logical condition to identify rows where all specified columns have a value of 0
rows_with_all_zeros <- apply(analysis_gse34826_log_fpkm_data[columns_of_interest], 1, function(row) all(row == 0))

# Print the number of rows identified
cat("Number of rows with all zeros in the specified columns:", sum(rows_with_all_zeros), "\n")
# Number of rows with all zeros in the specified columns: 76

# Filter out the rows with all zeros
analysis_gse34826_log_fpkm_data <- analysis_gse34826_log_fpkm_data[!rows_with_all_zeros, ]

cat("Dimensions after filtering:", dim(analysis_gse34826_log_fpkm_data), "\n")
#Dimensions after filtering: 17665 6

# Ensure gene_id is set as row names
analysis_gse34826_log_fpkm_data <- analysis_gse34826_log_fpkm_data %>% column_to_rownames("gene_id")

# Check the structure of the data
head(analysis_gse34826_log_fpkm_data)

# Transpose the data for WGCNA (genes as columns and samples as rows)
datExpr <- as.data.frame(t(analysis_gse34826_log_fpkm_data))

# Check the transposed data structure
str(datExpr)
dim(datExpr)

#6 17665

powers <- c(1:10, seq(from = 12, to = 50, by = 2))

sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
#Takes about a minute

write.csv(analysis_gse34826_log_fpkm_data, "analysis_gse34826_log_fpkm_data.csv", row.names = TRUE)

#Next step is to visualise the results from the pickSoftThreshold function

# Plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1, 2))

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type="n", 
     main = "Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=0.9, col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], 
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n", 
     main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9, col="red")

#--------------------------------------------------------------------

#Going with power 18 with SFT.R.sq as 0.896

# Set the soft threshold power
softPower <- 18

# Construct the adjacency matrix
adjacency <- adjacency(datExpr, power = softPower)
#This line takes about a minute
dim(adjacency)
#17665 17665
#calculates it for 17665 genes x 17665 genes hence why there are
#312,052,225 elements

# Turn adjacency into topological overlap matrix (TOM)
TOM <- TOMsimilarity(adjacency)
dim(TOM)
#17665 17665
#Takes about a minute

dissTOM <- 1 - TOM
#Takes about a minute
#in this matrix values of 0 mean very similar whereas 1 not similar at all

# Hierarchical clustering of the genes
geneTree <- hclust(as.dist(dissTOM), method = "average")
#Takes about a minute
plot(geneTree, main = "Gene clustering on TOM-based dissimilarity", sub = "", xlab = "")

# Dynamic tree cut to identify modules
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 0, pamRespectsDendro = FALSE)
#cutreeDynamic will calculate a cutHeight itself to give a reasonable
#number of modules
#set cutHeight to 0.991 ===> 99% of the (truncated) height range in dendro.

table(dynamicMods)
#There are 265 modules?

#assigns colours to the modules
dynamicColors <- labels2colors(dynamicMods)

# Plot the dendrogram and the module colors
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

#---------------------------------------------------------------------
#Other code to try and reduce number of modules

# Adjusted parameters for reducing the number of modules
deepSplitValue <- 0 # Lower sensitivity to module detection
cutHeightValue <- 0.991 # keep it to the one CutTreeDynamic chose
#because it gives a lower number of modules than when i used values closer to 1

# Dynamic tree cut to identify modules
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, 
                             deepSplit = deepSplitValue, 
                             pamRespectsDendro = FALSE,
                             cutHeight = cutHeightValue)

# Display the resulting module sizes
moduleSizes <- table(dynamicMods)
print(moduleSizes)

dynamicColors <- labels2colors(dynamicMods)

# Plotting the dendrogram with the dynamic tree cut
plotDendroAndColors(geneTree, dynamicMods, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


#---------------------------------------------------------------------
#After module identification you need to summarise the modules by
#calculating the module eigengenes

# Calculate module eigengenes
MEs <- moduleEigengenes(datExpr, colors = dynamicColors)$eigengenes

# Calculate the dissimilarity of module eigengenes
MEDiss <- 1 - cor(MEs)

# Cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average")

# Plot the module eigengene dendrogram
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
#---------------------------------------------------------------------

#Code to find module colour of genes of interest

# Assuming you have already run the WGCNA steps and have the dynamicMods and dynamicColors

# List of genes of interest
genes_of_interest <- c("Cre07.g317250", "Cre06.g270500", "Cre06.g273100")

# Ensure gene_id is set as row names
analysis_gse34826_log_fpkm_data <- analysis_gse34826_log_fpkm_data %>% column_to_rownames("gene_id")

# Check the structure of the data
head(analysis_gse34826_log_fpkm_data)

# Transpose the data for WGCNA (genes as columns and samples as rows)
datExpr <- as.data.frame(t(analysis_gse34826_log_fpkm_data))
#---------------------------------------------------------------------
#EROOR HERE - CODE RAN UP TO HERE
#---------------------------------------------------------------------
# Ensure that datExpr and dynamicMods have the same row names (gene IDs)
stopifnot(all(rownames(datExpr) == names(dynamicMods)))

# Create a data frame with gene IDs and their corresponding module colors
gene_module_assignment <- data.frame(
  gene_id = names(dynamicMods),
  module = dynamicColors
)

# Filter for the genes of interest
genes_of_interest_modules <- gene_module_assignment %>%
  filter(gene_id %in% genes_of_interest)

# Display the genes of interest with their assigned modules
print(genes_of_interest_modules)

#Next step is to relate modules to external/phenotypic traits - i think
#this is an optional step so unsure at this point to apply this to my
#network

#---------------------------------------------------------------------
#To export the network to Cytoscape
#Exporting all 5 modules to Cytoscape - there is also an option of 
#exporting just a single module

# Select all unique modules
modules <- unique(dynamicColors)
probes <- colnames(datExpr)

# Loop over each module to export the corresponding network data
#Takes about a minute
for (module in modules) {
  inModule <- (dynamicColors == module)
  modProbes <- probes[inModule]
  modTOM <- TOM[inModule, inModule]
  
  # Export the edge list for the module
  exportNetworkToCytoscape(
    modTOM,
    edgeFile = paste0("all-anaerobiosis-genes-edges-", module, ".txt"),
    nodeFile = paste0("all-anaerobiosis-genes-nodes-", module, ".txt"),
    weighted = TRUE,
    threshold = 0.02, # Adjust this threshold based on your preference
    nodeNames = modProbes,
    nodeAttr = dynamicColors[inModule]
  )
}

# Cre07.g317250 + Cre06.g270500 are in BLUE module
# Cre06.g273100 are in PURPLE module

#---------------------------------------------------------------------
#Here extracting the edges relating to Cre07.g317250 in the blue module

anaerobiosis_gse42035_blue_module <- read.table("all-anaerobiosis-genes-edges/all-anaerobiosis-genes-edges-blue.txt", header = TRUE, sep = "\t")

# Filter edges containing 'Cre07.g317250' in either 'fromNode' or 'toNode'
cre07.g317250_anaerobiosis_gse42035_blue_module <- anaerobiosis_gse42035_blue_module %>%
  filter(fromNode == "Cre07.g317250" | toNode == "Cre07.g317250")

head(cre07.g317250_anaerobiosis_gse42035_blue_module)
dim(cre07.g317250_anaerobiosis_gse42035_blue_module)
#2332 6

#Order by weight column in descending order
cre07.g317250_anaerobiosis_gse42035_blue_module <- cre07.g317250_anaerobiosis_gse42035_blue_module[order(-cre07.g317250_anaerobiosis_gse42035_blue_module$weight), ]

write.table(cre07.g317250_anaerobiosis_gse42035_blue_module, "cre07.g317250_anaerobiosis_gse42035_blue_module.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#---------------------------------------------------------------------
#Here extracting the edges relating to Cre06.g270500 in the blue module

anaerobiosis_gse42035_blue_module <- read.table("all-anaerobiosis-genes-edges/all-anaerobiosis-genes-edges-blue.txt", header = TRUE, sep = "\t")

# Filter edges containing 'Cre06.g270500' in either 'fromNode' or 'toNode'
cre06.g270500_anaerobiosis_gse42035_blue_module <- anaerobiosis_gse42035_blue_module %>%
  filter(fromNode == "Cre06.g270500" | toNode == "Cre06.g270500")

head(cre06.g270500_anaerobiosis_gse42035_blue_module)
dim(cre06.g270500_anaerobiosis_gse42035_blue_module)
#2332 6
#dim values are the same as Cre07.g317250 because they belong to the same module
#and there are 2332 + 1(gene of interest) nodes in the blue module

cre06.g270500_anaerobiosis_gse42035_blue_module <- cre06.g270500_anaerobiosis_gse42035_blue_module[order(-cre06.g270500_anaerobiosis_gse42035_blue_module$weight), ]

write.table(cre06.g270500_anaerobiosis_gse42035_blue_module, "cre06.g270500_anaerobiosis_gse42035_blue_module.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#---------------------------------------------------------------------
#Here extracting the edges relating to Cre06.g273100 in the purple module

anaerobiosis_gse42035_purple_module <- read.table("all-anaerobiosis-genes-edges/all-anaerobiosis-genes-edges-purple.txt", header = TRUE, sep = "\t")

# Filter edges containing 'Cre06.g273100' in either 'fromNode' or 'toNode'
cre06.g273100_anaerobiosis_gse42035_purple_module <- anaerobiosis_gse42035_purple_module %>%
  filter(fromNode == "Cre06.g273100" | toNode == "Cre06.g273100")

head(cre06.g273100_anaerobiosis_gse42035_purple_module)
dim(cre06.g273100_anaerobiosis_gse42035_purple_module)
#1009 6

cre06.g273100_anaerobiosis_gse42035_purple_module <- cre06.g273100_anaerobiosis_gse42035_purple_module[order(-cre06.g273100_anaerobiosis_gse42035_purple_module$weight), ]

write.table(cre06.g273100_anaerobiosis_gse42035_purple_module, "cre06.g273100_anaerobiosis_gse42035_purple_module.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#---------------------------------------------------------------------
#Extract edges of the top 100 weight

# Extract the top edges until we have X no. of unique genes in the 'toNode' column
unique_genes_toNode <- unique(cre06.g273100_anaerobiosis_gse42035_purple_module$toNode)

# Check if we have at least 100 unique genes
if (length(unique_genes_toNode) > 100) {
  unique_genes_toNode <- unique_genes_toNode[1:100]
}

# Export the list of unique genes to a new file
write.table(unique_genes_toNode, "Top100Genes_toNode.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)


