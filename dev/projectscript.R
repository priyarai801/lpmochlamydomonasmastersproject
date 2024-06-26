getwd()
setwd("/Users/priyarai/Documents/Researchproject/data")

#Wheat Thermotolerance performed Differential expression
#analysis using DESeq2 and the differential expression plots were visualised with volcano
#plots in the ggplot2 package also used ashr package to shrink the expression fold changes
#more info in paper) and had log2foldchanges and those with FDR-adjusted p-values were considered
#further for GO enrichment analysis

#will use merged_df_unique

#DESeq2 only works on raw counts but i have RPKM values - possible
#alternative is limma

#--------------------------------------------------------------------
#For WGCNA

library(WGCNA)
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)

#i've got the GSE42035 FPKM tracking data from Phytozome

gse42035_fpkm_data <- read_tsv("GSE42035_genes.fpkm_tracking")

head(gse42035_fpkm_data)
dim(gse42035_fpkm_data)

#in GSE42035_fpkm_data I have 17,741 genes
#I think next step is to figure out which columns I actually need
#the FPKMs are labelled with SRRs

#in gse42035_fpkm_data I want columns 4 (gene_id), 10 (SRR611223_FPKM CC-124 Light oxic)
#SRR611224 CC-124 Dark anoxic 0.5 hours column 50 and SRR611225 Dark anoxic 6 hours column 14


analysis_gse42035_fpkm_data <- gse42035_fpkm_data %>%
  select(4, 10, 50, 14)

head(analysis_gse42035_fpkm_data)
dim(analysis_gse42035_fpkm_data)
write.csv(analysis_gse42035_fpkm_data, "analysis_gse42035_fpkm_data.csv", row.names = FALSE)

#---------------------------------------------------------------------
# Log-transform the FPKM values (adding a small constant to avoid log(0))
analysis_gse42035_log_fpkm_data <- analysis_gse42035_fpkm_data %>%
  mutate(across(starts_with("SRR"), ~ log2(. + 1)))

head(analysis_gse42035_log_fpkm_data)
dim(analysis_gse42035_log_fpkm_data)
#---------------------------------------------------------------------

#below code is just the genes of interest to see if it's suitable
#to have a FPKM cutoff value of 1

genes_of_interest <- c("Cre07.g317250", "Cre06.g270500", "Cre06.g273100")

filtered_genes <- analysis_gse42035_fpkm_data %>%
  filter(gene_id %in% genes_of_interest)

head(filtered_genes)
dim(filtered_genes)

#I will remove the genes with 0 in at least 2 of the columns in either column 2, 3 or 4

has_two_zeros <- function(row) {
  sum(row == 0) >= 2
}

zero_filtered_gse42035_fpkm_data <- analysis_gse42035_fpkm_data %>%
  filter(!apply(select(., 2, 3, 4), 1, has_two_zeros))

head(zero_filtered_gse42035_fpkm_data)
dim(zero_filtered_gse42035_fpkm_data)

#Doing so reduced number of genes from 17741 -> 17353
#Number is still way too big

#To reduce it further did quantile reduction
#bear in mind threshold can not go over 0.796424
#Actually now that've logged the FPKM values the threshold can not be
#higher than 0.8451279

# Step 1: Calculate the mean FPKM for each gene across all samples
zero_filtered_gse42035_fpkm_data$mean_fpkm <- rowMeans(zero_filtered_gse42035_fpkm_data[ , -1], na.rm = TRUE)

# Step 2: Calculate the quantiles of the mean log FPKM
fpkm_quantiles <- quantile(zero_filtered_gse42035_fpkm_data$mean_fpkm, probs = seq(0, 1, 0.01))

# View the quantiles to choose a threshold
print(fpkm_quantiles)

#Looks like threshold of 15% is okay = 0.6692659
#for log it's the same too 15% = 0.7228476
#aka the 15th percentile quantile threshold
# Step 3: Choose a quantile threshold
threshold <- fpkm_quantiles["15%"]

# Step 4: Filter the genes using the chosen threshold
filtered_data <- zero_filtered_gse42035_fpkm_data %>%
  filter(mean_fpkm > threshold)
#Genes reduced to 14750 - still way too big

# Step 1: Filter the genes where all FPKM values are above the threshold in all samples
filtered_data <- zero_filtered_gse42035_fpkm_data %>%
  filter(SRR611223_FPKM > threshold &
           SRR611224_FPKM > threshold &
           SRR611225_FPKM > threshold)

# Remove the mean_fpkm column as it's no longer needed
filtered_data <- filtered_data %>%
  select(-mean_fpkm)

# View the filtered dataframe
dim(filtered_data)

#Genes reduced to 13775 - still way too big
#idk how logging made a difference but now there are 13825 genes

#---------------------------------------------------------------------
#Trying to reduce the 13825 genes with log  in filtered_data even further by only
#selecting the top 1000 highest variance genes as these are most likely
#to be biologically meaningful

write.csv(filtered_data, "delete.csv", row.names = FALSE)

# Calculate the variance for each gene across the FPKM columns
filtered_data <- zero_filtered_gse42035_fpkm_data %>%
  rowwise() %>%
  mutate(variance = var(c_across(SRR611223_FPKM:SRR611225_FPKM), na.rm = TRUE)) %>%
  ungroup()

# Select the top 1000 genes with the highest variance
top_genes <- filtered_data %>%
  arrange(desc(variance)) %>%
  slice(1:1000) %>%
  select(gene_id, SRR611223_FPKM, SRR611224_FPKM, SRR611225_FPKM)
#It's important you have the select line otherwise code will not reduce
#the number of rows to 1000

# View the top genes
print(head(top_genes))
dim(top_genes)
#yep selecting the top 1000 genes removes the genes of interest 

#--------------------------------------------------------------------
#So adding back in the genes of interest

# Step 3: Ensure inclusion of genes of interest
genes_of_interest <- c("Cre07.g317250", "Cre06.g270500", "Cre06.g273100")
additional_genes <- filtered_data %>%
  filter(gene_id %in% genes_of_interest)

# Combine the genes of interest with the top 1000 variance genes
combined_genes <- bind_rows(top_genes, additional_genes) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  arrange(desc(variance)) %>%
  slice(1:1004) %>%
  select(gene_id, SRR611223_FPKM, SRR611224_FPKM, SRR611225_FPKM)
  
# View the combined genes
head(combined_genes)
dim(combined_genes)
#No idea why after logging there are only 1002 genes here, but all the
#genes of interest are still here


write.csv(combined_genes,"combined_genes.csv", row.names=FALSE)
#--------------------------------------------------------------------
#now trying to make WGCNA

#I realised i did not log them, should i have??? - i'm not sure at what
#I should have logged them

library(WGCNA)
library(dplyr)
library(readr)
library(tibble)

# Ensure gene_id is set as row names
combined_genes <- combined_genes %>% column_to_rownames("gene_id")

# Check the structure of the data
head(combined_genes)

# Transpose the data for WGCNA (genes as columns and samples as rows)
datExpr <- as.data.frame(t(combined_genes))

# Check the transposed data structure
str(datExpr)
dim(datExpr)
#---------------------------------------------------------------------


#Says I have to check for outliers? wouldn't me deliberately adding back
#the genes of interest despite them not being the top 1000 variance genes
#make them more prone to being an outlier?
# sampleTree <- hclust(dist(datExpr), method = "average")
# plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "")


powers <- c(1:10, seq(from = 12, to = 50, by = 2))

sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

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

#Go with power 5, the scale-free score is 0.025 which is very low but
#number of samples is too low

# Set the soft threshold power
softPower <- 5

# Construct the adjacency matrix
adjacency <- adjacency(datExpr, power = softPower)

dim(adjacency)
#1002 1002
#calculates it for 1002 genes x 1002 genes hence why there are
#1004004 elements

# Turn adjacency into topological overlap matrix (TOM)
TOM <- TOMsimilarity(adjacency)
dim(TOM)
#1002 1002
dissTOM <- 1 - TOM
#in this matrix values of 0 mean very similar whereas 1 not similar at all

# Hierarchical clustering of the genes
geneTree <- hclust(as.dist(dissTOM), method = "average")
plot(geneTree, main = "Gene clustering on TOM-based dissimilarity", sub = "", xlab = "")


# Dynamic tree cut to identify modules
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE)
#cutreeDynamic will calculate a cutHeight itself to give a reasonable
#number of modules
#set cutHeight to 0.767 ===> 99% of the (truncated) height range in dendro.

table(dynamicMods)
#shows that there are 5 modules with the no. of genes in the corresponding
#columns

#assigns colours to the modules
dynamicColors <- labels2colors(dynamicMods)

# Plot the dendrogram and the module colors
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

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
for (module in modules) {
  inModule <- (dynamicColors == module)
  modProbes <- probes[inModule]
  modTOM <- TOM[inModule, inModule]
  
  # Export the edge list for the module
  exportNetworkToCytoscape(
    modTOM,
    edgeFile = paste0("CytoscapeInput-edges-", module, ".txt"),
    nodeFile = paste0("CytoscapeInput-nodes-", module, ".txt"),
    weighted = TRUE,
    threshold = 0.02, # Adjust this threshold based on your preference
    nodeNames = modProbes,
    nodeAttr = dynamicColors[inModule]
  )
}

#--------------------------------------------------------------------
#Searched the genes of interest, here are the module colours the
#genes belong to

# Cre07.g317250 is in the green module
# Cre06.g270500 is in turquoise module
# Cre06.g273100 is in the yellow module

#---------------------------------------------------------------------
#For now extract the top 2000 edges from the modules containing the
#genes of interest: green, turquiose, yellow

#Extract top 2000 edges from green module Cre07.g317250
library(dplyr)

green_edges <- read.table("CytoscapeInput-edges-green.txt", header = TRUE, sep = "\t")

# Sort edges by weight from highest to lowest
green_edges_sorted <- green_edges %>%
  arrange(desc(weight))

#Select the top 2000 edges
green_top_2000_edges <- green_edges_sorted %>%
  slice(1:2000)

#Export the filtered edges to a new file
write.table(green_top_2000_edges, "CytoscapeInput-edges-green-top2000.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#--------------------------------------------------------------------
#Extract top 2000 edges from turquoise module Cre06.g270500
#THIS DOES NOT WORK - THE TURQUOISE MODULE IS VERY BIG
#AND THE WEIGHT SCORE OF THE EDGES CORRESPONDING TO Cre06.g270500
#ARE NOWHERE NEAR THE TOP

turquoise_edges <- read.table("CytoscapeInput-edges-turquoise.txt", header = TRUE, sep = "\t")

# Sort edges by weight from highest to lowest
turquoise_edges_sorted <- turquoise_edges %>%
  arrange(desc(weight))

#Select the top 2000 edges
turquoise_top_2000_edges <- turquoise_edges_sorted %>%
  slice(1:2000)

#Export the filtered edges to a new file
write.table(turquoise_top_2000_edges, "CytoscapeInput-edges-turquoise-top2000.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#---------------------------------------------------------------------
#Extract top 2000 edges from yellow module Cre06.g273100

yellow_edges <- read.table("CytoscapeInput-edges-yellow.txt", header = TRUE, sep = "\t")

# Sort edges by weight from highest to lowest
yellow_edges_sorted <- yellow_edges %>%
  arrange(desc(weight))

#Select the top 2000 edges
yellow_top_2000_edges <- yellow_edges_sorted %>%
  slice(1:2000)

#Export the filtered edges to a new file
write.table(yellow_top_2000_edges, "CytoscapeInput-edges-yellow-top2000.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#---------------------------------------------------------------------
#I am foregoing the top 1000 weight edges because they remove the
#data for the genes of interest

#Here extracting the edges relating to Cre07.g317250 in the green module

library(dplyr)

green_edges <- read.table("CytoscapeInput-edges-green.txt", header = TRUE, sep = "\t")

# Filter edges containing 'Cre07.g317250' in either 'fromNode' or 'toNode'
green_gene_of_interest_edges <- green_edges %>%
  filter(fromNode == "Cre07.g317250" | toNode == "Cre07.g317250")

head(green_gene_of_interest_edges)
dim(green_gene_of_interest_edges)
#129 6

write.table(green_gene_of_interest_edges, "anaerobiosis_gse42035_edges_Cre07.g317250.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#---------------------------------------------------------------------
#Here extracting the edges relating to Cre06.g270500 in the turquoise module

library(dplyr)

turquoise_edges <- read.table("CytoscapeInput-edges-turquoise.txt", header = TRUE, sep = "\t")

# Filter edges containing 'Cre06.g270500' in either 'fromNode' or 'toNode'
turquoise_gene_of_interest_edges <- turquoise_edges %>%
  filter(fromNode == "Cre06.g270500" | toNode == "Cre06.g270500")

head(turquoise_gene_of_interest_edges)
dim(turquoise_gene_of_interest_edges)
#313 6

write.table(turquoise_gene_of_interest_edges, "anaerobiosis_gse42035_edges_Cre06.g270500.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#---------------------------------------------------------------------
#Here extracting the edges relating to Cre06.g273100 in the yellow module

library(dplyr)

yellow_edges <- read.table("CytoscapeInput-edges-yellow.txt", header = TRUE, sep = "\t")

# Filter edges containing 'Cre06.g270500' in either 'fromNode' or 'toNode'
yellow_gene_of_interest_edges <- yellow_edges %>%
  filter(fromNode == "Cre06.g273100" | toNode == "Cre06.g273100")

head(yellow_gene_of_interest_edges)
dim(yellow_gene_of_interest_edges)
#143 6

write.table(yellow_gene_of_interest_edges, "anaerobiosis_gse42035_edges_Cre06.g273100.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#---------------------------------------------------------------------

#Get list of coexpressed genes wth Cre07.g317250 for GO enrichment analysis

# Load necessary libraries
library(dplyr)

# Read the edge file
edges <- read.table("anaerobiosis_gse42035_edges_Cre06.g273100.txt", header = TRUE, sep = "\t")

# Extract unique genes from the 'toNode' column
unique_genes_toNode <- unique(edges$toNode)

# Check the list of unique genes
head(unique_genes_toNode)
length(unique_genes_toNode)

# Export the list of unique genes to a new file
write.table(unique_genes_toNode, "UniqueGenes_toNode.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#Seeing if it is possible to create the coexpression edges list again
#but with no filtering

library(readr)
library(dplyr)

gse42035_fpkm_data <- read_tsv("GSE42035_genes.fpkm_tracking")

head(gse42035_fpkm_data)
dim(gse42035_fpkm_data)

#in GSE42035_fpkm_data I have 17,741 genes

#in gse42035_fpkm_data I want columns 4 (gene_id), 10 (SRR611223_FPKM CC-124 Light oxic)
#SRR611224 CC-124 Dark anoxic 0.5 hours column 50 and SRR611225 Dark anoxic 6 hours column 14

analysis_gse42035_fpkm_data <- gse42035_fpkm_data %>%
  select(4, 10, 50, 14)

head(analysis_gse42035_fpkm_data)
dim(analysis_gse42035_fpkm_data)

# Log-transform the FPKM values (adding a small constant to avoid log(0))
analysis_gse42035_log_fpkm_data <- analysis_gse42035_fpkm_data %>%
  mutate(across(starts_with("SRR"), ~ log2(. + 1)))

head(analysis_gse42035_log_fpkm_data)
dim(analysis_gse42035_log_fpkm_data)
#17741 4

write.csv(analysis_gse42035_log_fpkm_data, "analysis_gse42035_log_fpkm_data.csv", row.names = TRUE)

#filtered_genes just has the log_fpkm of the lpmo genes
genes_of_interest <- c("Cre07.g317250", "Cre06.g270500", "Cre06.g273100")
filtered_genes <- analysis_gse42035_log_fpkm_data %>%
  filter(gene_id %in% genes_of_interest)
head(filtered_genes)
dim(filtered_genes)

#Now i have to create a WGCNA with the 17741 genes in analysis_gse42035_log_fpkm_data
library(WGCNA)
library(dplyr)
library(readr)
library(tibble)

# Ensure gene_id is set as row names
no_zero_analysis_gse42035_log_fpkm_data <- no_zero_analysis_gse42035_log_fpkm_data %>% column_to_rownames("gene_id")

# Check the structure of the data
head(analysis_gse42035_log_fpkm_data)

# Transpose the data for WGCNA (genes as columns and samples as rows)
datExpr <- as.data.frame(t(analysis_gse42035_log_fpkm_data))

# Check the transposed data structure
str(datExpr)
dim(datExpr)
# 3 17741

powers <- c(1:10, seq(from = 12, to = 50, by = 2))

sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
#Warning message with this line "1: In eval(xpr, envir = envir) :
#Some correlations are NA in block 1 : 2521 .
#still does perform the pickSoftThreshold function but gives the
#warning messages above

#I'm pretty sure that the error here is that there are genes with
#different formats of 0 in all 3 columns and R is not able to tell this

# Display column names to verify
print(names(analysis_gse42035_log_fpkm_data))

# Define the columns of interest based on the actual column names
columns_of_interest <- c("SRR611223_FPKM", "SRR611224_FPKM", "SRR611225_FPKM")

# Verify that the columns exist in the dataframe
missing_columns <- setdiff(columns_of_interest, names(analysis_gse42035_log_fpkm_data))
if (length(missing_columns) > 0) {
  stop("The following columns are missing in the dataframe: ", paste(missing_columns, collapse = ", "))
}

# Convert the specified columns to numeric
analysis_gse42035_log_fpkm_data[columns_of_interest] <- lapply(analysis_gse42035_log_fpkm_data[columns_of_interest], as.numeric)

# Create a logical condition to identify rows where all specified columns have a value of 0
rows_with_all_zeros <- apply(analysis_gse42035_log_fpkm_data[columns_of_interest], 1, function(row) all(row == 0))

# Print the number of rows identified
cat("Number of rows with all zeros in the specified columns:", sum(rows_with_all_zeros), "\n")
# Number of rows with all zeros in the specified columns: 202

# Filter out the rows with all zeros
no_zero_analysis_gse42035_log_fpkm_data <- analysis_gse42035_log_fpkm_data[!rows_with_all_zeros, ]


# Print the dimensions of the data before and after filtering
cat("Dimensions before filtering:", dim(analysis_gse42035_log_fpkm_data), "\n")
#Dimensions before filtering: 17741 3
cat("Dimensions after filtering:", dim(no_zero_analysis_gse42035_log_fpkm_data), "\n")
#Dimensions after filtering: 17539 3

write.csv(no_zero_analysis_gse42035_log_fpkm_data, "no_zero_analysis_gse42035_log_fpkm_data.csv", row.names = TRUE)

#Run the pickSoftThreshold function again now that the genes causing the
#NA values should be rid of

# Ensure gene_id is set as row names
no_zero_analysis_gse42035_log_fpkm_data <- no_zero_analysis_gse42035_log_fpkm_data %>% column_to_rownames("gene_id")
#It already was hence why the error comes up

# Check the structure of the data
head(no_zero_analysis_gse42035_log_fpkm_data)

# Transpose the data for WGCNA (genes as columns and samples as rows)
datExpr <- as.data.frame(t(no_zero_analysis_gse42035_log_fpkm_data))

# Check the transposed data structure
str(datExpr)
dim(datExpr)
# 3 17539
powers <- c(1:10, seq(from = 12, to = 50, by = 2))

sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

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

#Going with power 12 with SFT.R.sq as 0.0545 , slope -1.08,
#truncated.R.sq 0.818 m mean.k. 4060

# Set the soft threshold power
softPower <- 12

# Construct the adjacency matrix
adjacency <- adjacency(datExpr, power = softPower)
#This line takes about a minute
dim(adjacency)
#17539 17539
#calculates it for 17539 genes x 17539 genes hence why there are
#307,616,521 elements

# Turn adjacency into topological overlap matrix (TOM)
TOM <- TOMsimilarity(adjacency)
dim(TOM)
#17539 17539
#Takes about a minute

dissTOM <- 1 - TOM
#Takes about a minute
#in this matrix values of 0 mean very similar whereas 1 not similar at all

# Hierarchical clustering of the genes
geneTree <- hclust(as.dist(dissTOM), method = "average")
#Takes about a minute
plot(geneTree, main = "Gene clustering on TOM-based dissimilarity", sub = "", xlab = "")

# Dynamic tree cut to identify modules
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE)
#cutreeDynamic will calculate a cutHeight itself to give a reasonable
#number of modules
#set cutHeight to 0.906 ===> 99% of the (truncated) height range in dendro.

table(dynamicMods)
#shows that there are 11 modules with the no. of genes in the corresponding
#columns

#assigns colours to the modules
dynamicColors <- labels2colors(dynamicMods)

# Plot the dendrogram and the module colors
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

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
unique_genes_toNode <- unique(cre07.g317250_anaerobiosis_gse42035_blue_module$toNode)

# Check if we have at least 100 unique genes
if (length(unique_genes_toNode) > 400) {
  unique_genes_toNode <- unique_genes_toNode[1:400]
}

# Export the list of unique genes to a new file
write.table(unique_genes_toNode, "Top400Genes_toNode.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)


