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

# Transpose the data for WGCNA (genes as columns and samples as rows)s
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
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE)
#cutreeDynamic will calculate a cutHeight itself to give a reasonable
#number of modules
#set cutHeight to 0.993 ===> 99% of the (truncated) height range in dendro.

table(dynamicMods)
#There are 263 modules?

#assigns colours to the modules
dynamicColors <- labels2colors(dynamicMods)

# Plot the dendrogram and the module colors
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

#---------------------------------------------------------------------
#Other code to try and reduce number of modules - not running this part
#i'm scared of modules being too big and non-specific

# # Adjusted parameters for reducing the number of modules
# deepSplitValue <- 0 # Lower sensitivity to module detection
# cutHeightValue <- 0.991 # keep it to the one CutTreeDynamic chose
# #because it gives a lower number of modules than when i used values closer to 1
# 
# # Dynamic tree cut to identify modules
# dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, 
#                              deepSplit = deepSplitValue, 
#                              pamRespectsDendro = FALSE,
#                              cutHeight = cutHeightValue)
# 
# # Display the resulting module sizes
# moduleSizes <- table(dynamicMods)
# print(moduleSizes)
# 
# #---------------------------------------------------------------------
# dynamicColors <- labels2colors(dynamicMods)
# 
# # Plotting the dendrogram with the dynamic tree cut
# plotDendroAndColors(geneTree, dynamicMods, "Dynamic Tree Cut",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)
# 
# 
#---------------------------------------------------------------------
#After module identification you need to summarise the modules by
#calculating the module eigengenes

#Continue running code here
dynamicColors <- labels2colors(dynamicMods)

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

#Next step is to relate modules to external/phenotypic traits - i think
#this is an optional step so unsure at this point to apply this to my
#network

#---------------------------------------------------------------------
#Trying to find out which module the LPMO genes are in

# Create a data frame with gene IDs and their corresponding module colors
gene_module_df <- data.frame(
  gene_id = rownames(analysis_gse34826_log_fpkm_data),
  module = dynamicColors
)

# Filter for genes of interest
genes_of_interest_modules <- gene_module_df %>%
  filter(gene_id %in% genes_of_interest)

# Print the results
print(genes_of_interest_modules)

#Cre07.g317250 lightslateblue module
#Cre06.g270500 lightcyan module
#Cre06.g273100 grey60 module

#---------------------------------------------------------------------

#Create dataframe of all the genes in the dataset that are designated
#to the module of interest

#Cre07.g317250 lightslateblue module
lightslateblue_module_genes <- gene_module_df %>%
  filter(module == "lightslateblue")

#Add in the log FPKM values too
lightslateblue_module_genes_with_fpkm <- lightslateblue_module_genes %>%
  left_join(analysis_gse34826_log_fpkm_data %>% rownames_to_column("gene_id"), by = "gene_id")

write.csv(lightslateblue_module_genes_with_fpkm, "lightslateblue_module_genes_with_fpkm.csv", row.names = FALSE)

#---------------------------------------------------------------------

#Do the same for the other 2 modules

#Cre06.g270500 lightcyan module
lightcyan_module_genes <- gene_module_df %>%
  filter(module == "lightcyan")

lightcyan_module_genes_with_fpkm <- lightcyan_module_genes %>%
  left_join(analysis_gse34826_log_fpkm_data %>% rownames_to_column("gene_id"), by = "gene_id")

write.csv(lightcyan_module_genes_with_fpkm, "lightcyan_module_genes_with_fpkm.csv", row.names = FALSE)

#---------------------------------------------------------------------

#Cre06.g273100 grey60 module
grey60_module_genes <- gene_module_df %>%
  filter(module == "grey60")

grey60_module_genes_with_fpkm <- grey60_module_genes %>%
  left_join(analysis_gse34826_log_fpkm_data %>% rownames_to_column("gene_id"), by = "gene_id")

write.csv(grey60_module_genes_with_fpkm, "grey60_module_genes_with_fpkm.csv", row.names = FALSE)

#---------------------------------------------------------------------
#For Cre07.g317250 black module can i also have a list of the edges and weight in wgcna

#---------------------------------------------------------------------
#Just to get a list of gene_ids for a quick GO enrichment analysis

#Cre07.g317250 lightslateblue module
lightslateblue_module_gene_ids <- gene_module_df %>%
  filter(module == "lightslateblue") %>%
  pull(gene_id)

lightslateblue_module_gene_ids <- data.frame(gene_id = lightslateblue_module_gene_ids)

write.csv(lightslateblue_module_gene_ids, "lightslateblue_module_genes_ids.csv", row.names = FALSE)

#---------------------------------------------------------------------
#seeing if I can do GO analysis here as gprofiler website has limit of
#100 genes

install.packages("gprofiler2")
library(gprofiler2)

#Get a list of just the genes for GO analysis
lightslateblue_module_gene_ids_list <- lightslateblue_module_gene_ids$gene_id

gost_res <- gost(query = lightslateblue_module_gene_ids_list,
                 organism = "creinhardtii",
                 ordered_query = FALSE,
                 multi_query = FALSE,
                 significant = TRUE,
                 exclude_iea = FALSE,
                 measure_underrepresentation = FALSE,
                 evcodes = TRUE,
                 user_threshold = 0.05,
                 correction_method = "g_SCS")

print(gost_res$result)
#---------------------------------------------------------------------
#Cre06.g270500 lightcyan module

lightcyan_module_gene_ids <- gene_module_df %>%
  filter(module == "lightcyan") %>%
  pull(gene_id)

lightcyan_module_gene_ids <- data.frame(gene_id = lightcyan_module_gene_ids)

write.csv(lightcyan_module_gene_ids, "lightcyan_module_genes_ids.csv", row.names = FALSE)

lightcyan_module_gene_ids_list <- lightcyan_module_gene_ids$gene_id

gost_res <- gost(query = lightcyan_module_gene_ids_list,
                 organism = "creinhardtii",
                 ordered_query = FALSE,
                 multi_query = FALSE,
                 significant = TRUE,
                 exclude_iea = FALSE,
                 measure_underrepresentation = FALSE,
                 evcodes = TRUE,
                 user_threshold = 0.05,
                 correction_method = "g_SCS")

print(gost_res$result)

#---------------------------------------------------------------------
#Cre06.g273100 grey60 module

grey60_module_gene_ids <- gene_module_df %>%
  filter(module == "grey60") %>%
  pull(gene_id)

grey60_module_gene_ids <- data.frame(gene_id = grey60_module_gene_ids)

write.csv(grey60_module_gene_ids, "grey60_module_genes_ids.csv", row.names = FALSE)

grey60_module_gene_ids_list <- grey60_module_gene_ids$gene_id

gost_res <- gost(query = grey60_module_gene_ids_list,
                 organism = "creinhardtii",
                 ordered_query = FALSE,
                 multi_query = FALSE,
                 significant = TRUE,
                 exclude_iea = FALSE,
                 measure_underrepresentation = FALSE,
                 evcodes = TRUE,
                 user_threshold = 0.05,
                 correction_method = "g_SCS")

print(gost_res$result)

#---------------------------------------------------------------------

#Trying to get a list of all the edges and modules and weights in 
# analysis_gse34826_log_fpkm_data so that they can be entered into
#cytoscape to create a weighted gene coexpression network

# 1. Extracting the adjacency matrix and module assignments
#this part's already been run previously - but wouldn't hurt to run
#again I suppose
adjacency <- adjacency(datExpr, power = softPower)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

# 2. Convert the adjacency matrix to an edge list
threshold <- 0.1  # You can adjust this threshold to filter edges by weight
adjacency[adjacency < threshold] <- 0
#takes a minute

# Create edge list from adjacency matrix
library(reshape2)
edge_list <- melt(adjacency)
colnames(edge_list) <- c("Source", "Target", "Weight")
edge_list <- edge_list[edge_list$Weight > 0, ]

# Save edge list to CSV - takes a minute
write.csv(edge_list, "edge_list.csv", row.names = FALSE)

head(edge_list)
str(edge_list)

#Network is way too unprecedently big to visualise in Cytoscape so now
#going to visualise each module with the LPMO gene in

#---------------------------------------------------------------------

#Getting list of edges and weights for visualisation in cytoscape
#for lightslateblue module for Cre07.g317250

#---------------------------------------------------------------------
# Ensure gene_id is set as row names
analysis_gse34826_log_fpkm_data <- analysis_gse34826_log_fpkm_data %>% 
  column_to_rownames("gene_id")

# Transpose the data for WGCNA (genes as columns and samples as rows)
datExpr <- as.data.frame(t(analysis_gse34826_log_fpkm_data))

# Check the structure of the transposed data
str(datExpr)
dim(datExpr)

# Construct adjacency matrix
softPower <- 18
adjacency <- adjacency(datExpr, power = softPower)

# Convert adjacency to topological overlap matrix (TOM)
TOM <- TOMsimilarity(adjacency)

# Ensure that the dimnames for TOM are set correctly
gene_ids <- colnames(datExpr)  # gene_ids should now be the column names of datExpr
dimnames(TOM) <- list(gene_ids, gene_ids)

#---------------------------------------------------------------------
# Extract genes in the lightslateblue module
lightslateblue_module_genes <- gene_module_df %>%
  filter(module == "lightslateblue")

# Extract the gene IDs in the lightslateblue module
lightslateblue_gene_ids <- lightslateblue_module_genes$gene_id

# Subset the TOM for the lightslateblue genes
TOM_lightslateblue <- TOM[lightslateblue_gene_ids, lightslateblue_gene_ids]

# Check if TOM_lightslateblue is not empty
if (dim(TOM_lightslateblue)[1] == 0) {
  stop("The TOM subset is empty. Check the gene IDs and TOM matrix.")
}

# Convert the TOM to an edge list
threshold <- 0.1  # Set a threshold to filter edges by weight
TOM_lightslateblue[TOM_lightslateblue < threshold] <- 0

library(reshape2)
edge_list_lightslateblue <- melt(TOM_lightslateblue)
colnames(edge_list_lightslateblue) <- c("Source", "Target", "Weight")
edge_list_lightslateblue <- edge_list_lightslateblue[edge_list_lightslateblue$Weight > 0, ]

edge_list_lightslateblue <- edge_list_lightslateblue %>%
  filter(Weight > 0 & Weight < 1)

# Save the edge list to CSV
write.csv(edge_list_lightslateblue, "lightslateblue_edge_list.csv", row.names = FALSE)

#-----------------
#list of edges related to Cre07.g317250 as well as weights

filtered_edge_list_lightslateblue <- edge_list_lightslateblue %>%
  filter(Source == "Cre07.g317250" | Target == "Cre07.g317250")

filtered_edge_list_lightslateblue <- filtered_edge_list_lightslateblue %>%
  arrange(desc(Weight))

write.csv(filtered_edge_list_lightslateblue, "filtered_lightslateblue_edge_list.csv", row.names = FALSE)

#--------------------------------------------------------------------
#Cre06.g270500 lightcyan module 
#list of edges related to Cre06.g270500 as well as weights

lightcyan_module_genes <- gene_module_df %>%
  filter(module == "lightcyan")

lightcyan_gene_ids <- lightcyan_module_genes$gene_id

TOM_lightcyan <- TOM[lightcyan_gene_ids, lightcyan_gene_ids]

if (dim(TOM_lightcyan)[1] == 0) {
  stop("The TOM subset is empty. Check the gene IDs and TOM matrix.")
}


threshold <- 0.1
TOM_lightcyan[TOM_lightcyan < threshold] <- 0

library(reshape2)
edge_list_lightcyan <- melt(TOM_lightcyan)
colnames(edge_list_lightcyan) <- c("Source", "Target", "Weight")
edge_list_lightcyan <- edge_list_lightcyan[edge_list_lightcyan$Weight > 0, ]

edge_list_lightcyan <- edge_list_lightcyan %>%
  filter(Weight > 0 & Weight < 1)

write.csv(edge_list_lightcyan, "lightcyan_edge_list.csv", row.names = FALSE)

#list of edges related to Cre06.g270500 as well as weights

filtered_edge_list_lightcyan <- edge_list_lightcyan %>%
  filter(Source == "Cre06.g270500" | Target == "Cre06.g270500")

filtered_edge_list_lightcyan <- filtered_edge_list_lightcyan %>%
  arrange(desc(Weight))

write.csv(filtered_edge_list_lightcyan, "filtered_lightcyan_edge_list.csv", row.names = FALSE)

#---------------------------------------------------------------------
#Cre06.g273100 grey60 module

grey60_module_genes <- gene_module_df %>%
  filter(module == "grey60")

grey60_gene_ids <- grey60_module_genes$gene_id

TOM_grey60 <- TOM[grey60_gene_ids, grey60_gene_ids]

if (dim(TOM_grey60)[1] == 0) {
  stop("The TOM subset is empty. Check the gene IDs and TOM matrix.")
}

threshold <- 0.1
TOM_grey60n[TOM_grey60 < threshold] <- 0

library(reshape2)
edge_list_grey60 <- melt(TOM_grey60)
colnames(edge_list_grey60) <- c("Source", "Target", "Weight")
edge_list_grey60 <- edge_list_grey60[edge_list_grey60$Weight > 0, ]

edge_list_grey60 <- edge_list_grey60 %>%
  filter(Weight > 0 & Weight < 1)

write.csv(edge_list_grey60, "grey60_edge_list.csv", row.names = FALSE)

#list of edges related to Cre06.g273100 as well as weights

filtered_edge_list_grey60 <- edge_list_grey60 %>%
  filter(Source == "Cre06.g273100" | Target == "Cre06.g273100")

filtered_edge_list_grey60 <- filtered_edge_list_grey60 %>%
  arrange(desc(Weight))

write.csv(filtered_edge_list_grey60, "filtered_grey60_edge_list.csv", row.names = FALSE)
