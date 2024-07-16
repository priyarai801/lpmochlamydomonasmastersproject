getwd()
setwd("/Users/priyarai/Documents/Researchproject/data/nconcentration")

#--------------------------------------------------------------------
#For WGCNA

library(WGCNA)
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)

#i've got the GSE34585 FPKM tracking data from Phytozome

gse34585_fpkm_data <- read_tsv("GSE34585_genes.fpkm_tracking")

head(gse34585_fpkm_data)
dim(gse34585_fpkm_data)

#in gse34585_fpkm_data I have 17,741 genes

#in gse34585_fpkm_data I want columns: 
# 4 (gene_id)

# 82 SRR393780_FPKM Nitrogen deficiency t=0 (1)
# 34 SRR393781_FPKM Nitrogen deficiency t=0 (2)
# 
# 86 SRR393782_FPKM Nitrogen deficiency t=30m (1)
# 38 SRR393783_FPKM Nitrogen deficiency t=30m (2)
# 
# 90 SRR393784_FPKM Nitrogen deficiency t=4h (1)
# 42 SRR393785_FPKM Nitrogen deficiency t=4h (2)
# 
# 94 SRR393786_FPKM Nitrogen deficiency t=8h(1)
# 46 SRR393787_FPKM Nitrogen deficiency t=8h (2)

analysis_gse34585_fpkm_data <- gse34585_fpkm_data %>%
  select(gene_id = 4, SRR393780_FPKM = 82, SRR393781_FPKM = 34, SRR393782_FPKM = 86, SRR393783_FPKM = 38, 
         SRR393784_FPKM = 90, SRR393785_FPKM = 42, SRR393786_FPKM = 94, SRR393787_FPKM = 46)

head(analysis_gse34585_fpkm_data)
dim(analysis_gse34585_fpkm_data)
#17741 9

write.csv(analysis_gse34585_fpkm_data, "analysis_gse34585_fpkm_data.csv", row.names = FALSE)

#---------------------------------------------------------------------
# Log-transform the FPKM values (adding a small constant to avoid log(0))
analysis_gse34585_log_fpkm_data <- analysis_gse34585_fpkm_data %>%
  mutate(across(starts_with("SRR"), ~ log2(. + 1)))

head(analysis_gse34585_log_fpkm_data)
dim(analysis_gse34585_log_fpkm_data)
#17741 9
#---------------------------------------------------------------------

#Just a table of genes of interest

genes_of_interest <- c("Cre07.g317250", "Cre06.g270500", "Cre06.g273100")

filtered_genes <- analysis_gse34585_log_fpkm_data %>%
  filter(gene_id %in% genes_of_interest)

head(filtered_genes)
dim(filtered_genes)

#---------------------------------------------------------------------
#create plot of filtered genes

# Calculate averages and standard errors - they're not yet defined 
filtered_genes <- filtered_genes %>%
  rowwise() %>%
  mutate(
    avg_0h = mean(c(SRR393780_FPKM, SRR393781_FPKM)),
    avg_0_5h = mean(c(SRR393782_FPKM, SRR393783_FPKM)),
    avg_4h = mean(c(SRR393784_FPKM, SRR393785_FPKM)),
    avg_8h = mean(c(SRR393786_FPKM, SRR393787_FPKM)),
    se_0h = sd(c(SRR393780_FPKM, SRR393781_FPKM)) / sqrt(2),
    se_0_5h = sd(c(SRR393782_FPKM, SRR393783_FPKM)) / sqrt(2),
    se_4h = sd(c(SRR393784_FPKM, SRR393785_FPKM)) / sqrt(2),
    se_8h = sd(c(SRR393786_FPKM, SRR393787_FPKM)) / sqrt(2)
  )

# Reshape data for plotting
plot_filtered_genes <- filtered_genes %>%
  select(gene_id, avg_0h, avg_0_5h, avg_4h, avg_8h, se_0h, se_0_5h, se_4h, se_8h) %>%
  pivot_longer(cols = starts_with("avg_"), names_to = "timepoint", values_to = "average") %>%
  pivot_longer(cols = starts_with("se_"), names_to = "se_timepoint", values_to = "se") %>%
  filter(
    (timepoint == "avg_0h" & se_timepoint == "se_0h") |
      (timepoint == "avg_0_5h" & se_timepoint == "se_0_5h") |
      (timepoint == "avg_4h" & se_timepoint == "se_4h") |
      (timepoint == "avg_8h" & se_timepoint == "se_8h")
  )

# Plot the data
ggplot(plot_filtered_genes, aes(x = gene_id, y = average, fill = timepoint)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = average - se, ymax = average + se), position = position_dodge(1), width = 0.25) +
  labs(
    title = "Expression of LPMO-containing Genes",
    x = "Gene ID",
    y = "Log2(FPKM + 1) Gene Expression",
    fill = "Timepoint"
  ) +
  theme_minimal()

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
print(names(analysis_gse34585_log_fpkm_data))

# Define the columns of interest based on the actual column names
columns_of_interest <- c("SRR393780_FPKM", "SRR393781_FPKM", "SRR393782_FPKM", "SRR393783_FPKM", "SRR393784_FPKM", 
                         "SRR393785_FPKM", "SRR393786_FPKM", "SRR393787_FPKM")

# Verify that the columns exist in the dataframe
missing_columns <- setdiff(columns_of_interest, names(analysis_gse34585_log_fpkm_data))
if (length(missing_columns) > 0) {
  stop("The following columns are missing in the dataframe: ", paste(missing_columns, collapse = ", "))
}

# Convert the specified columns to numeric
analysis_gse34585_log_fpkm_data[columns_of_interest] <- lapply(analysis_gse34585_log_fpkm_data[columns_of_interest], as.numeric)

# Create a logical condition to identify rows where all specified columns have a value of 0
rows_with_all_zeros <- apply(analysis_gse34585_log_fpkm_data[columns_of_interest], 1, function(row) all(row == 0))

# Print the number of rows identified
cat("Number of rows with all zeros in the specified columns:", sum(rows_with_all_zeros), "\n")
# Number of rows with all zeros in the specified columns: 529

# Filter out the rows with all zeros
analysis_gse34585_log_fpkm_data <- analysis_gse34585_log_fpkm_data[!rows_with_all_zeros, ]

cat("Dimensions after filtering:", dim(analysis_gse34585_log_fpkm_data), "\n")
#Dimensions after filtering: 17212 9

# Ensure gene_id is set as row names
analysis_gse34585_log_fpkm_data <- analysis_gse34585_log_fpkm_data %>% column_to_rownames("gene_id")

# Check the structure of the data
head(analysis_gse34585_log_fpkm_data)

# Transpose the data for WGCNA (genes as columns and samples as rows)s
datExpr <- as.data.frame(t(analysis_gse34585_log_fpkm_data))

# Check the transposed data structure
str(datExpr)
dim(datExpr)

#8 17212

powers <- c(1:10, seq(from = 12, to = 50, by = 2))

sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
#Takes about a minute

write.csv(analysis_gse34585_log_fpkm_data, "analysis_gse34585_log_fpkm_data.csv", row.names = TRUE)

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

#Going with power 8 with SFT.R.sq as 0.93, first above 0.9

# Set the soft threshold power
softPower <- 8

# Construct the adjacency matrix
adjacency <- adjacency(datExpr, power = softPower)
#This line takes about a minute
dim(adjacency)
#17212 17212
#calculates it for 17212 genes x 17212 genes hence why there are
#nearly 300,000,000 elements

# Turn adjacency into topological overlap matrix (TOM)
TOM <- TOMsimilarity(adjacency)
dim(TOM)
#17212 17212
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
#set cutHeight to 0.985 ===> 99% of the (truncated) height range in dendro.

table(dynamicMods)
#There are 197 modules

#assigns colours to the modules
dynamicColors <- labels2colors(dynamicMods)

# Plot the dendrogram and the module colors
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)


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
analysis_gse34585_log_fpkm_data <- analysis_gse34585_log_fpkm_data %>% column_to_rownames("gene_id")

# Check the structure of the data
head(analysis_gse34585_log_fpkm_data)

# Transpose the data for WGCNA (genes as columns and samples as rows)
datExpr <- as.data.frame(t(analysis_gse34585_log_fpkm_data))
#---------------------------------------------------------------------

#Next step is to relate modules to external/phenotypic traits - i think
#this is an optional step so unsure at this point to apply this to my
#network

#---------------------------------------------------------------------
#Trying to find out which module the LPMO genes are in

# Create a data frame with gene IDs and their corresponding module colors
gene_module_df <- data.frame(
  gene_id = rownames(analysis_gse34585_log_fpkm_data),
  module = dynamicColors
)

# Filter for genes of interest
genes_of_interest_modules <- gene_module_df %>%
  filter(gene_id %in% genes_of_interest)

# Print the results
print(genes_of_interest_modules)

#---------------------------------------------------------------------
#Cre07.g317250 brown
#Cre06.g270500 darkviolet
#Cre06.g273100 blue

#Create dataframe of all the genes in the dataset that are designated
#to the module of interest

#Cre07.g317250 brown module
brown_module_genes <- gene_module_df %>%
  filter(module == "brown")

#Add in the log FPKM values too
brown_module_genes_with_fpkm <- brown_module_genes %>%
  left_join(analysis_gse34585_log_fpkm_data %>% rownames_to_column("gene_id"), by = "gene_id")

write.csv(brown_module_genes_with_fpkm, "brown_module_genes_with_fpkm.csv", row.names = FALSE)

#---------------------------------------------------------------------
#Cre06.g270500 darkviolet
darkviolet_module_genes <- gene_module_df %>%
  filter(module == "darkviolet")

#Add in the log FPKM values too
darkviolet_module_genes_with_fpkm <- darkviolet_module_genes %>%
  left_join(analysis_gse34585_log_fpkm_data %>% rownames_to_column("gene_id"), by = "gene_id")

write.csv(darkviolet_module_genes_with_fpkm, "darkviolet_module_genes_with_fpkm.csv", row.names = FALSE)

#---------------------------------------------------------------------
#Cre06.g273100 blue
blue_module_genes <- gene_module_df %>%
  filter(module == "blue")

#Add in the log FPKM values too
blue_module_genes_with_fpkm <- blue_module_genes %>%
  left_join(analysis_gse34585_log_fpkm_data %>% rownames_to_column("gene_id"), by = "gene_id")

write.csv(blue_module_genes_with_fpkm, "blue_module_genes_with_fpkm.csv", row.names = FALSE)

#---------------------------------------------------------------------
#Just to get a list of gene_ids for a quick GO enrichment analysis

#Cre07.g317250 brown module GO
brown_module_gene_ids <- gene_module_df %>%
  filter(module == "brown") %>%
  pull(gene_id)

brown_module_gene_ids <- data.frame(gene_id = brown_module_gene_ids)

write.csv(brown_module_gene_ids, "brown_module_genes_ids.csv", row.names = FALSE)

#---------------------------------------------------------------------
#Cre06.g270500 darkviolet GO

darkviolet_module_gene_ids <- gene_module_df %>%
  filter(module == "darkviolet") %>%
  pull(gene_id)

darkviolet_module_gene_ids <- data.frame(gene_id = darkviolet_module_gene_ids)

write.csv(darkviolet_module_gene_ids, "darkviolet_module_genes_ids.csv", row.names = FALSE)

#---------------------------------------------------------------------
#Cre06.g273100 blue GO

blue_module_gene_ids <- gene_module_df %>%
  filter(module == "blue") %>%
  pull(gene_id)

blue_module_gene_ids <- data.frame(gene_id = blue_module_gene_ids)

write.csv(blue_module_gene_ids, "blue_module_genes_ids.csv", row.names = FALSE)

#---------------------------------------------------------------------
#get list of edges related to Cre07.g317250 brown module

#---------------------------------------------------------------------
# Ensure gene_id is set as row names
analysis_gse34585_log_fpkm_data <- analysis_gse34585_log_fpkm_data %>% 
  column_to_rownames("gene_id")

# Transpose the data for WGCNA (genes as columns and samples as rows)
datExpr <- as.data.frame(t(analysis_gse34585_log_fpkm_data))

# Check the structure of the transposed data
str(datExpr)
dim(datExpr)

# Construct adjacency matrix
softPower <- 8
adjacency <- adjacency(datExpr, power = softPower)

# Convert adjacency to topological overlap matrix (TOM)
TOM <- TOMsimilarity(adjacency)

# Ensure that the dimnames for TOM are set correctly
gene_ids <- colnames(datExpr)  # gene_ids should now be the column names of datExpr
dimnames(TOM) <- list(gene_ids, gene_ids)

# Extract genes in the brown module
brown_module_genes <- gene_module_df %>%
  filter(module == "brown")

# Extract the gene IDs in the module
brown_gene_ids <- brown_module_genes$gene_id

# Subset the TOM for the genes
TOM_brown <- TOM[brown_gene_ids, brown_gene_ids]

# Check if TOM_brown is not empty
if (dim(TOM_brown)[1] == 0) {
  stop("The TOM subset is empty. Check the gene IDs and TOM matrix.")
}

# Convert the TOM to an edge list
threshold <- 0.1  # Set a threshold to filter edges by weight
TOM_brown[TOM_brown < threshold] <- 0

library(reshape2)
edge_list_brown <- melt(TOM_brown)
colnames(edge_list_brown) <- c("Source", "Target", "Weight")
edge_list_brown <- edge_list_brown[edge_list_brown$Weight > 0, ]

edge_list_brown <- edge_list_brown %>%
  filter(Weight > 0 & Weight < 1)

# Save the edge list to CSV
write.csv(edge_list_brown, "brown_edge_list.csv", row.names = FALSE)

#-----------------
#list of edges related to Cre07.g317250 as well as weights

filtered_edge_list_brown <- edge_list_brown %>%
  filter(Source == "Cre07.g317250" | Target == "Cre07.g317250")

filtered_edge_list_brown <- filtered_edge_list_brown %>%
  arrange(desc(Weight))

write.csv(filtered_edge_list_brown, "filtered_brown_edge_list.csv", row.names = FALSE)

#---------------------------------------------------------------------
# Cre06.g270500 darkviolet module
darkviolet_module_genes <- gene_module_df %>%
  filter(module == "darkviolet")

# Extract the gene IDs in the module
darkviolet_gene_ids <- darkviolet_module_genes$gene_id

# Subset the TOM for the genes
TOM_darkviolet <- TOM[darkviolet_gene_ids, darkviolet_gene_ids]

# Check if TOM_darkviolet is not empty
if (dim(TOM_darkviolet)[1] == 0) {
  stop("The TOM subset is empty. Check the gene IDs and TOM matrix.")
}

# Convert the TOM to an edge list
threshold <- 0.1  # Set a threshold to filter edges by weight
TOM_darkviolet[TOM_darkviolet < threshold] <- 0

library(reshape2)
edge_list_darkviolet <- melt(TOM_darkviolet)
colnames(edge_list_darkviolet) <- c("Source", "Target", "Weight")
edge_list_darkviolet <- edge_list_darkviolet[edge_list_darkviolet$Weight > 0, ]

edge_list_darkviolet <- edge_list_darkviolet %>%
  filter(Weight > 0 & Weight < 1)

# Save the edge list to CSV
write.csv(edge_list_darkviolet, "darkviolet_edge_list.csv", row.names = FALSE)

#-----------------
#list of edges related to Cre06.g270500 as well as weights

filtered_edge_list_darkviolet <- edge_list_darkviolet %>%
  filter(Source == "Cre06.g270500" | Target == "Cre06.g270500")

filtered_edge_list_darkviolet <- filtered_edge_list_darkviolet %>%
  arrange(desc(Weight))

write.csv(filtered_edge_list_darkviolet, "filtered_darkviolet_edge_list.csv", row.names = FALSE)


#---------------------------------------------------------------------
#code to see what modules the genes from 2012 paper are in -
#REMEMBER to convert back the genes_of_interest to the 3 LPMO genes

# List of genes of interest
genes_of_interest <- c("Cre07.g317250", "Cre06.g270500", "Cre06.g273100")

genes_of_interest <- c("Cre07.g317250", "Cre06.g270500", "Cre06.g273100",
                       "Cre07.g320850", "Cre03.g171050", "Cre17.g730600",
                       "Cre17.g730550", "Cre08.g379450", "Cre19.g752200",
                       "Cre01.g048350", "Cre02.g098000", "Cre07.g343950",
                       "Cre12.g492800", "Cre07.g337750", "Cre03.g170700",
                       "Cre03.g173100", "Cre02.g115950", "Cre01.g026250",
                       "Cre12.g556350", "Cre03.g194700", "Cre03.g190500",
                       "Cre12.g488000", "Cre12.g488050", "Cre12.g501900",
                       "Cre12.g501850", "Cre12.g507029", "Cre13.g582300",
                       "Cre10.g437950", "Cre02.g141600", "Cre07.g336600",
                       "Cre06.g301600", "Cre13.g579750", "Cre12.g513400",
                       "Cre07.g314850")

# Filter for genes of interest
genes_of_interest_modules <- gene_module_df %>%
  filter(gene_id %in% genes_of_interest)

# Print the results
print(genes_of_interest_modules)

#---------------------------------------------------------------------
#code to make a heatmap to see how similar module eigengenes are
# Assuming you have the module eigengenes (MEs) calculated from WGCNA


#---------------------------------------------------------------------
#with the fastqfiles trying to diy

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("ShortRead", "Rsamtools"))

#---------------------------------------------------------------------
#added list of genes from phycocosm to see if these genes are also in
#the same modules as LPMO genes bcos GO annotation wasn't good
#at picking out LPMO-function related modules

genes_of_interest <- c(
  "Cre01.g002787", "Cre01.g026250", "Cre01.g031500", "Cre01.g044100", "Cre01.g048350",
  "Cre02.g087700", "Cre02.g091750", "Cre02.g095126", "Cre02.g096100", "Cre02.g098000",
  "Cre02.g098350", "Cre02.g112400", "Cre02.g115950", "Cre02.g141600", "Cre03.g146727",
  "Cre03.g146747", "Cre03.g146767", "Cre03.g148201", "Cre03.g151000", "Cre03.g152050",
  "Cre03.g155001", "Cre03.g158050", "Cre03.g163050", "Cre03.g163150", "Cre03.g170700",
  "Cre03.g171050", "Cre03.g173100", "Cre03.g181500", "Cre03.g189050", "Cre03.g190500",
  "Cre03.g194700", "Cre03.g195600", "Cre03.g200655", "Cre03.g207713", "Cre05.g233304",
  "Cre05.g243450", "Cre05.g245352", "Cre06.g267050", "Cre06.g269601", "Cre06.g269650",
  "Cre06.g270100", "Cre06.g270350", "Cre06.g270500", "Cre06.g273100", "Cre06.g278252",
  "Cre06.g282000", "Cre06.g283400", "Cre06.g285150", "Cre06.g289850", "Cre06.g301600",
  "Cre06.g306000", "Cre06.g307150", "Cre06.g307650", "Cre07.g314850", "Cre07.g314866",
  "Cre07.g317250", "Cre07.g319300", "Cre07.g320850", "Cre07.g332300", "Cre07.g336600",
  "Cre07.g337750", "Cre07.g338550", "Cre07.g339600", "Cre07.g343933", "Cre07.g343950",
  "Cre08.g362450", "Cre08.g373450", "Cre08.g379450", "Cre08.g384750", "Cre08.g385500",
  "Cre09.g386137", "Cre09.g387100", "Cre09.g393765", "Cre09.g394473", "Cre09.g394510",
  "Cre09.g394547", "Cre09.g396451", "Cre09.g401886", "Cre09.g407501", "Cre09.g415600",
  "Cre10.g437950", "Cre10.g444700", "Cre10.g447550", "Cre10.g450500", "Cre10.g451600",
  "Cre10.g455950", "Cre10.g456000", "Cre10.g456050", "Cre10.g456100", "Cre10.g457500",
  "Cre10.g458350", "Cre11.g467538", "Cre11.g467539", "Cre11.g467540", "Cre11.g467779",
  "Cre11.g476650", "Cre11.g478184", "Cre12.g488000", "Cre12.g488050", "Cre12.g492750",
  "Cre12.g492800", "Cre12.g492851", "Cre12.g501850", "Cre12.g501900", "Cre12.g507029",
  "Cre12.g507051", "Cre12.g509200", "Cre12.g513400", "Cre12.g514200", "Cre12.g551200",
  "Cre12.g556350", "Cre13.g570700", "Cre13.g577300", "Cre13.g579582", "Cre13.g579734",
  "Cre13.g579750", "Cre13.g582250", "Cre13.g582270", "Cre13.g582300", "Cre13.g587600",
  "Cre14.g631300", "Cre16.g653350", "Cre16.g657200", "Cre16.g657250", "Cre16.g666334",
  "Cre16.g667100", "Cre16.g677300", "Cre16.g693950", "Cre17.g696900", "Cre17.g698850",
  "Cre17.g703000", "Cre17.g719900", "Cre17.g730550", "Cre17.g730600", "Cre17.g732350",
  "Cre17.g732600", "Cre19.g752200")

# Filter for genes of interest
genes_of_interest_modules <- gene_module_df %>%
  filter(gene_id %in% genes_of_interest)

print(genes_of_interest_modules)

write.csv(genes_of_interest_modules, "genes_of_interest_modules_137.csv", row.names = FALSE)

# Filter out the genes that belong to the specified modules
filtered_genes_of_interest_modules <- genes_of_interest_modules %>%
  filter(module %in% c("blue", "brown", "darkviolet"))

print(filtered_genes_of_interest_modules)

write.csv(filtered_genes_of_interest_modules, "filtered_genes_of_interest_modules.csv", row.names = FALSE)
