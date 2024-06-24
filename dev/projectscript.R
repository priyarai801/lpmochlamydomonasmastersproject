getwd()
setwd("/Users/priyarai/Documents/Researchproject/data")

#get all three files and save it into separate vectors, but only the column 
#with the FPKM value for CC-124 Light oxic column 4 of the original files

#Read in the data and save as dataframe
anaerobiosiscre07g317250 <- read.table("anaerobiosisfpkmcre07g317250coexpressed.txt", header = TRUE, sep = "\t")


anaerobiosiscre06g270500 <- read.table("anaerobiosisfpkmcre06g270500coexpressed.txt", header = TRUE, sep = "\t")


anaerobiosiscre06g273100 <- read.table("anaerobiosisfpkmcre06g273100coexpressed.txt", header = TRUE, sep = "\t")

#for wgcna there needs to be a weight - i assume wgcna will calculate this
#for a wgcna network for differentially expressed genes - read the 
#thermotolerant genes paper

#they first did pca to see if there were different clusters


# X (1st column) + CC.124.Light.oxic (5th column) + CC.124.Dark.anoxic.6.hours (4th column)
library(dplyr)
selectedcolumnsanaerobiosiscre07g317250 <- anaerobiosiscre07g317250 %>% select(1, 5, 3, 4)
selectedcolumnsanaerobiosiscre06g270500 <- anaerobiosiscre06g270500 %>% select(1, 5, 3, 4)
selectedcolumnsanaerobiosiscre06g273100 <- anaerobiosiscre06g273100 %>% select(1, 5, 3, 4)

#Merge 3 dataframes with multiple columns
merged_df <- bind_rows(selectedcolumnsanaerobiosiscre07g317250, selectedcolumnsanaerobiosiscre06g270500, selectedcolumnsanaerobiosiscre06g273100)
#in merged_df idk why the values are the e+, i guess there must have been 
#very small values but 7.6e+00 is legit the same thing as 7.6
#but values like 8.7e+01 means 87


#Get rid of duplicates
merged_df_unique <- merged_df %>% distinct(X, .keep_all = TRUE)
#Reduced from 153 to 141 genes

#---------------------------------------------------------------------
#Make PCA to see if clusters for cc-124 light oxic and cc-124 dark anoxic 6 hours are different 
#NOTE CODE DOES NOT WORK - PC1 and PC2 values keep appearing the same for different conditions

#I need the columns that refers to CC-124 Light Oxic and CC-124 Dark anoxic hours
data_for_pca <- merged_df_unique %>% select('CC.124.Light.oxic', 'CC.124.Dark.anoxic.6.hours')

#Does the PCA - it made a list of 5?
pca_comparelightoxicanddarkanoxic <- prcomp(data_for_pca, scale. = TRUE)

#Stores PCA scores with PC1 and PC2 into a dataframe
pca_scores <- as.data.frame(pca_comparelightoxicanddarkanoxic$x)

# PCA scores are integrated back with the original dataframe so that you
#don't lost the context aka creates dataframe of PC1, PC2, X(gene), Light oxic, Dark Oxic
pca_scores_merged_df_unique <- cbind(pca_scores, merged_df_unique)


#In order to use ggplot2 to plot PCA, u have to convert the data from a wide to long format
#By this, I mean a long format is where instead of light and dark being in 2 different columns
#u have one column that says whether the gene expression value belongs to light or dark
plot_data_long <- pca_scores_merged_df_unique %>%
  pivot_longer(cols = c('CC.124.Light.oxic', 'CC.124.Dark.anoxic.6.hours'), names_to = "Sample", values_to = "Expression")


# Plot PCA with colors based on the sample
library(ggplot2)
library(tidyr)
#Trying to plot the PCA - the points overlap with each other. I realised that
#the PC1 and PC2 values are the exact same for the same gene across the two different samples - why?

p <- ggplot(plot_data_long, aes(x = PC1, y = PC2, color = Sample)) +
  geom_point(size = 1) +
  labs(title = "PCA Plot Colored by Gene Expression in Sample1 and Sample2",
       x = "PC1",
       y = "PC2") +
  theme_minimal() +
  scale_color_manual(values = c("CC.124.Light.oxic" = "red", "CC.124.Dark.anoxic.6.hours" = "blue"))

print(p)

#---------------------------------------------------------------------
#Anyways next step of thermotolerance in wheat was carrying out Differential expression
#analysis using DESeq2 and the differential expression plots were visualised with volcano
#plots in the ggplot2 package also used ashr package to shrink the expression fold changes
#more info in paper) and had log2foldchanges and those with FDR-adjusted p-values were considered
#further for GO enrichment analysis

#will use merged_df_unique

#DESeq2 only works on raw counts but i have RPKM
#will use limma package for DE analysis

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
library(limma)
library(dplyr)

#Changes merged_df_unique to be in correct format for DE analysis
#sets rownames of merged_df_unique to be the gene names in X column
#Then the select-X will remove the X column and instead just have a column
#corresponding to the genes
#This is a useful step for me
rownames(merged_df_unique) <- merged_df_unique$X
merged_df_unique <- merged_df_unique %>% select(-X)

#Logtransform the RPKM gene expression values
log_rpkm <- log2(merged_df_unique + 1)

#Create a design matrix to represent the different conditions and helps with
#linear modelling i.e. DE analysis
condition <- factor(c(rep("Light", ncol(merged_df_unique) / 2), rep("Dark", ncol(merged_df_unique) / 2)))
design <- model.matrix(~ condition)

# Fit the linear model - error here
fit <- lmFit(log_rpkm, design)
fit <- eBayes(fit)

#Not sure how to fix error here - but i'll either try DESeq2

#---------------------------------------------------------------------
#Figure out how to do WGCNA

BiocManager::install("WGCNA")
library(WGCNA)

#there's an error with my WGCNA -it's not installing completely so i'm missing functions
#i can't even set soft threshold power - installing it via the tabs seemed
#to have maybe worked because i don't get the same error

install.packages("RSQLite")
library(RSQLite)


options(repos = c(CRAN = "https://cran.rstudio.com/",
                  BioCsoft = "https://bioconductor.org/packages/3.17/bioc",
                  BioCann = "https://bioconductor.org/packages/3.17/data/annotation",
                  BioCexp = "https://bioconductor.org/packages/3.17/data/experiment",
                  BioCworkflows = "https://bioconductor.org/packages/3.17/workflows"))

# Check current repositories
print(getOption("repos"))

# Install WGCNA
BiocManager::install("WGCNA")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.19")

#--------------------------------------------------------------------
#Trying WGCNA again

# Set row names to gene names
rownames(merged_df_unique) <- merged_df_unique$X
merged_df_unique <- merged_df_unique %>% select(-X)

log_rpkm <- log2(merged_df_unique + 1)

# Transpose data for WGCNA
datExpr <- as.data.frame(t(log_rpkm))

# Remove zero-variance genes
zero_variance_genes <- which(apply(datExpr, 2, var) == 0)
if (length(zero_variance_genes) > 0) {
  datExpr <- datExpr[, -zero_variance_genes]
}

# Remove rows with any NA values
datExpr <- na.omit(datExpr)


library(WGCNA)
library(dplyr)
library(ggplot2)
library(tidyr)

#This line creates potential candidate values for the soft-threshold
#power for the WGCNA so the options are 1-10 and 12,14,16,18,20
powers <- c(1:10, seq(from = 12, to = 20, by = 2))
#CODE HAS BEEN RUN UP TO THIS POINT

#With the candidate values in powers, below code will choose which
#value should be chosen as the soft threshold power
#the pickSoftThreshold function does this
#datExpr is the input data
#powerVector are the values in powers
#verbose is the level of info u want to get
#5 is the highest it tells you the connectivity + free topology fit index
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

#ERROR HERE??????????????? - this is what it says
# > sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# pickSoftThreshold: will use block size 141.
# pickSoftThreshold: calculating connectivity for given powers...
# ..working on genes 1 through 141 of 141
# Error in summary(lm1)$coefficients[2, 1] : subscript out of bounds


# Write datExpr to a CSV file
write.csv(datExpr, file = "datExpr.csv", row.names = TRUE)

#was it bcos the data was not numeric?
datExpr <- as.data.frame(lapply(datExpr, as.numeric))
#nope not it

#says to try doing it manually?
# Manually run the steps inside pickSoftThreshold
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

# Initialize variables to hold the fit indices
fitIndices <- matrix(0, nrow = length(powers), ncol = 4)
colnames(fitIndices) <- c("Power", "SFT.R.sq", "slope", "truncated.R.sq")

for (i in 1:length(powers)) {
  power <- powers[i]
  
  # Calculate adjacency
  adjacency <- adjacency(datExpr, power = power)
  
  # Check the adjacency matrix
  cat("Adjacency matrix for power:", power, "\n")
  print(head(adjacency))
  
  # Calculate connectivity
  k <- apply(adjacency, 1, sum) - 1
  
  # Check the connectivity values
  cat("Connectivity for power:", power, "\n")
  print(head(k))
  
  # Scale-free topology fit index
  fit <- scaleFreeFitIndex(k)
  
  # Store the fit indices
  fitIndices[i, 1] <- power
  fitIndices[i, 2] <- fit$Rsquared.SFT
  fitIndices[i, 3] <- fit$slope
  fitIndices[i, 4] <- fit$truncated.R.sq
}

# Print the fit indices
print(fitIndices)

# Select the optimal power
optimalPower <- fitIndices[which.max(fitIndices[, 2]), 1]
print(optimalPower)

#running the for i gives same error as with the sft code line 
#Error in summary(lm1)$coefficients[2, 1] : subscript out of bounds
#yeah the fitIndices should be not zero - that's where the error is 

#when running the for i, for the power 1 all the adjacency matrix values are 1 and
#the connectivity values are 140 - this is an error so 
#trying to see if the calculation for the correlation and adjacency matrix are wrong
#correlation values are 1 or -1 which should not be the case

str(datExpr)

# If your data frame needs to be transposed

datExpr <- as.data.frame(t(datExpr))


datExpr <- t(datExpr)
datExpr <- as.data.frame(datExpr)

# Check the structure again
str(datExpr)

#the SFT.R.sq values are too low, for the power 3, the value is 0.000378
#scale free topology value is too low should be close to 0
#or else it's not suitable for WGCNA
#say the way to overcome this is to increase the number of samples
#so either i include the other conditions
# or i can try starting from scratch and generating counts matrix myself

#--------------------------------------------------------------------
#i've got the GSE42035 FPKM tracking data from Phytozome

install.packages("readr")
library(readr)

gse42035_fpkm_data <- read_tsv("GSE42035_genes.fpkm_tracking")

head(gse42035_fpkm_data)
dim(gse42035_fpkm_data)

#in GSE42035_fpkm_data I have 17,741 genes
#I think next step is to figure out which columns I actually need
#the FPKMs are labelled with SRRs

#in gse42035_fpkm_data I want columns 4 (gene_id), 10 (SRR611223_FPKM CC-124 Light oxic)
#SRR611224 CC-124 Dark anoxic 0.5 hours column 50 and SRR611225 Dark anoxic 6 hours column 14

library(dplyr)

analysis_gse42035_fpkm_data <- gse42035_fpkm_data %>%
  select(4, 10, 50, 14)

head(analysis_gse42035_fpkm_data)
write.csv(analysis_gse42035_fpkm_data, "analysis_gse42035_fpkm_data.csv", row.names = FALSE)

#Code here is to log them which i should have done--------------------
#will probably have to change later dataframe names

# Log-transform the FPKM values (adding a small constant to avoid log(0))
analysis_gse42035_fpkm_data <- analysis_gse42035_fpkm_data %>%
  mutate(across(starts_with("SRR"), ~ log2(. + 1)))
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

#CODE RAN UP TO HERE 

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
