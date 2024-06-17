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
selectedcolumnsanaerobiosiscre07g317250 <- anaerobiosiscre07g317250 %>% select(1, 5, 4)
selectedcolumnsanaerobiosiscre06g270500 <- anaerobiosiscre06g270500 %>% select(1, 5, 4)
selectedcolumnsanaerobiosiscre06g273100 <- anaerobiosiscre06g273100 %>% select(1, 5, 4)

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
#the PC! and PC2 values are the exact same for the same gene across the two different samples - why?

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
#---------------------------------------------------------------------

# Transpose the data for WGCNA so genes are now columns, and conditions are rows
#so that it is easier for correlation calculations
datExpr <- as.data.frame(t(log_rpkm))


# Choose a set of soft-thresholding powers
powers = c(1:10, seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results:
par(mfrow = c(1,2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.90, col = "red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", 
     type = "n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1, col = "red")

softPower = sft$powerEstimate

# Create adjacency matrix
adjacency = adjacency(datExpr, power = softPower)

# Turn adjacency into Topological Overlap Matrix (TOM)
TOM = TOMsimilarity(adjacency)
dissTOM = 1 - TOM

# Perform hierarchical clustering
geneTree = hclust(as.dist(dissTOM), method = "average")

# Module identification using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = 30)

# Convert labels to colors for plotting
dynamicColors = labels2colors(dynamicMods)

# Plot the dendrogram and the module colors
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut")

# Define the trait data
trait_data <- data.frame(
  Light = c(rep(1, 3), rep(0, 3)),
  Dark = c(rep(0, 3), rep(1, 3))
)
rownames(trait_data) <- colnames(datExpr)

# Module eigengenes
MEs = moduleEigengenes(datExpr, colors = dynamicColors)$eigengenes

# Correlate modules with traits
moduleTraitCor = cor(MEs, trait_data, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))

# Plot the correlations
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(trait_data),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#--------------------------------------------------------------------
#Trying WGCNA again

# Set row names to gene names
rownames(merged_df_unique) <- merged_df_unique$X
merged_df_unique <- merged_df_unique %>% select(-X)

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

#going to write out datExpr as a file and attach it if it can find
#anything wrong

# Write datExpr to a CSV file
write.csv(datExpr, file = "datExpr.csv", row.names = TRUE)
