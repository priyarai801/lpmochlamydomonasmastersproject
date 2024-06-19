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

#in GSE42035_fpkm_data I have 17,741 genes
#I think next step is to figure out which columns I actually need
#the FPKMs are labelled with SRRs

#in gse42035_fpkm_data I want columns 4 (gene_id), 10 (SRR611223_FPKM CC-124 Light oxic)
#SRR611224 CC-124 Dark anoxic 0.5 hours column 50 and SRR611225 Dark anoxic 6 hours column 14

library(dplyr)

analysis_gse42035_fpkm_data <- gse42035_fpkm_data %>%
  select(4, 10, 50, 14)

write.csv(analysis_gse42035_fpkm_data, "analysis_gse42035_fpkm_data.csv", row.names = FALSE)

#below code is just the genes of interest to see if it's suitable
#to have a FPKM cutoff value of 1

genes_of_interest <- c("Cre07.g317250", "Cre06.g270500", "Cre06.g273100")

filtered_genes <- analysis_gse42035_fpkm_data %>%
  filter(gene_id %in% genes_of_interest)

head(filtered_genes)

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

# Step 1: Calculate the mean FPKM for each gene across all samples
zero_filtered_gse42035_fpkm_data$mean_fpkm <- rowMeans(zero_filtered_gse42035_fpkm_data[ , -1], na.rm = TRUE)

# Step 2: Calculate the quantiles of the mean FPKM
fpkm_quantiles <- quantile(zero_filtered_gse42035_fpkm_data$mean_fpkm, probs = seq(0, 1, 0.01))

# View the quantiles to choose a threshold
print(fpkm_quantiles)

#Looks like threshold of 15% is okay = 0.6692659
# Step 3: Choose a quantile threshold
threshold <- fpkm_quantiles["15%"]

# Step 4: Filter the genes using the chosen threshold
filtered_data <- zero_filtered_gse42035_fpkm_data %>%
  filter(mean_fpkm > threshold)

# Remove the mean_fpkm column as it's no longer needed
filtered_data <- filtered_data %>%
  select(-mean_fpkm)

# View the filtered dataframe
dim(filtered_data)

#Genes reduced to 14750 - still way too big

#---------------------------------------------------------------------
#CODE RAN UP TO HERE

# Load necessary libraries
library(dplyr)

# Assuming filtered_data is already loaded

# Step 1: Calculate the variance for each gene across all samples
filtered_data <- filtered_data %>%
  rowwise() %>%
  mutate(variance = var(c_across(-gene_id)))

# Step 2: Select the top 1,000 genes with the highest variance
top_genes <- filtered_data %>%
  arrange(desc(variance)) %>%
  slice(1:1000)

# Step 3: Ensure inclusion of genes of interest
genes_of_interest <- c("Cre07.g317250", "Cre06.g270500", "Cre06.g273100")
top_genes <- filtered_data %>%
  filter(gene_id %in% genes_of_interest | gene_id %in% top_genes$gene_id) %>%
  select(-variance)  # Remove the variance column if not needed

# View the top genes
head(top_genes)
dim(top_genes)

# Prepare data for WGCNA
datExpr <- as.data.frame(t(top_genes[,-1]))
rownames(datExpr) <- top_genes$gene_id

# Check the final dimensions
dim(datExpr)

# Write the filtered dataframe to a CSV file
write.csv(top_genes, "top_var_filtered_gse42035_fpkm_data.csv", row.names = FALSE)


