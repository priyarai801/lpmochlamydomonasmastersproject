
#Load the following libraries for analysis
library(WGCNA)
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)


#Read in the FPKM data obtained from Phytozome
gse34585_fpkm_data <- read_tsv("GSE34585_genes.fpkm_tracking")

head(gse34585_fpkm_data)
dim(gse34585_fpkm_data)

#The gse34585_fpkm_data contains 17,741 genes

#In gse34585_fpkm_data I want columns: 
# columns 4 (gene_id)

# columns 82 SRR393780_FPKM Nitrogen deficiency t=0 (1)
# columns 34 SRR393781_FPKM Nitrogen deficiency t=0 (2)
# 
# columns 86 SRR393782_FPKM Nitrogen deficiency t=30m (1)
# columns 38 SRR393783_FPKM Nitrogen deficiency t=30m (2)
# 
# columns 90 SRR393784_FPKM Nitrogen deficiency t=4h (1)
# columns 42 SRR393785_FPKM Nitrogen deficiency t=4h (2)
# 
# columns 94 SRR393786_FPKM Nitrogen deficiency t=8h(1)
# columns 46 SRR393787_FPKM Nitrogen deficiency t=8h (2)

analysis_gse34585_fpkm_data <- gse34585_fpkm_data %>%
  select(gene_id = 4, SRR393780_FPKM = 82, SRR393781_FPKM = 34, SRR393782_FPKM = 86, SRR393783_FPKM = 38, 
         SRR393784_FPKM = 90, SRR393785_FPKM = 42, SRR393786_FPKM = 94, SRR393787_FPKM = 46)

head(analysis_gse34585_fpkm_data)
dim(analysis_gse34585_fpkm_data)
#17741 9

write.csv(analysis_gse34585_fpkm_data, "analysis_gse34585_fpkm_data.csv", row.names = FALSE)

#---------------------------------------------------------------------
# Log-transform the FPKM values
analysis_gse34585_log_fpkm_data <- analysis_gse34585_fpkm_data %>%
  mutate(across(starts_with("SRR"), ~ log2(. + 1)))

head(analysis_gse34585_log_fpkm_data)
dim(analysis_gse34585_log_fpkm_data)
#17741 9
#---------------------------------------------------------------------

#Create a vector containing just the LPMO-containing proteins

genes_of_interest <- c("Cre07.g317250", "Cre06.g270500", "Cre06.g273100")

filtered_genes <- analysis_gse34585_log_fpkm_data %>%
  filter(gene_id %in% genes_of_interest)

head(filtered_genes)
dim(filtered_genes)

#---------------------------------------------------------------------
#Create plot to look at gene expression levels of the LPMO-containing proteins

# Calculate the average and standard error values
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

plot_filtered_genes <- plot_filtered_genes %>%
  mutate(timepoint = factor(timepoint, levels = c("avg_0h", "avg_0_5h", "avg_4h", "avg_8h"),
                            labels = c("0h", "0.5h", "4h", "8h")))

#Code to plot the data
ggplot(plot_filtered_genes, aes(x = gene_id, y = average, fill = timepoint)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = average - se, ymax = average + se), position = position_dodge(1), width = 0.25) +
  labs(
    title = expression("Expression Of LPMO-Containing Genes in " * italic("C. reinhardtii")),
    x = "Gene",
    y = "Normalised Gene Expression (Log FPKM)",
    fill = "Time After Nitrogen Deprivation"
  ) +
  theme_minimal()

#---------------------------------------------------------------------
#Code to carry out pairwise t-test testing


perform_pairwise_t_test <- function(gene_data, gene_name) {
  expression_values <- unlist(gene_data[2:9])
  time_points <- factor(rep(c("0h", "0h", "0.5h", "0.5h", "4h", "4h", "8h", "8h"), times = 1))
  
  pairwise_results <- pairwise.t.test(expression_values, time_points, p.adjust.method = "bonferroni")
  
  cat("Pairwise t-test results for", gene_name, ":\n")
  print(pairwise_results)
}

#Conducts parwise t-test on the LPMO-containing genes
perform_pairwise_t_test(filtered_genes[filtered_genes$gene_id == "Cre06.g270500", ], "Cre06.g270500")
perform_pairwise_t_test(filtered_genes[filtered_genes$gene_id == "Cre06.g273100", ], "Cre06.g273100")
perform_pairwise_t_test(filtered_genes[filtered_genes$gene_id == "Cre07.g317250", ], "Cre07.g317250")


#---------------------------------------------------------------------
#Load libraries for WGCNA analysis if not already done so

library(WGCNA)
library(dplyr)
library(readr)
library(tibble)

#----------------------------------------------------------------------
#Begin WGCNA Analysis by removing rows with 0 values

#Check column names in the analysis_gse34585_log_fpkm_data
print(names(analysis_gse34585_log_fpkm_data))


columns_of_interest <- c("SRR393780_FPKM", "SRR393781_FPKM", "SRR393782_FPKM", "SRR393783_FPKM", "SRR393784_FPKM", 
                         "SRR393785_FPKM", "SRR393786_FPKM", "SRR393787_FPKM")


missing_columns <- setdiff(columns_of_interest, names(analysis_gse34585_log_fpkm_data))
if (length(missing_columns) > 0) {
  stop("The following columns are missing in the dataframe: ", paste(missing_columns, collapse = ", "))
}


analysis_gse34585_log_fpkm_data[columns_of_interest] <- lapply(analysis_gse34585_log_fpkm_data[columns_of_interest], as.numeric)


rows_with_all_zeros <- apply(analysis_gse34585_log_fpkm_data[columns_of_interest], 1, function(row) all(row == 0))


cat("Number of rows with all zeros in the specified columns:", sum(rows_with_all_zeros), "\n")
# Number of rows with all zeros in the specified columns: 529

# Filter out the rows with all zeros
analysis_gse34585_log_fpkm_data <- analysis_gse34585_log_fpkm_data[!rows_with_all_zeros, ]

cat("Dimensions after filtering:", dim(analysis_gse34585_log_fpkm_data), "\n")
#Dimensions after filtering: 17212 9


analysis_gse34585_log_fpkm_data <- analysis_gse34585_log_fpkm_data %>% column_to_rownames("gene_id")


head(analysis_gse34585_log_fpkm_data)


datExpr <- as.data.frame(t(analysis_gse34585_log_fpkm_data))


str(datExpr)
dim(datExpr)
#8 17212

#---------------------------------------------------------------------
#Continue WGCNA analysis by selecting a soft threshold power

powers <- c(1:10, seq(from = 12, to = 50, by = 2))

sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)


write.csv(analysis_gse34585_log_fpkm_data, "analysis_gse34585_log_fpkm_data.csv", row.names = TRUE)

#--------------------------------------------------------------------
# Set the soft threshold power as 8
softPower <- 8

# Create the adjacency matrix
adjacency <- adjacency(datExpr, power = softPower)
dim(adjacency)
#17212 17212


# Convert the adjacency matrix into topological overlap matrix (TOM)
TOM <- TOMsimilarity(adjacency)
dim(TOM)
#17212 17212

#Create a dissimilarity TOM matrix
dissTOM <- 1 - TOM
#Takes about a minute


# Hierarchical clustering of the genes
geneTree <- hclust(as.dist(dissTOM), method = "average")

# Dynamic tree cut to identify modules
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE)
#set cutHeight to 0.985 ===> 99% of the (truncated) height range in dendro.

table(dynamicMods)
#There are 197 modules in the network

#assigns unique colours to each module
dynamicColors <- labels2colors(dynamicMods)


#---------------------------------------------------------------------
#Calculate the module eigengene which is a summary of the overall gene expression profile of all the genes
#within a module


MEs <- moduleEigengenes(datExpr, colors = dynamicColors)$eigengenes

# Calculate the dissimilarity of module eigengenes
MEDiss <- 1 - cor(MEs)

# Cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average")

#---------------------------------------------------------------------

#Code to find the colour of the modules associated with the LPMO-containing proteinsm


genes_of_interest <- c("Cre07.g317250", "Cre06.g270500", "Cre06.g273100")


analysis_gse34585_log_fpkm_data <- analysis_gse34585_log_fpkm_data %>% column_to_rownames("gene_id")

head(analysis_gse34585_log_fpkm_data)

datExpr <- as.data.frame(t(analysis_gse34585_log_fpkm_data))

# Create a data frame with gene IDs and their corresponding module colors
gene_module_df <- data.frame(
  gene_id = rownames(analysis_gse34585_log_fpkm_data),
  module = dynamicColors
)

# Filter for genes of interest
genes_of_interest_modules <- gene_module_df %>%
  filter(gene_id %in% genes_of_interest)

print(genes_of_interest_modules)

#Here are the colour of the modules of interest:
#Cre07.g317250 brown
#Cre06.g270500 darkviolet
#Cre06.g273100 blue


#Create a dataframe of all the genes in the brown module associated with the
#LPMO-containing protein Cre07.g317250 with the log10 FPKM values

brown_module_genes <- gene_module_df %>%
  filter(module == "brown")

brown_module_genes_with_fpkm <- brown_module_genes %>%
  left_join(analysis_gse34585_log_fpkm_data %>% rownames_to_column("gene_id"), by = "gene_id")

write.csv(brown_module_genes_with_fpkm, "brown_module_genes_with_fpkm.csv", row.names = FALSE)

#---------------------------------------------------------------------
#Create a dataframe of all the genes in the darkviolet module associated with the
#LPMO-containing protein Cre06.g270500 with the log10 FPKM values

darkviolet_module_genes <- gene_module_df %>%
  filter(module == "darkviolet")

darkviolet_module_genes_with_fpkm <- darkviolet_module_genes %>%
  left_join(analysis_gse34585_log_fpkm_data %>% rownames_to_column("gene_id"), by = "gene_id")

write.csv(darkviolet_module_genes_with_fpkm, "darkviolet_module_genes_with_fpkm.csv", row.names = FALSE)

#---------------------------------------------------------------------
#Create a dataframe of all the genes in the blue module associated with the
#LPMO-containing protein Cre06.g273100with the log10 FPKM values

blue_module_genes <- gene_module_df %>%
  filter(module == "blue")

blue_module_genes_with_fpkm <- blue_module_genes %>%
  left_join(analysis_gse34585_log_fpkm_data %>% rownames_to_column("gene_id"), by = "gene_id")

write.csv(blue_module_genes_with_fpkm, "blue_module_genes_with_fpkm.csv", row.names = FALSE)

#---------------------------------------------------------------------
#Code to obtain just a list of the genes in each module associated with
#the LPMO-containing proteins for GO Enrchment Analysis

#Cre07.g317250 brown module
brown_module_gene_ids <- gene_module_df %>%
  filter(module == "brown") %>%
  pull(gene_id)

brown_module_gene_ids <- data.frame(gene_id = brown_module_gene_ids)

write.csv(brown_module_gene_ids, "brown_module_genes_ids.csv", row.names = FALSE)

#Cre06.g270500 darkviolet GO

darkviolet_module_gene_ids <- gene_module_df %>%
  filter(module == "darkviolet") %>%
  pull(gene_id)

darkviolet_module_gene_ids <- data.frame(gene_id = darkviolet_module_gene_ids)

write.csv(darkviolet_module_gene_ids, "darkviolet_module_genes_ids.csv", row.names = FALSE)

#Cre06.g273100 blue GO

blue_module_gene_ids <- gene_module_df %>%
  filter(module == "blue") %>%
  pull(gene_id)

blue_module_gene_ids <- data.frame(gene_id = blue_module_gene_ids)

write.csv(blue_module_gene_ids, "blue_module_genes_ids.csv", row.names = FALSE)

#---------------------------------------------------------------------

#Code to get list of edges of each module associated with LPMO-containing protein
#to be visualised in Cytoscape

#Code to obtain edge list for Cre07.g317250
datExpr <- as.data.frame(t(analysis_gse34585_log_fpkm_data))

str(datExpr)
dim(datExpr)

softPower <- 8
adjacency <- adjacency(datExpr, power = softPower)

TOM <- TOMsimilarity(adjacency)

gene_ids <- colnames(datExpr)
dimnames(TOM) <- list(gene_ids, gene_ids)

# Extract genes in the brown module
brown_module_genes <- gene_module_df %>%
  filter(module == "brown")

# Extract the gene IDs in the brown module
brown_gene_ids <- brown_module_genes$gene_id

TOM_brown <- TOM[brown_gene_ids, brown_gene_ids]

if (dim(TOM_brown)[1] == 0) {
  stop("The TOM subset is empty. Check the gene IDs and TOM matrix.")
}

# Set a threshold to filter edges by weight as 0.1
threshold <- 0.1 
TOM_brown[TOM_brown < threshold] <- 0

library(reshape2)
edge_list_brown <- melt(TOM_brown)
colnames(edge_list_brown) <- c("Source", "Target", "Weight")
edge_list_brown <- edge_list_brown[edge_list_brown$Weight > 0, ]

edge_list_brown <- edge_list_brown %>%
  filter(Weight > 0 & Weight < 1)

write.csv(edge_list_brown, "brown_edge_list.csv", row.names = FALSE)

#---------------------------------------------------------------------
#Code to get edge list of edges only related to Cre07.g317250 in the module

filtered_edge_list_brown <- edge_list_brown %>%
  filter(Source == "Cre07.g317250" | Target == "Cre07.g317250")

filtered_edge_list_brown <- filtered_edge_list_brown %>%
  arrange(desc(Weight))

write.csv(filtered_edge_list_brown, "filtered_brown_edge_list.csv", row.names = FALSE)

#---------------------------------------------------------------------
#Code for getting just top 20 nodes with highest weight in Cre07.g317250 brown module

library(dplyr)

filtered_edge_list_brown_removed_duplicates <- filtered_edge_list_brown %>%
  mutate(Edge = paste(pmin(Source, Target), pmax(Source, Target), sep = "-")) %>%
  distinct(Edge, .keep_all = TRUE) %>%
  select(-Edge)

filtered_edge_list_brown_top_20 <- filtered_edge_list_brown_removed_duplicates %>%
  arrange(desc(Weight)) %>%
  slice(1:20)

filtered_edge_list_brown_top_20$Source <- as.character(filtered_edge_list_brown_top_20$Source)
filtered_edge_list_brown_top_20$Target <- as.character(filtered_edge_list_brown_top_20$Target)

write.csv(filtered_edge_list_brown_top_20, "filtered_edge_list_brown_top_20_nodes.csv", row.names = FALSE)


#---------------------------------------------------------------------
# #Code to obtain edge list for Cre06.g270500 darkviolet module
darkviolet_module_genes <- gene_module_df %>%
  filter(module == "darkviolet")


darkviolet_gene_ids <- darkviolet_module_genes$gene_id


TOM_darkviolet <- TOM[darkviolet_gene_ids, darkviolet_gene_ids]


if (dim(TOM_darkviolet)[1] == 0) {
  stop("The TOM subset is empty. Check the gene IDs and TOM matrix.")
}


threshold <- 0.1
TOM_darkviolet[TOM_darkviolet < threshold] <- 0

library(reshape2)
edge_list_darkviolet <- melt(TOM_darkviolet)
colnames(edge_list_darkviolet) <- c("Source", "Target", "Weight")
edge_list_darkviolet <- edge_list_darkviolet[edge_list_darkviolet$Weight > 0, ]

edge_list_darkviolet <- edge_list_darkviolet %>%
  filter(Weight > 0 & Weight < 1)


write.csv(edge_list_darkviolet, "darkviolet_edge_list.csv", row.names = FALSE)

#---------------------------------------------------------------------
#Code to get edge list of edges only related to Cre06.g270500 in the module

filtered_edge_list_darkviolet <- edge_list_darkviolet %>%
  filter(Source == "Cre06.g270500" | Target == "Cre06.g270500")

filtered_edge_list_darkviolet <- filtered_edge_list_darkviolet %>%
  arrange(desc(Weight))

write.csv(filtered_edge_list_darkviolet, "filtered_darkviolet_edge_list.csv", row.names = FALSE)

#---------------------------------------------------------------------
#Code for getting just top 20 nodes in Cre06.g270500 darkviolet module

filtered_edge_list_darkviolet_removed_duplicates <- filtered_edge_list_darkviolet %>%
  mutate(Edge = paste(pmin(Source, Target), pmax(Source, Target), sep = "-")) %>%
  distinct(Edge, .keep_all = TRUE) %>%
  select(-Edge)

filtered_edge_list_darkviolet_top_20 <- filtered_edge_list_darkviolet_removed_duplicates %>%
  arrange(desc(Weight)) %>%
  slice(1:20)

filtered_edge_list_darkviolet_top_20$Source <- as.character(filtered_edge_list_darkviolet_top_20$Source)
filtered_edge_list_darkviolet_top_20$Target <- as.character(filtered_edge_list_darkviolet_top_20$Target)

write.csv(filtered_edge_list_darkviolet_top_20, "filtered_edge_list_darkviolet_top_20_nodes.csv", row.names = FALSE)

#---------------------------------------------------------------------
## #Code to obtain edge list for blue module for Cre06.g273100

blue_module_genes <- gene_module_df %>%
  filter(module == "blue")

blue_gene_ids <- blue_module_genes$gene_id

TOM_blue <- TOM[blue_gene_ids, blue_gene_ids]

if (dim(TOM_blue)[1] == 0) {
  stop("The TOM subset is empty. Check the gene IDs and TOM matrix.")
}

threshold <- 0.1
TOM_blue[TOM_blue < threshold] <- 0

library(reshape2)
edge_list_blue <- melt(TOM_blue)
colnames(edge_list_blue) <- c("Source", "Target", "Weight")
edge_list_blue <- edge_list_blue[edge_list_blue$Weight > 0, ]

edge_list_blue <- edge_list_blue %>%
  filter(Weight > 0 & Weight < 1)

write.csv(edge_list_blue, "blue_edge_list.csv", row.names = FALSE)

#---------------------------------------------------------------------
#Code to get edge list of edges only related to Cre06.g273100 in the module

filtered_edge_list_blue <- edge_list_blue %>%
  filter(Source == "Cre06.g273100" | Target == "Cre06.g273100")

filtered_edge_list_blue <- filtered_edge_list_blue %>%
  arrange(desc(Weight))

write.csv(filtered_edge_list_blue, "filtered_blue_edge_list.csv", row.names = FALSE)
#---------------------------------------------------------------------
#Code for getting just top 20 nodes in Cre06.g273100 blue module

filtered_edge_list_blue_removed_duplicates <- filtered_edge_list_blue %>%
  mutate(Edge = paste(pmin(Source, Target), pmax(Source, Target), sep = "-")) %>%
  distinct(Edge, .keep_all = TRUE) %>%
  select(-Edge)

filtered_edge_list_blue_top_20 <- filtered_edge_list_blue_removed_duplicates %>%
  arrange(desc(Weight)) %>%
  slice(1:20)

filtered_edge_list_blue_top_20$Source <- as.character(filtered_edge_list_blue_top_20$Source)
filtered_edge_list_blue_top_20$Target <- as.character(filtered_edge_list_blue_top_20$Target)

write.csv(filtered_edge_list_blue_top_20, "filtered_edge_list_blue_top_20_nodes.csv", row.names = FALSE)

#---------------------------------------------------------------------
#Code to see if any polysaccahride metabolism related genes from Phycocosm
#are in the same module as the LPMO-containg proteins

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

genes_of_interest_modules <- gene_module_df %>%
  filter(gene_id %in% genes_of_interest)

print(genes_of_interest_modules)

write.csv(genes_of_interest_modules, "genes_of_interest_modules_137.csv", row.names = FALSE)

filtered_genes_of_interest_modules <- genes_of_interest_modules %>%
  filter(module %in% c("blue", "brown", "darkviolet"))

print(filtered_genes_of_interest_modules)

write.csv(filtered_genes_of_interest_modules, "filtered_genes_of_interest_modules.csv", row.names = FALSE)

#--------------------------------------------------------------------
