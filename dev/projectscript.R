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

#Anyways next step of thermotolerance in wheat 