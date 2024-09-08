
library(dplyr)

#---------------------------------------------------------------------
#Code to obtain FPKM values of sample of Chlamydomonas reinhardtii monoculture
raw_counts <- read.table("crcounts.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

head(raw_counts)
str(raw_counts)

total_counts <- sum(raw_counts$...aligned_reads.aligned_reads_sorted.bam)
raw_counts <- raw_counts %>%
  mutate(FPKM = ...aligned_reads.aligned_reads_sorted.bam / (Length / 1000) / (total_counts / 1e6))
head(raw_counts)

write.table(raw_counts, "fpkm_crcounts.txt", sep = "\t", quote = FALSE, row.names = FALSE)

filtered_fpkm <- raw_counts %>%
  filter(Geneid %in% c("Cre07.g317250.v5.5", "Cre06.g270500.v5.5", "Cre06.g273100.v5.5")) %>%
  dplyr::select(Geneid, FPKM) %>%
  mutate(Sample = "C.reinhardtii")


#---------------------------------------------------------------------
##Code to obtain FPKM values of sample of C. reinhardtii + S. cerevisiae coculture

raw_counts_coculture <- read.table("counts_DRR513084.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

head(raw_counts_coculture)
str(raw_counts_coculture)

total_counts_coculture <- sum(raw_counts_coculture$`X.users.pr801.aligned_reads.aligned_reads_DRR513084_sorted.bam`)

raw_counts_coculture <- raw_counts_coculture %>%
  mutate(FPKM = `X.users.pr801.aligned_reads.aligned_reads_DRR513084_sorted.bam` / (Length / 1000) / (total_counts_coculture / 1e6))

filtered_fpkm_coculture <- raw_counts_coculture %>%
  filter(Geneid %in% c("Cre07.g317250.v5.5", "Cre06.g270500.v5.5", "Cre06.g273100.v5.5")) %>%
  dplyr::select(Geneid, FPKM) %>%
  mutate(Sample = "C.reinhardtii+S.cerevisae")

#---------------------------------------------------------------------
#Code to make plots to look for significant differences between
#filtered_fpkm (C. reinhardtii alone) and filtered_fpkm_coculture
#(coculture of C. reinhardtii + S. cerevisae)

library(ggplot2)

combined_fpkm <- rbind(filtered_fpkm, filtered_fpkm_coculture)

combined_fpkm_plot_bars <- combined_fpkm %>%
  group_by(Geneid, Sample) %>%
  summarise(
    mean_FPKM = mean(FPKM),
    se_FPKM = sd(FPKM) / sqrt(n())
  )

print(combined_fpkm)

combined_fpkm$Geneid <- sub("\\.v5\\.5", "", combined_fpkm$Geneid)
combined_fpkm_plot_bars <- combined_fpkm_plot_bars %>%
  mutate(Geneid = sub("\\.v5\\.5", "", Geneid))

ggplot(combined_fpkm_plot_bars, aes(x = Geneid, y = mean_FPKM, fill = Sample)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean_FPKM - se_FPKM, ymax = mean_FPKM + se_FPKM),
                width = 0.2, position = position_dodge(0.7)) +
  theme_minimal() +
  labs(title = expression("FPKM Comparison between " * italic("C. reinhardtii") * " monoculture and Coculture with " * italic("S. cerevisiae")),
       x = "Gene ID",
       y = "Mean FPKM",
       fill = " Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
