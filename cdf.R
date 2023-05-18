# Import libraries
library(rio)
library(tidyverse)

# Read in files
dir <- "/projectnb/bf528/users/group_3/Project_4/Data_Curator/"
cdf4 <- rio::import(paste0(dir, "SRR3879604_barcode_counts.txt"), col.names = c("barcodes", "counts"))
cdf5 <- rio::import(paste0(dir, "SRR3879605_barcode_counts.txt"), col.names = c("barcodes", "counts"))
cdf6 <- rio::import(paste0(dir, "SRR3879606_barcode_counts.txt"), col.names = c("barcodes", "counts"))

cdf4 <- cdf4 %>%
  dplyr::arrange(desc(counts)) %>%
  dplyr::mutate(ratio = counts / sum(counts)) %>%
  dplyr::group_by(counts) %>%
  dplyr::summarize(occurance = n(),
                   ratio = sum(ratio)) %>%
  dplyr::mutate(sum_ratio = cumsum(ratio),
                id = "SRR3879604")
cdf5 <- cdf5 %>%
  dplyr::arrange(desc(counts)) %>%
  dplyr::mutate(ratio = counts / sum(counts)) %>%
  dplyr::group_by(counts) %>%
  dplyr::summarize(occurance = n(),
                   ratio = sum(ratio)) %>%
  dplyr::mutate(sum_ratio = cumsum(ratio),
                id = "SRR3879605")
cdf6 <- cdf6 %>%
  dplyr::arrange(desc(counts)) %>%
  dplyr::mutate(ratio = counts / sum(counts)) %>%
  dplyr::group_by(counts) %>%
  dplyr::summarize(occurance = n(),
                   ratio = sum(ratio)) %>%
  dplyr::mutate(sum_ratio = cumsum(ratio),
                id = "SRR3879606")
cdf_combined <- cdf4 %>%
  dplyr::bind_rows(list(cdf5, cdf6))

# ggplot
plot4 <- ggplot(cdf4, aes(x = counts, y = sum_ratio)) +
  geom_point(color = "#66C2A5") +  
  geom_line(color = "#66C2A5") + 
  theme_minimal()
  labs(x = "Counts", y = "Cumulative Ratio", title = "Cumulative Distribution Function") 
  
plot5 <- ggplot(cdf5, aes(x = counts, y = sum_ratio)) +
  geom_point(color = "#FC8D62") +  
  geom_line(color = "#FC8D62") + 
  theme_minimal()
labs(x = "Counts", y = "Cumulative Ratio", title = "Cumulative Distribution Function") 

plot6 <- ggplot(cdf6, aes(x = counts, y = sum_ratio)) +
  geom_point(color = "#8DA0CB") +  
  geom_line(color = "#8DA0CB") + 
  theme_minimal()
labs(x = "Counts", y = "Cumulative Ratio", title = "Cumulative Distribution Function") 

plot_combined <- ggplot(cdf_combined, aes(x = counts, y = sum_ratio, color = id)) +
  geom_point() +  
  geom_line() + 
  scale_color_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB"), name = "Samples") +
  theme_minimal() + 
  geom_hline(yintercept = 0.05) +
  geom_vline(xintercept = 250) +
  labs(x = "Counts", y = "Cumulative Ratio", title = "Cumulative Distribution Function")

plot_combined_log <- ggplot(cdf_combined, aes(x = counts, y = sum_ratio, color = id)) +
  geom_point() +  
  geom_line() + 
  scale_x_log10() +
  scale_color_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB"), name = "Samples") +
  theme_minimal() + 
  geom_hline(yintercept = 0.05) +
  geom_vline(xintercept = 250) +
  labs(x = "Log10(Counts)", y = "Cumulative Ratio", title = "Cumulative Distribution Function") 


ggsave("plot4.png", plot = plot4, units="in", width=8, height=6, dpi = 600)
ggsave("plot5.png", plot = plot5, units="in", width=8, height=6, dpi = 600)
ggsave("plot6.png", plot = plot6, units="in", width=8, height=6, dpi = 600)
ggsave("plot_combined.png", plot = plot_combined, units="in", width=8, height=6, dpi = 600)
ggsave("plot_combined_log.png", plot = plot_combined_log, units="in", width=8, height=6, dpi = 600)

