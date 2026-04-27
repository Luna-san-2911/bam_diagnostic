# Load necessary libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gridExtra)
})

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Ensure 5 arguments are provided
if (length(args) != 5) {
  stop("Please provide exactly 5 arguments: 1) flagstat file, 2) mapq file, 3) mosdepth global dist file, 4) mosdepth regions bed, and 5) output prefix")
}

# Assign arguments to variables
flagstat_file <- args[1]
mapq_file <- args[2]
mosdepth_dist_file <- args[3]
mosdepth_regions_file <- args[4]
output_prefix <- args[5]

# 1. Mapped reads stats
data <- read.delim(flagstat_file, sep = " ", header = FALSE)
pmr1 <- substring(data[7,5], 2)
pmr2 <- substr(data[7,7], 1, nchar(data[7,7]) - 1)
data[6,4] <- 'primary_duplicates'
data[8,4] <- 'primary_mapped'
data <- data[c(2:8),c(1,3,4)]
data$V1 <- as.numeric(data$V1)
data$V3 <- as.numeric(data$V3)

data_long <- data %>%
  rename(QC_passed_reads = V1, QC_failed_reads = V3) %>%
  pivot_longer(
    cols = c(QC_passed_reads, QC_failed_reads),
    names_to = "type",
    values_to = "values")

annotation_text <- paste0("QC-passed mapped reads: ", pmr1, "\n",
                          "QC-failed mapped reads: ", pmr2)

p1 <- ggplot(data_long, aes(y = values, x = '', fill = type)) +
  geom_col(position = "dodge", color='black') +
  scale_fill_manual(values = c("#d6ad2aff", "#642e8cff"))+
  labs(x = '',
       y = 'Read count',
       fill = '',
       title='Read Alignment Categories',
       tag = annotation_text) +          
  theme_bw() +
  facet_wrap(~ V4, scales = "free_y") +
  theme(legend.position = 'top',
        axis.ticks.x = element_blank(),
        plot.tag.position = c(0.7, 0.2),
        plot.tag = element_text(size = 10, hjust = 0.5))

# 2. Histogram of quality
mapq_data <- read.delim(mapq_file, header = FALSE)
mapq_data$V1 <- as.numeric(mapq_data$V1)
mapq_data$V2 <- as.numeric(mapq_data$V2)
wm <- weighted.mean(mapq_data$V1, w = mapq_data$V2)

p2 <- ggplot(mapq_data, aes(x = V1, y = V2)) +
  geom_col(color = 'black', fill = '#b81d1bff') +
  labs(x = "MQ", y = "Reads",
       title = paste0('Mapping Quality (mean: ', round(wm, 2), ')')) +
  theme_bw()

# 3. Coverage per chromosome
data3 <- read.delim(mosdepth_dist_file, header=F)
data3 <- data3[data3$V2 == 1, ]
data3$V4 <- data3$V3*100
data3$V1 <- as.factor(data3$V1)

p3 <- ggplot(data3, aes(x = V1,y=V4)) +
  geom_col(color='black', fill='#58bb9dff') +
  labs(x = "", y = "Covered by at least 1 read (%)",
       title = 'Breadth of Coverage') +
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust = 1))

# 4. Proportion of genome at coverage / Depth chromosomes
data4 <- read.delim(mosdepth_regions_file, header=F)
data4$V1 <- as.factor(data4$V1)
data4$V2 <- as.numeric(data4$V2)

p4 <- ggplot(data4, aes(x = V2, y = V4)) +
  geom_area(fill = "darkblue") +
  facet_wrap(~ V1, ncol = 1, scales = "free", strip.position = "right") +
  scale_x_continuous(limits = c(0, NA), expand = c(0, 0)) +
  labs(x = "Genomic Position", y = "Depth") +
  theme_bw()

mean_depth_chr <- data4 %>%
  mutate(window_size = V3 - V2) %>%
  group_by(V1) %>%
  summarize(mean_depth = weighted.mean(V4, w = window_size)) %>%
  mutate(V1 = as.character(V1)) %>% 
  add_row(
    V1 = "total",
    mean_depth = weighted.mean(data4$V4, w = data4$V3 - data4$V2)
  )

p5 <- ggplot(mean_depth_chr, aes(x = V1,y=mean_depth)) +
  geom_col(color='black', fill='#b35623ff') +
  labs(x = "", y = "Mean depth",
       title = 'Depth of Aligned Reads') +
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust = 1))

# Save plots
plot1 <- arrangeGrob(p1,p2,p3,p5,ncol=2)

out_file_1 <- paste0(output_prefix, "_alignment_summary_1.png")
out_file_2 <- paste0(output_prefix, "_alignment_summary_2.png")

ggsave(plot1,
       filename=out_file_1,
       dpi=600,
       height = 10,
       width = 10)

ggsave(p4,
       filename=out_file_2,
       dpi=600,
       height = 10,
       width = 10)
