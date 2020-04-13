
library(tidyverse)
library(ggrepel)

source('R001_common.R')


## GET SOURCE DATA ------------------------------------------------------------

microcol_areas_counts <- read_tsv(files$microcol_areas_counts)

## Make Table -----------------------------------------------------------------

loci_order <- c(
  'nodna' = 'wt',
  'atpB-C' = 'atpB-C',
  'nrdAB' = 'nrdAB',
  'rpsA' = 'rpsA',
  'dnaE' = 'dnaE',
  'rpoA' = 'rpoA',
  'lpxB' = 'lpxB',
  'bamA' = 'bamA',
  'ftsZ' = 'ftsZ'
)

# manually entering depletion times summary:
#  - fastest dep time of the genes deleted, and
#  - median dep time for the genes of the category
depl_times_summary <- tibble(
  Gene = loci_order,
  deptm_fastest = c(NaN, 1.84, 2.30, 3.32, 3.37, 2.68, 4.47, 4.75, 5.44),
) %>% mutate(Gene = factor(Gene, levels = loci_order))

microcol_areas_counts <- microcol_areas_counts %>%
  mutate(
    Gene = factor(Gene, levels = loci_order),
    nCells = as.integer(nCells)
  )

microcols_summary <- microcol_areas_counts %>%
  select(
    -path_ID
  ) %>%
  group_by(
    Gene
  ) %>%
  summarize(
    n = n(),
    Min_A = min(Area),
    Max_A = max(Area),
    Mean_A = mean(Area),
    SD_A = sd(Area),
    SE_A = sd(Area) / sqrt(n()),
    CI95_A = SE_A * qt(0.95/2 + 0.5, n()-1),  # 95% confidence interval
    Median_A = median(Area),
    Q25_A = quantile(Area, probs = 0.25, na.rm = TRUE),
    Q75_A = quantile(Area, probs = 0.75, na.rm = TRUE),
    
    Min_N = min(nCells),
    Max_N = max(nCells),
    Mean_N = mean(nCells),
    SD_N = sd(nCells),
    SE_N = sd(nCells) / sqrt(n()),
    CI95_N = SE_N * qt(0.95/2 + 0.5, n()-1),  # 95% confidence interval
    Median_N = median(nCells),
    Q25_N = quantile(nCells, probs = 0.25, na.rm = TRUE),
    Q75_N = quantile(nCells, probs = 0.75, na.rm = TRUE)
  )  %>% 
  left_join(
    depl_times_summary, by = "Gene"
  )



## Simple barplot -----------------------------------------------------------

# p <- microcol_areas_counts %>%
#   select(
#     -path_ID
#   ) %>%
#   ggplot(
#     aes(Gene, Area)
#   ) +
#   geom_boxplot(
#     outlier.shape = NA
#   ) +
#   geom_jitter(
#     width = 0.25,
#     height = 0,
#     shape = 21,
#     size = 2,    # 3 for full-size
#     color = '#FFFFFF00',
#     fill = '#00000077'
#   ) +
#   theme_classic() +
#   labs(
#     x = "Locus",
#     y = "Microcolony area (um2)"
#   )

# save_file = '070_plots/colony_areas_barplot.pdf'
# pdf(file = save_file, width = 176 / 25.4, height = 150 / 25.4)
# p
# dev.off()
  


## Plot median microcolony area vs depletion time (fastest) -------------------
p1 <- microcols_summary %>%
  ggplot(
    aes(deptm_fastest, Median_A)
  ) +
  geom_point(
    # size = 1
  ) +
  geom_errorbar(
    aes(ymin = Q25_A, ymax = Q75_A)
  ) +
  # use geom_text_repel to ID points, but add text to final figure manually
  # geom_text_repel(
  #   aes(label = Gene), point.padding = 0.25
  # ) +
  labs(
    x = "Mutant depletion time (hours)",
    y = "Microcolony area (um2)"
  ) +
  theme_classic()

# Save the plot
pdf(file = '070_plots/microcol_area_v_depl_time.pdf', 
    width = 90 / 25.4, height = 90 / 25.4, useDingbats = FALSE)
p1
dev.off()

cor(microcols_summary$Median_A, microcols_summary$deptm_fastest, use = "complete.obs")
# 0.7225604


## Plot cell count vs microcol area -------------------------------------------
p2 <- microcols_summary %>%
  ggplot(
    aes(Median_A, Median_N)
  ) + 
  geom_point() +
  # use geom_text_repel to ID points, but add text manually to final figure using Illustrator
  # geom_text_repel(
  #   aes(label = Gene)
  # ) +
  labs(
    x = "Microcolony area (um2)",
    y = "Microcolony cell count"
  ) +
  theme_classic()

# Save the plot
pdf(file = '070_plots/cell_count_v_microcol_area.pdf', 
    width = 90 / 25.4, height = 90 / 25.4, useDingbats = FALSE)
p2
dev.off()
