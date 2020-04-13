## For making processes boxplot and associated supplementary table


library(tidyverse)

source('R001_common.R')


## GET SOURCE DATA ------------------------------------------------------------

subprocess_loci <- read_tsv(files$subproc_loci_list)
d_t_r_essentials <- read_tsv(files$depl_times_ranks_essls)
subproc_order_inc <- read_tsv(files$subproc_order_incl)


## MAKE SUMMARY TABLE ---------------------------------------------------------

# Make grouped data table
dat <- subprocess_loci %>%
  left_join(
    select(
      annot,
      Locus, Gene, Label
    ),
    by = "Locus"
  ) %>%
  left_join(
    select(
      d_t_r_essentials,
      Locus, dt_T1                # alternatively, try with dt_T2
    ),
    by = "Locus"
  ) %>%
  rename(
    dt = starts_with('dt')
  ) %>%
  select(
    -Locus, -Gene
  )

# Get lists of the genes by SubProcess and Class (group 'C', 'D', etc together)
gene_lists_lower_classes <- dat %>%
  filter(
    !(Class %in% c('A', 'B'))
  ) %>%
  group_by(
    SubProcess
  ) %>%
  summarize(
    Poss_related = paste(sort(Label), collapse = ", ")
  )

gene_counts_class_A <- dat %>%
  filter(
    Class == 'A'
  ) %>%
  group_by(
    SubProcess
  ) %>%
  summarize(
    nA = n()
  )

gene_lists_by_group <- dat %>%
  filter(
    Class %in% c('A', 'B')
  ) %>%
  group_by(
    SubProcess, Class
  ) %>%
  summarize(
    genes = paste(sort(Label), collapse = ", ")
  ) %>%
  spread(
    key = Class, value = genes
  ) %>%
  rename(
    Central_genes = A, Other_genes = B
  ) %>%
  full_join(
    gene_lists_lower_classes,
    by = "SubProcess"
  )

noness_B_gene_lists <- dat %>%
  left_join(
    select(
      annot,
      Label, Locus
    ),
    by = "Label"
  ) %>%
  filter(
    Class == 'B',
    !(Locus %in% d_t_r_essentials$Locus)
  ) %>%
  select(
    SubProcess, Label
  ) %>%
  group_by(
    SubProcess
  ) %>%
  summarize(
    noness_B = paste(sort(unique(Label)), collapse = ", ")
  )

# Get the 'A' group depletion summary statistics by SubProcess
depl_stats_by_group <- dat %>%
  ungroup() %>%
  filter(
    Class == 'A'
  ) %>%
  select(
    SubProcess, dt
  ) %>%
  group_by(
    SubProcess
  ) %>%
  summarize(
    dt_median = median(dt, na.rm = TRUE),
    dt_25q = quantile(dt, probs = 0.25, na.rm = TRUE),
    dt_75q = quantile(dt, probs = 0.75, na.rm = TRUE),
    dt_min = min(dt, na.rm = TRUE),
    dt_max = max(dt, na.rm = TRUE)
  )
  
# Compile summary table data
summary_table <- gene_lists_by_group %>%
  left_join(
    noness_B_gene_lists,
    by = "SubProcess"
  ) %>%
  full_join(
    depl_stats_by_group,
    by = "SubProcess"
  ) %>%
  full_join(
    gene_counts_class_A,
    by = "SubProcess"
  ) %>%
  full_join(
    select(
      subproc_order_inc,
      SubProcess, Process, Order_Incl
    ),
    by = "SubProcess"
  ) %>%                        # to sort by median depl time  
  arrange(
    dt_median
  )

# Write output
write_tsv(summary_table, files$processes_summary_table)


## BOX PLOT -------------------------------------------------------------------

# option to order by median or specified order (see commented lines below)

mk_new_lab <- function(lab, n) {
  newlab <- paste(lab, " (", as.character(n), ")", sep = "")
  newlab
}

gene_counts_class_A$new_Label <- 
  mapply(mk_new_lab, gene_counts_class_A$SubProcess, gene_counts_class_A$nA)

plot_dat <- dat %>%
  ungroup() %>%
  left_join(
    select(
      summary_table,
      SubProcess, dt_median
    ),
    by = "SubProcess"
  ) %>%
  left_join(
    select(
      subproc_order_inc,
      SubProcess, Process, Order_Incl, Plot, Color
    ),
    by = "SubProcess"
  ) %>%
  filter(
    Class == 'A',
    !is.na(Order_Incl)           # use this line for full plot
    #!is.na(Order_Incl), Plot    # use this line for abridged plot
  ) %>%
  mutate(
    Category = factor(Process)
  ) %>%
  select(
    SubProcess, Category, dt, dt_median, Order_Incl, Color
  ) %>%
  left_join(
    select(
      gene_counts_class_A,
      SubProcess, new_Label
    ),
    by = "SubProcess"
  ) %>%
  mutate(
    Process = fct_reorder(new_Label, dt_median, .desc = TRUE)   # to order by median value
    #Process = fct_reorder(SubProcess, Order_Incl, .desc = TRUE)  # to order by specified order
  ) %>%
  select(
    -new_Label
  )

# set up key order and inclusion (manually adjust here as needed)
key_list = levels(plot_dat$Category)[c(5, 7, 9, 3, 2, 4, 10, 6, 8)]
# or set it up automatically
cat_min_meds <- plot_dat %>%
  ungroup() %>%
  select(Category, dt_median) %>%
  group_by(Category) %>%
  summarize(min_dt_med = min(dt_median)) %>%
  filter(Category != levels(plot_dat$Category)[1]) %>%
  arrange(min_dt_med)
keys_by_first_med <- cat_min_meds$Category

# set up color palette
col_palette <- c("#F0F0F0",
                 "#D55E00", "#CC79A7", "#009E73", "#56B4E9", "#945ED2",
                 "#7E6148", "#0072B2", "#E69F00", "#F0E542")
names(col_palette) <- c(levels(plot_dat$Category)[1], as.character(cat_min_meds$Category))

## Make the plot
set.seed(95)
p <- plot_dat %>%
  ggplot(aes(x = Process, y = dt)) +
  geom_point(aes(fill = Category), shape = NA, color = NA, show.legend = TRUE) +
  geom_boxplot(
    aes(fill = Category), show.legend = FALSE, outlier.shape = NA
    # outlier.shape = 21, outlier.color = NA, outlier.fill = "Black"
  ) +
  geom_jitter(width = 0.2, 
              shape = 21, 
              size = 1.25, 
              fill = "#FFFFFF00", 
              color = "#35353565") +
  scale_fill_manual(
    values = col_palette, name = NULL,
    breaks = keys_by_first_med
  ) +
  coord_flip() +
  labs(x = "", y = "Depletion time (hours)") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    legend.position = c(0.84, 0.84),
    legend.key.height = unit(0.03, "npc"),
    plot.margin = margin(1, 4, 1, 0),
    plot.background = element_blank()
    #axis.ticks.y = element_blank()
    #text = element_text(family = "Arial")
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = 22, size = 4
      )
    )
  )


save_file = '070_plots/boxplot.pdf'
pdf(file = save_file, width = 186 / 25.4, useDingbats = FALSE)
p
dev.off()




## version without individual points
p2 <- plot_dat %>%
  ggplot(aes(x = Process, y = dt)) +
  geom_point(aes(fill = Category), shape = NA, color = NA, show.legend = TRUE) +
  geom_boxplot(
    aes(fill = Category), show.legend = FALSE
    # outlier.shape = 21, outlier.color = NA, outlier.fill = "Black"
  ) +
  scale_fill_manual(
    values = col_palette, name = NULL,
    breaks = keys_by_first_med
  ) +
  coord_flip() +
  labs(x = "", y = "Depletion time (hours)") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    legend.position = c(0.84, 0.84),
    legend.key.height = unit(0.03, "npc"),
    plot.margin = margin(1, 4, 1, 0),
    plot.background = element_blank()
    #axis.ticks.y = element_blank()
    #text = element_text(family = "Arial")
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = 22, size = 4
      )
    )
  )









