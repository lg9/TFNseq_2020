## Validation plots
## - Trial2 v Trial1 depletion times
## - qPCR depletion rates v depletion times for selected deletions
##    (depletion rate = least-sqared-fit slopes of deletion qPCR starting quantities across time points) 

library(tidyverse)
library(ggrepel)

source('R001_common.R')


## GET DATA -------------------------------------------------------------------

essentials_final <- read_tsv(files$essentials_final)
d_t_r_all <- read_tsv(files$depl_times_ranks_all)
qPCR_SQ_table <- read_csv(files$qPCR_SQ_table)


## Plot T2 v T1 depl times for essentials  ------------------------------------

cor_tbl <- essentials_final %>%
  left_join(
    d_t_r_all,
    by = "Locus"
  ) %>%
  select(
    Locus, wbp, Gene, Product, dt_T1, dt_T2
  )

# Make the plot
p1 <- cor_tbl %>%
  mutate(
    wbp_cat = wbp < 500
  ) %>%
  ggplot(
    aes(dt_T1, dt_T2, color = wbp_cat)
  ) +
  geom_point(size = 1) +
  scale_color_manual(values = c("Black", "Black")) +
  # scale_color_manual(values = c("Black", "Red2")) +   # to demarcate small genes
  theme_classic() +
  labs(
    x = "Primary trial depletion time (hours)",
    y = "Secondary trial depletion time (hours)"
  ) +
  theme(legend.position = "none")

# Save plot
pdf(file = '070_plots/T2_v_T1_depl_times.pdf', 
    width = 90 / 25.4, height = 90 / 25.4, useDingbats = FALSE)
p1
dev.off()


# Calculate correlation
mdl <- lm(cor_tbl$dt_T2 ~ cor_tbl$dt_T1)

cor(cor_tbl$dt_T1, cor_tbl$dt_T2, use = "complete")
# 0.8324942
cor(cor_tbl$dt_T1[cor_tbl$wbp > 500], cor_tbl$dt_T2[cor_tbl$wbp > 500], use = "complete")
# 0.893196



## Plot deletion qPCR depletion rate vs TFNseq depletion time -----------------

# choose runs to label (manual here)
runs_to_label <- c('d102 190429c', 'd133 190507c', 'd120 190507b',
                   'd111 190430a', 'd112 190515b', 'd110 190513a',
                   'd103 190429a', 'd126 190514a',                  # 'd125 190513c', 
                   'd101 190515a', 'd106 190522a', 'd109 190717a',
                   'd107 190522b', 'd127 190514b', 'd132 190510a',
                   'd116 190430a', 'd116 190722a',
                   'd117 190722b')

# choose runs to exclude from graph (manual here)
runs_to_exclude <- c('d116 190430a', 'd117 190722b')  # removing duplicate qPCRs of atp locus

## Tidy qPCR data

qPCR_time_points <- tibble(
  tp = c("t_inoc", "t_2h", "t_3h", "t_4h", "t_5.5h", "t_7h", "t_9h"),
  hrs = c(0.0, 2.18, 3.10, 4.07, 5.55, 7.03, 9.03)
)

qPCR_tidy <- qPCR_SQ_table %>%
  pivot_longer(
    cols = starts_with("t_"),
    names_to = "tp",
    values_to = "delta_log2_SQ"
  ) %>%
  left_join(
    qPCR_time_points,
    by = "tp"
  ) %>%
  select(
    Og_run, loci, genes, rep_loci, hrs, delta_log2_SQ
  )

dt_v_slopes <- qPCR_tidy %>%
  group_by(
    Og_run, loci, genes, rep_loci
  ) %>%
  do(
    mod = lm(delta_log2_SQ ~ hrs, data = .)
  ) %>%
  do(
    data.frame(
      Og_run = .$Og_run,
      loci = .$loci,
      genes = .$genes,
      rep_locus = .$rep_loci,
      lm_slope = coef(.$mod)[2]
    )
  ) %>%
  left_join(
    d_t_r_all,
    by = c("rep_locus" = "Locus")
  ) %>%
  select(
    Og_run, loci, genes, rep_locus, dt_T1, lm_slope
  )


# Make the plot
p2 <- dt_v_slopes %>%
  mutate(
    Label = ifelse(Og_run %in% runs_to_label, genes, "")
  ) %>%
  filter(
    !(Og_run %in% runs_to_exclude),
    lm_slope < (-0.0001)
  ) %>%
  ggplot(
    aes(dt_T1, lm_slope),  # add aes(... color = genes) to distinguish points
  ) +
  geom_point() +
  theme_classic() +
  labs(
    x = "Tn-seq depletion time (hours)",
    y = "qPCR depletion rate (log2 per hour)"
  ) +
  # geom_text_repel(
  #   aes(label = Label)      # change to aes(label = genes) to show all
  # ) +
  theme(legend.position = "none")


pdf(file = '070_plots/qPCR_depl_rate_v_TFNseq_depl_time.pdf',
    width = 90 / 25.4, height = 90 / 25.4, useDingbats = FALSE)
p2
dev.off()

# Calculate correlation
sub_dt_v_slopes <- dt_v_slopes %>%
  mutate(
    Label = ifelse(Og_run %in% runs_to_label, genes, "")
  ) %>%
  filter(
    !(Og_run %in% runs_to_exclude),
    lm_slope < (-0.0001)
  )

cor(sub_dt_v_slopes$lm_slope, sub_dt_v_slopes$dt_T1, use = "complete")
# 0.6944288


