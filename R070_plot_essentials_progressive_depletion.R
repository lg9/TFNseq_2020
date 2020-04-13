
library(tidyverse)
library(colorspace)

source('R001_common.R')

## GET DATA -------------------------------------------------------------------

TBGpkNm <- read_tsv(files$TBGpkNm)
essentials_final <- read_tsv(files$essentials_final)

## PARAMS ---------------------------------------------------------------------

levels_order <- c('T1_2h', 'T1_4h', 'T1_6h', 'T1_8h', 'AyL4')
minrep <- 1
sample_labels <- c(
  'T1_2h' = '2h',
  'T1_4h' = '4h',
  'T1_6h' = '6h',
  'T1_8h' = '8h',
  'AyL4' = '18h'
)

## Make table of data for plotting --------------------------------------------
pdat <- TBGpkNm %>%
  select(
    Locus, Sample, rpkNm
  ) %>%
  group_by(
    Locus
  ) %>%
  mutate(
    MRg_rpkNm = first(rpkNm)
  ) %>%
  filter(
    Sample %in% levels_order
  ) %>%
  mutate(
    essl = ifelse(Locus %in% essentials_final$Locus, "essential gene", "non-essential gene"),
    Sample = factor(Sample, levels = levels_order),
    limd_rpk = ifelse(rpkNm < minrep, minrep, rpkNm)
  ) %>%
  select(
    Locus, essl, Sample, MRg_rpkNm, rpkNm, limd_rpk
  )

pdat$limd_rpk[pdat$limd_rpk == minrep] <- jitter(pdat$limd_rpk[pdat$limd_rpk == minrep], factor = 2.8)

essl_colors <- c('#CC79A79F', '#30303055')
names(essl_colors) <- levels(factor(pdat$essl))

# Make the plot
p <- ggplot(pdat) +
  geom_point(
    aes(MRg_rpkNm, limd_rpk, color = essl),
    shape = 21, size = 2,
    fill = NA
    #color = "#FFFFFF00",
  ) +
  scale_color_manual(
    values = essl_colors,
    name = NULL
  ) +
  scale_x_log10() + 
  scale_y_log10(
    limits = c(min(pdat$limd_rpk), max(pdat$rpkNm)),
    expand = expand_scale(mult = 0.025),
    breaks = c(1, 10, 100, 1000, 10000),
    labels = c("<= 1", "10", "100", "1000", "10000")
  ) +
  facet_grid(
    cols = vars(Sample),
    labeller = as_labeller(sample_labels)
  ) +
  theme_classic() +
  theme(
    strip.background = element_rect(color = "#FFFFFF00", fill = "#C7C7C7"),
    axis.title = element_text(size = 10),
    legend.position = c(0.2, 0.14),
    legend.text = element_text(size = 11)
  ) +
  labs(
    x = 'Reads per gene kb (mutagenized DNA)',
    y = 'Reads per gene kb - after transformation and growth'
  )

p

save_file = '070_plots/progressive_depletion.pdf'
pdf(file = save_file, width = 183 / 25.4, height = 112 / 25.4)
p
dev.off()
