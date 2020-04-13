## Spread and compile all the data for all-genes summary table
##   and for essentials summary table (Dataset S1)


library(tidyverse)

source('R001_common.R')


## GET SOURCE DATA ------------------------------------------------------------

subprocess_loci <- read_tsv(files$subproc_loci_list)
TBGpkNm <- read_tsv(files$TBGpkNm)
essl_lists <- read_tsv(files$essentiality_by_sample)
essentials_final <- read_tsv(files$essentials_final)
d_t_r_essentials <- read_tsv(files$depl_times_ranks_essls)


## SPREAD HIT AND READ DATA BY SAMPLE -----------------------------------------

hits_spread <- TBGpkNm %>%
  select(
    -Portion, -rpkNm
  ) %>%
  spread(
    key = Sample, value = hits
  ) %>%
  mutate_at(2:12, as.integer)

rpk_spread <- TBGpkNm %>%
  select(
    -Portion, -hits
  ) %>%
  spread(
    key = Sample, value = rpkNm
  )

# Make the table of spread hits and reads
h_r_spread <- 
  full_join(
    hits_spread, rpk_spread,
    by = "Locus",
    suffix = c(".h", ".r")
  ) %>%
  select(1, rep(1:11, each = 2) + rep(c(1, 12), 11)) %>%     # interspace the hits and reads
  select(1, 6:15, 4, 5, 16:23, 2, 3)                         # manually reorder

# Make spread table of delta log read values
delta_l2r_spread <- TBGpkNm %>%
  select(
    -Portion, -hits
  ) %>%
  left_join(
      TBGpkNm %>%
        filter(
          Sample == "MRg1a"
        ) %>%
        select(
          Locus, rpkNm
        ) %>%
        rename(
          ref_rpkNm = rpkNm
        ),
      by = "Locus"
    ) %>%
  mutate(
    delta_lr = log(rpkNm, 2) - log(ref_rpkNm, 2)
  ) %>%
  select(
    Locus, Sample, delta_lr
  ) %>%
  spread(
    key = Sample, value = delta_lr
  ) %>%
  select(
    Locus, T1_2h, T1_4h, T1_6h, T1_8h, AyL4,
    T2_2h, T2_4h, T2_6h, T2_8h, AyL3
  )

rm(TBGpkNm, hits_spread, rpk_spread)

# Make table listing functional classes
loci_functions <- annot %>%
  select(
    Locus
  ) %>%
  left_join(
    group_by(
      subprocess_loci,
      Locus
    ) %>%
      summarize(
        SubProcesses = paste(sort(SubProcess), collapse = ", ")
    ),
    by = "Locus"
  )
loci_functions$SubProcesses[is.na(loci_functions$SubProcesses)] <- ""  # change NA to blank

## MAKE GeneTable DATA TABLE (spread combined data for all genes) -------------

spread_gene_data_all <- annot %>%
  select(
    -Replicon, -Label
  ) %>%
  left_join(
    loci_functions,
    by = "Locus"
  ) %>%
  left_join(
    select(
      h_r_spread,
      Locus, MRg1a.h, MRg1a.r,
      T1_2h.h, T1_2h.r, T1_4h.h, T1_4h.r, T1_6h.h, T1_6h.r, T1_8h.h, T1_8h.r, AyL4.h, AyL4.r,
      T2_2h.h, T2_2h.r, T2_4h.h, T2_4h.r, T2_6h.h, T2_6h.r, T2_8h.h, T2_8h.r, AyL3.h, AyL3.r
    ),
    by = "Locus"
  ) %>%
  left_join(
    delta_l2r_spread,
    by = "Locus"
  ) %>%
  left_join(
    select(
      essl_lists,
      Locus,
      T1_2h, T1_4h, T1_6h, T1_8h, AyL4,
      T2_2h, T2_4h, T2_6h, T2_8h, AyL3
    ),
    by = "Locus",
    suffix = c(".dlr", ".ess")
  ) %>%
  left_join(
    select(
      d_t_r_essentials,
      Locus,
      dt_T1, rank_T1, dt_T2, rank_T2
    ),
    by = "Locus"
  )

# Write output
write_tsv(spread_gene_data_all, files$spread_gene_data_all)


## MAKE ESSENTIALS SPREAD DATA TABLE (Dataset S1) -------

spread_gene_data_ess <- spread_gene_data_all %>%
  select(
    -bp, -wbp,
    -T1_2h.ess, -T1_4h.ess, -T1_6h.ess, -T1_8h.ess, -AyL4.ess,
    -T2_2h.ess, -T2_4h.ess, -T2_6h.ess, -T2_8h.ess, -AyL3.ess
  ) %>%
  filter(
    Locus %in% essentials_final$Locus
  )

# Write output
write_tsv(spread_gene_data_ess, files$spread_gene_data_ess)







