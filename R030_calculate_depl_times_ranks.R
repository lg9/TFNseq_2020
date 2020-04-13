## Calculate depletion times and ranks

library(tidyverse)

source('R001_common.R')


## Parameters ------------------------------------------------------------------

depl_thresh <- (-1.0)  # threshold delta-log-reads value defining depletion time


## FUNCTIONS -------------------------------------------------------------------
locus_depl_time <- function(times, log_reads, thresh = depl_thresh) {
  # get the number of the last lrpk value above the threshold
  aboves <- which(log_reads >= thresh)
  last_above <- last(aboves)
  # return NaN if it's the last time point
  if (last_above == length(log_reads)) { return( NaN ) }
  # calculate the depletion time
  dt <- times[last_above] + ( times[last_above + 1] - times[last_above] ) *
    ( thresh - log_reads[last_above] ) / ( log_reads[last_above + 1] - log_reads[last_above] )
  return(dt)
}

## GET DATA -------------------------------------------------------------------

TBGpkNmNp <- read_tsv(files$TBGpkNmNp)

essentials_final <- read_tsv(files$essentials_final)


## CALCULATE DEPL TIMES AND RANKS ---------------------------------------------

## Make data table showing times and Log-transform rpkNN values
depl_dat <- TBGpkNmNp %>%
  mutate(
    lrpk = log2(rpkNmNpre)
  ) %>%
  select(
    -rpkNmNpre
  ) %>%
  left_join(
    outgr_times, by = "Sample"
  ) %>%
  select(-Portion, -hits, -lrpk, lrpk)


## Trial_1 depletion times all genes
depl_times_T1 <- depl_dat %>%
  filter(
    Sample %in% c('MRg1a', 'T1_2h', 'T1_4h', 'T1_6h', 'T1_8h', 'AyL4')
  ) %>%
  group_by(
    Locus
  ) %>%
  summarize(
    dt_T1 = locus_depl_time(time, lrpk, depl_thresh)
  ) %>%
  filter(
    !is.na(dt_T1)
  )

## Trial_2 depletion times all genes
depl_times_T2 <- depl_dat %>%
  filter(
    Sample %in% c('MRg1a', 'T2_2h', 'T2_4h', 'T2_6h', 'T2_8h', 'AyL3')
  ) %>%
  group_by(
    Locus
  ) %>%
  summarize(
    dt_T2 = locus_depl_time(time, lrpk, depl_thresh)
  ) %>%
  filter(
    !is.na(dt_T2)
  )

## Rank all genes per trial (independent of final essentiality call)
d_t_r_all <- annot %>%
  select(
    Locus
  ) %>%
  left_join(
    depl_times_T1, by = "Locus"
  ) %>%
  mutate(
    rank_T1_all = min_rank(dt_T1)
  ) %>%
  left_join(
    depl_times_T2, by = "Locus"
  ) %>%
  mutate(
    rank_T2_all = min_rank(dt_T2)
  )

# Write output
write_tsv(d_t_r_all, files$depl_times_ranks_all)

## Rank only final essentials
d_t_r_essentials <- 
  full_join(
    depl_times_T1, depl_times_T2,
    by = "Locus"
  ) %>%
  filter(
    Locus %in% essentials_final$Locus
  ) %>%
  mutate(
    rank_T1 = min_rank(dt_T1),
    rank_T2 = min_rank(dt_T2)
  ) %>%
  select(
    Locus, dt_T1, rank_T1, dt_T2, rank_T2
  )

# Write output
write_tsv(d_t_r_essentials, files$depl_times_ranks_essls)
