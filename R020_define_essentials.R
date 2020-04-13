## Defining the essential genes

library(tidyverse)

source('R001_common.R')


## FUNCTIONS ------------------------------------------------------------------
get_limd_peak_m_s <- function(lvals, breaks = bin_breaks, rettype = NULL) {
  h <- hist(lvals, breaks = breaks, plot = FALSE)
  counts2 <- h$counts; counts2[1] <- 0; counts2[length(counts2)] <- 0
  highest_bin <- match(max(counts2), counts2)
  high_bins <- rank(desc(counts2)) <= 0.02 * length(counts2)
  avg_high_bin_val <- mean(counts2[high_bins])
  peak_calc_bins <- counts2 > 0.1 * avg_high_bin_val  # 0.1 is somewhat arbitrary; could try 0.05 or other values
  left_peak_bin <- max(seq(1, highest_bin)[!peak_calc_bins[1:highest_bin]]) + 1  # the last T left of the peak before an F
  right_peak_bin <- min(seq(highest_bin, length(counts2))[
    !peak_calc_bins[highest_bin:length(counts2)]]) - 1     # the last T right of the peak before an F
  peak_bounds <- c( h$breaks[left_peak_bin], h$breaks[right_peak_bin + 1])
  peak_bin_set <- h$mids > peak_bounds[1] & h$mids < peak_bounds[2]
  ms_calc_set <- lvals[!is.na(lvals) & lvals >= peak_bounds[1] & lvals <= peak_bounds[2]]
  wt = sum(h$counts[peak_bin_set]) / sum(peak_bin_set)
  m <- mean(ms_calc_set)
  s <- sd(ms_calc_set)
  ms <- list(mean = m, sd = s, wt = wt, bound_lo = peak_bounds[1], bound_hi = peak_bounds[2])
  if (is.null(rettype)) return(ms)
  if (rettype == "m") return(m)
  if (rettype == "s") return(s)
  if (rettype == "w") return(wt)
  if (rettype == "l") return(peak_bounds[1])
  if (rettype == "h") return(peak_bounds[2])
}

## LOAD DATA -------------------------------------------------------------------

TBGpkNm <- read_tsv(files$TBGpkNm)


## Define Essentials -----------------------------------------------------------

bin_breaks <- seq(from = (-3), to = 1, length.out = 201)

# make table with columns:
#  - lrpkNm (log of size-factor-normalized reads per gene per gene length)
#  - lrpkNmNp (difference in lrpkNm values from MRg1a sample for each gene)
#  - lim_lNN (as lrpkNmNp, but values constrained to within bin_breaks limits)
dat <- TBGpkNm %>%
  left_join(TBGpkNm %>%
              filter(as.character(Sample) == 'MRg1a') %>%
              select(Locus, t0rpkNm = rpkNm),
            by = "Locus") %>%
  mutate(lrpkNm = log2(rpkNm), 
         lt0rpkNm = log2(t0rpkNm),
         lrpkNmNp = lrpkNm - lt0rpkNm) %>%
  select(-t0rpkNm, -lt0rpkNm, -Portion) %>%
  mutate(lim_lNN = ifelse(lrpkNmNp < first(bin_breaks), first(bin_breaks),
                          ifelse(lrpkNmNp > last(bin_breaks), last(bin_breaks),
                                 lrpkNmNp)))

# get statistical and other data per sample
msdat <- dat %>%
  group_by(Sample) %>%
  summarise(n = n(),
            m = get_limd_peak_m_s(lim_lNN)$mean,
            s = get_limd_peak_m_s(lim_lNN)$sd,
            w = get_limd_peak_m_s(lim_lNN)$wt,
            bl = get_limd_peak_m_s(lim_lNN)$bound_lo,
            bh = get_limd_peak_m_s(lim_lNN)$bound_hi) %>%
  mutate(qn3 = qnorm(0.001, m, s),
         qn5 = qnorm(0.00001, m, s),
         pn_min1 = pnorm(-1, m, s))

# Define "essentials" for each sample (i.e., depleted below 1e-5 p-value for fitted normal curve)
essl_lists <- dat %>%
  filter(!(Sample %in% c('MRg1a'))) %>%
  left_join(select(msdat, Sample, qn5), by = "Sample") %>%
  mutate(essEm5 = lrpkNmNp < qn5) %>%
  select(Locus, Sample, essEm5) %>%
  left_join(select(annot, Locus, wbp, Gene, Product), by = "Locus") %>%
  spread(key = Sample, value = essEm5) %>%
  mutate(B_2h = T1_2h & T2_2h,
         B_4h = T1_4h & T2_4h,
         B_6h = T1_6h & T2_6h,
         B_8h = T1_8h & T2_8h,
         B_out = AyL3 & AyL4)  # e.g., TRUE in Both of the AyL samples

# Define the final list of essentials
exclude_loci <- c('ACIAD0560', 'ACIAD1220', 'ACIAD1310', 'ACIAD1322', 'ACIAD1371',
                  'ACIAD2029', 'ACIAD2432', 'ACIAD3063', 'ACIAD3107', 'ACIAD3368',
                  'ACIAD3570')  # genes excluded because of hit bias implying polarity of anti insertions

essentials_final <- essl_lists %>%
  select(Locus, wbp, Gene, Product, B_8h, B_out) %>%
  mutate(Ess = B_8h & B_out) %>%
  filter(Ess) %>%
  select(-B_8h, -B_out, -Ess) %>%
  filter(!(Locus %in% exclude_loci))

## Write output files
# Final list of genes defined as essential
write_tsv(essentials_final, files$essentials_final)
# Table of 'essentiality' (e.g., significantly depleted) per sample
write_tsv(essl_lists, files$essentiality_by_sample)

