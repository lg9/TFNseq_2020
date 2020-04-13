## Normalize the reads per gene (per gene kb) by size factors
## (median rpkg relative to a pseudo-reference)


library(tidyverse)

source('R001_common.R')


## FUNCTIONS -------------------------------------------------------------------

geo.mean <- function(x, na.rm=TRUE, zeros.method=1) {
  # zeros methods:
  #     1 - add 1 to each value and subtract 1 from the result
  #     2 - ignore zeros
  #     3 - convert zeros to 1 for the calculation
  #     else - do nothing
  if(na.rm) { x <- x[!is.na(x)] }
  n <- length(x)
  if (zeros.method == 1) {
    x <- x + 1
    (prod(x)^(1/n) - 1)
  } else if (zeros.method == 2) {
    x <- x[x > 0]
    n <- length(x)
    prod(x)^(1/n)
  } else if (zeros.method == 3) {
    x[x==0] <- 1
    prod(x)^(1/n)
  } else {
    prod(x)^(1/n)
  }
}

make_rpk_column <- function(tbg = TBG, annot_table = annot) {
  tbg_rpk <- tbg %>%
    left_join(annot_table %>% select(Locus, bp, wbp), by = "Locus")
  tbg_rpk$portion_bp <- 0
  tbg_rpk$portion_bp[ tbg_rpk$Portion == 'ORF' ] <- tbg_rpk$bp[ tbg_rpk$Portion == 'ORF' ]
  tbg_rpk$portion_bp[ tbg_rpk$Portion == 'window' ] <- tbg_rpk$wbp[ tbg_rpk$Portion == 'window' ]
  tbg_rpk$rpk <- tbg_rpk$reads * 1000 / tbg_rpk$portion_bp
  tbg_rpk[,c('reads', 'bp', 'wbp', 'portion_bp')] <- NULL
  invisible(tbg_rpk)
}

calc_size_factors <- function(tbgpk = TBGpk, use_reads = "rpk", use_portion = "window") {
  # use_reads - change to "reads" to use reads per gene (as opposed to reads per gene per kb)
  # use_portion - "window" or "ORF"
  tbg_sub <- tbgpk %>% 
    filter(Portion == use_portion) %>%
    select(Locus, Sample, use_reads)
  
  # get a pseudo reference sample - the geometric mean of all rpk (read) values for the
  #  subset of genes which have reads in all samples
  pseudo_ref <- tbg_sub %>%
    group_by(Locus) %>%
    summarize(gm = geo.mean(rpk, na.rm = TRUE), n0 = sum(rpk == 0)) %>%
    filter(n0 == 0) %>%
    select(Locus, gm)
  
  # for each sample, compute the size factor as the median of: 
  #  (the sample gene rpk) / (the pseudo ref rpk) for all genes in the pseudo ref sample
  size_factors <- tbg_sub %>% group_by(Sample) %>% summarise()
  size_factors$sf <- NA
  sfno = 1
  for (s in size_factors$Sample) {
    sample_rpk <- tbg_sub %>% filter(Sample == s)
    size_factor <- pseudo_ref %>%
      left_join(sample_rpk, by = "Locus") %>%
      mutate(ratio = rpk / gm) %>%
      summarize(sf = median(ratio, na.rm = TRUE))
    size_factors$sf[sfno] <- size_factor$sf
    sfno <- sfno + 1
  }
  size_factors
}

make_rNm_column <- function(tbg = TBGpk, use_reads = "rpk", use_portion = "window") {
  size_factors <- calc_size_factors(tbg, use_reads = use_reads, use_portion = use_portion)
  tbgNm <- tbg %>%
    left_join(size_factors, by = "Sample")
  if (use_reads == "rpk") {
    tbgNm <- tbgNm %>%
      mutate(rpkNm = rpk / sf) %>%
      select(-rpk, -sf)
  } else if (use_reads == "reads") {
    tbgNm <- tbgNm %>%
      mutate(readsNm = rpk / sf) %>%
      select(-reads, -sf)
  }
  tbgNm <- tbgNm %>%
    filter(Portion == use_portion)
  invisible(tbgNm)
}

make_Npre_column <- function(tbg, pre_sample = "MRg1a",
                             use_reads = "rpkNm", use_portion = "window") {
  tbg_Pre <- tbg %>%
    filter(Sample == pre_sample) %>%
    rename(pre_vals = matches(use_reads)) %>%
    select(Locus, pre_vals)
  tbgNpre <- tbg %>%
    left_join(tbg_Pre, by = "Locus") %>%
    rename(r = eval(use_reads)) %>%
    mutate(Npre = r / pre_vals) %>%
    select(-r, -pre_vals) %>%
    rename(!!(paste(use_reads, "Npre", sep="")) := Npre)
  invisible(tbgNpre)
}

## CREATE TBGpkNm AND TBGpkNmNpre TABLES --------------------------------------

TBG_tidy <- read_tsv(files$TBG_tidy)

# TBGpkNm is table of reads normalized by sample size factors
TBGpkNm <- TBG_tidy %>%
  make_rpk_column() %>%
  make_rNm_column()

## Write intermediate file (Normalized by size factors, but not yet by pre-transformation sample)
write_tsv(TBGpkNm, files$TBGpkNm)

# TBGpkNmNp is further normalized by the pre-transformation sample (per gene)
TBGpkNmNp <- TBGpkNm %>%
  make_Npre_column()

## Write TBGpkNmNpre file
write_tsv(TBGpkNmNp, files$TBGpkNmNp)




## NORMALIZE rcmp TABLES BY SAMPLE SIZE FACTORS FOR MAKING WiG FILES ----------------------

# Calculate size factors
size_factors <- TBG_tidy %>%
  make_rpk_column() %>%
  calc_size_factors()

# Get rcmp data
rcmp_tidy <- read_tsv(files$rcmp_tidy)

# Normalize to size factors
rcmp_Nm <- rcmp_tidy %>%
  left_join(size_factors, by = "Sample") %>%
  mutate(reads_Nm = reads / sf) %>%
  select(-reads, -sf)

# Write output
write_tsv(rcmp_Nm, files$rcmp_Nm)

