## Make tidy rcmp file - from rcmp_simple_aQonly.xls file

library(tidyverse)

source('R001_common.R')

## ----------------------------------------------------------------------------

# get rcmp data
rcmp_simple <- read_tsv(files$rcmp_simple)

# make stripped-down tidy table
stripped <- rcmp_simple %>%
  select(
      EffPos, Dir, contains("_aQ")
    )
colnames(stripped) <- str_replace(colnames(stripped), "_aQ", "")

rcmp_tidy <- stripped %>%
  gather(-EffPos, -Dir, key = "Sample", value = "reads") %>%
  filter(reads > 0)

## Save the tidy file
write_tsv(rcmp_tidy, files$rcmp_tidy)
