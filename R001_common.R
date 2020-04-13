## common parameters and tables


## PATHS AND FILES ------------------------------------------------------------

files = list(
  
  annot = '040_other_source_tables/ADP1_annotation.xls',
  
  rcmp_simple = '020_simple_rcmp_tbgn_files/rcmp_simple_aQonly.xls',
  
  rcmp_tidy = '030_tidy_rcmp_TBG_files/rcmp_tidy.xls',
  
  TBG_tidy = '030_tidy_rcmp_TBG_files/TBG_tidy.xls',
  
  TBGpkNm = '050_transformed_data_tables/TBGpkNm.xls',
  TBGpkNmNp = '050_transformed_data_tables/TBGpkNmNp.xls',
  
  rcmp_Nm = '050_transformed_data_tables/rcmp_Nm.xls',
  
  essentials_final = '050_transformed_data_tables/essentials_final.xls',
  essentiality_by_sample = '050_transformed_data_tables/essentiality_by_sample.xls',
  
  depl_times_ranks_all = '050_transformed_data_tables/depl_times_ranks_all_genes.xls',
  depl_times_ranks_essls = '050_transformed_data_tables/depl_times_ranks_essentials.xls',
  
  subproc_loci_list = '040_other_source_tables/subprocesses_loci_list.xls',
  subproc_order_incl = '040_other_source_tables/processes_codes_order.xls',
  
  spread_gene_data_all = '060_output_tables/spread_gene_data_rpkNm_all.xls',
  spread_gene_data_ess = '060_output_tables/spread_gene_data_rpkNm_essentials.xls',
  
  processes_summary_table = '060_output_tables/processes_summary_table.xls',
  
  qPCR_SQ_table = '040_other_source_tables/qPCR_relSQ_data.csv',
  
  microcol_areas_counts = '040_other_source_tables/microcolony_areas_counts.xls'
)

## PARAMETERS -----------------------------------------------------------------

outgr_times <- tibble(
  Sample = c('MRg1a', 'T1_2h', 'T1_4h', 'T1_6h', 'T1_8h',
             'T2_2h', 'T2_4h', 'T2_6h', 'T2_8h', 'AyL3', 'AyL4'),
  time = c(0, 2.27, 4.17, 6.10, 8.35,
           2.37, 4.23, 6.13, 8.33, 12.21, 15.40)
  # For AyL4 and AyL3, used extrapolated outgrowth times based on total doublings, rather than 18 h
)

## LOAD COMMON TABLES ---------------------------------------------------------

## Annotation Table
annot <- read_tsv(files$annot)
annot <- annot %>%
  mutate_at(
    3:8, as.integer
  ) %>%
  mutate(
    Label = ifelse(nchar(annot$Gene) > 1, annot$Gene, as.character(annot$Locus))
  )
