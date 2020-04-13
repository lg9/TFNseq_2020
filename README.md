# TFNseq_2020
Bioinformatic scripts used for TFNseq analysis of essential gene depletion sensitivities in A. baylyi ADP1.  See Gallagher-LA, Bailey-J and Manoil-C. 2020. Ranking essential bacterial processes by speed of mutant death.

TFNseq data processing and analysis workflow

	1. Preliminary processing of sample fastq files.
	(Filter reads for those with transposon end sequence, trim and map remainder of each read against genome. Count reads per unique insertion location. Tabulate read counts for each sample at each insertion location. Tabulate total hits and reads per gene, both within the entire ORF and within its functional window.)
		
		Use python analysis scripts available at https://github.com/lg9/Tn-seq
		Annotation files (available in 005_annot_files_for_prelim_processing/):
			CR543861.fna (genome nucleotide sequence)
			CR543861.ptt (ORF annotation file)
			CR543861w_180809.ptt (functional windows annotation file)
		
		a. For each sample fastq file, run process_map.py
		
			E.g.:  $ python (path-to)/Tn-seq/python/process_map.py -r CR543861.fna -j -t AGACAG -n 0 -s -w work (path-to)/MRgDNA_sra.fastq
			
			Output metrics should match those reported in Table S1
			
		b.	Run process_tabulate_annotate.py twice, as follows:
		
			i. First, using annotation of the 'functional windows'
				
			$ python (path-to)/Tn-seq/python/process_annotate_tabulate.py -r CR543861.fna -a CR543861w_180809.ptt -o rcmp.xls -p TbGW.xls -w work work/MRgDNA_sra_trim_sum_mg.txt work/primary_2h_sra_trim_sum_mg.txt work/primary_4h_sra_trim_sum_mg.txt work/primary_6h_sra_trim_sum_mg.txt work/primary_8h_sra_trim_sum_mg.txt work/primary_18h_AyL4_sra_trim_sum_mg.txt work/secondary_2h_sra_trim_sum_mg.txt work/secondary_4h_sra_trim_sum_mg.txt work/secondary_6h_sra_trim_sum_mg.txt work/secondary_8h_sra_trim_sum_mg.txt work/secondary_18h_AyL3_sra_trim_sum_mg.txt
				
			ii. Second, using the annotation of the ORFs
			(This overwrites the rcmp.xls file created previously, but the TbGn.xls file will not overwrite the previous TbGW.xls file.)
				
			$ python (path-to)/Tn-seq/python/process_annotate_tabulate.py -r CR543861.fna -a CR543861.ptt -o rcmp.xls -p TbGn.xls -w work work/MRgDNA_sra_trim_sum_mg.txt work/primary_2h_sra_trim_sum_mg.txt work/primary_4h_sra_trim_sum_mg.txt work/primary_6h_sra_trim_sum_mg.txt work/primary_8h_sra_trim_sum_mg.txt work/primary_18h_AyL4_sra_trim_sum_mg.txt work/secondary_2h_sra_trim_sum_mg.txt work/secondary_4h_sra_trim_sum_mg.txt work/secondary_6h_sra_trim_sum_mg.txt work/secondary_8h_sra_trim_sum_mg.txt work/secondary_18h_AyL3_sra_trim_sum_mg.txt
				
	2. Make tidy data tables of the output files from step 1.
	
		a. Make simple versions of the output files
		
			i. run simplify_rcmp.py
				Input rcmp file name: 010_raw_rcmp_tbgn_files/rcmp.xls
				Output q0 columns [y/n]? n
				Output qn0 columns [y/n]? n
				Output file name: 020_simple_rcmp_tbgn_files/rcmp_simple_aQonly.xls 
				(Short names:  MRg1a, __ post_tfm __, T1_2h, T1_4h, T1_6h, T1_8h, AyL4, T2_2h, T2_4h, T2_6h, T2_8h, AyL3)
				
			ii. run simplify_tbgn.py twice
				a.	Input TbGn file name: 010_raw_rcmp_tbgn_files/TbGn.xls
					Output qn0 columns [y/n]? n
					Output file name: 020_simple_rcmp_tbgn_files/TbGn_simple_aQonly.xls
				b.	Input TbGn file name: 010_raw_rcmp_tbgn_files/TbGW.xls
					Output qn0 columns [y/n]? n
					Output file name: 020_simple_rcmp_tbgn_files/TbGW_simple_aQonly.xls
					
				(Short names as above)
				
		b. Make tidy tables from the simple versions
		
			i. run R005_make_tidy_rcmp_table.R
			
			ii. run make_tidy_TBG_table.py
				TbGn_simple_aQonly file name: 020_simple_rcmp_tbgn_files/TbGn_simple_aQonly.xls
				TbGW_simple_aQonly file name: 020_simple_rcmp_tbgn_files/TbGW_simple_aQonly.xls
				Output file name: 030_tidy_rcmp_TBG_files/TBG_tidy.xls
				
				(The output TBG file combines the TbGn and the TbGW data to produce a single table of hits and reads for both ORFs and functional windows.  The overall workflow followed to obtain the TBG file was used for legacy reasons; there are certainly more efficient ways to achieve the same result.)
				
	3. Normalize reads per sample by size factors as described in Methods.
	
		a. run R010_normalization.R
		
			Output tables:
				TBGpkNm - read counts per gene (per kb) normalized by sample size factors
				TBGpkNmNp - read counts per gene further normalized by the pre-transformation sample counts for the same gene
				rcmp_Nm - read counts per site normalized by the same sample size factors; useful for making .wig files and plotting reads per site (e.g. Fig 2b).
				
	4. Define essential genes.
	
		a. run R020_define_essentials.R
		
			Output tables:
				essentials_final - list of genes defined as essential (using criteria described in Methods)
				essentiality_by_sample - tabulation of 'essentiality' (significant depletion) in each sample (time point)
	
	5. Calculate depletion times and ranks.
	
		a. run R030_calculate_depl_times_ranks.R
		
			Output tables:
				d_t_r_all - depletion times and ranks (Trial 1 and Trial 2) for all genes
				d_t_r_essentials - depletion times and ranks (Trial 1 and Trial 2) for essential genes
	
	6. Make spread data tables for reporting.
	
		a. run R040_make_spread_tables.R
		
			Output tables:
				spread_gene_data_all - spread tabulation of hits, normalized rpk, delta log reads values, 'essentiality' (significant depletion) and depletion times and ranks for all samples and all genes.
				spread_gene_data_ess - as above but of essential genes only and lacking 'essentiality' calls.  This is the data for Dataset S1.
				
	7. Make boxplot and associated supplementary table (for Fig 4 and Table S2).

		a. run R050_make_processes_boxplot_and_table.R
		
			Output:
				processes_summary_table.xls - containing data for Table S2
				boxplot.pdf - core plot for Fig 4
				
	8. Make kinetic depletion plots for specific processes (for Fig 3).
	
		a. run R060_make_kinetic_plots.R
		
			Output:  individual panels used to make Fig 3.  
			Kinetic plots of other processes can be created using alternative plot_def_files.

	9. Make plots of progressive loss of essentials and hits-per-loci (for Fig 2).
	
		a. run R070_plot_essentials_progressive_depletion.R
		
			Output:  progressive depletion plot used to make Fig 2a
			
		b. run R075_plot_hits_at_locus.R
		
			Output:  plots of hits at specific loci - used to make Fig 2b
		
	10. Make validation plots
		
		a. run R080_make_validation_plots.R
		
			Output:
				plot of Trial2 vs Trial1 depletion times - for Fig S1a
				plot of depletion time vs qPCR depletion rate for selected deletions - for Fig S1b
				
	11. Make plots of microcolony area and cell count measurements
	
		a. run R090_plot_microcolony_areas.R
		
			Output:
				plot of microcolony median areas vs depletion time - for Fig S2a
				plot of microcolony cell count vs median area - for Fig S2b
