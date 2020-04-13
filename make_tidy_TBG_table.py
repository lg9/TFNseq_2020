'''
Make tidy TBG table

Takes as input:
- a TbGn_simple_aQonly.xls file
- a TbGW_simple_aQonly.xls file

Produces an output table to be read as a tidy tibble in R.  The output table
has the following variables:
- locus (e.g., ACIAD0001, aciad2345)
- sample (e.g., MRg1a, T3_4h)
- ORF portion ('ORF', 'window')
- hits (integer)
- reads (double)

'''

## PARAMETERS

## PROGRAM

tbgn_file = input('TbGn_simple_aQonly file name: ')
tbgW_file = input('TbGW_simple_aQonly file name: ')
# choose output file name
out_filename = input('Output file name: ')

# Process input files
with open(tbgn_file, 'r') as fin_n, open(tbgW_file, 'r') as fin_W:
    h = fin_n.readline()
    fin_W.readline()
    h_fields = h.rstrip().split('\t')
    # get the field numbers for locus and first sample
    locus_hf_no = h_fields.index('Locus')
    sample1_hf_no = [ fno for fno in range(len(h_fields)) if h_fields[fno].startswith('hTot') ][0]
    # get the sample names from h
    hTot_aQ_col_names = [ hf for hf in h_fields if hf.startswith('hTot') and hf.endswith('_aQ') ]
    sample_names = [ cn[5:-3] for cn in hTot_aQ_col_names ]
    sample_count = len(sample_names)

    # Now read the TbGn and TbGW files to get the data and write to output file
    with open(out_filename, 'w') as fout:
        lout = '\t'.join(['Locus', 'Sample', 'Portion', 'hits', 'reads']) + '\n'
        fout.write(lout)
        while True:
            tbgn_line = fin_n.readline()
            tbgW_line = fin_W.readline()
            if not tbgW_line:
                break
            tbgn_fields = tbgn_line.rstrip().split('\t')
            tbgW_fields = tbgW_line.rstrip().split('\t')

            loc = tbgW_fields[locus_hf_no]
            if loc.endswith('w'):
                loc = loc[:-1]

            for sn in range(sample_count):
                samp_hT_fno = sample1_hf_no + 4 * sn
                out_fields1 = [loc, sample_names[sn], 'ORF'] + tbgn_fields[samp_hT_fno : samp_hT_fno + 2]
                out_fields2 = [loc, sample_names[sn], 'window'] + tbgW_fields[samp_hT_fno : samp_hT_fno + 2]
                lout1 = '\t'.join(out_fields1) + '\n'
                lout2 = '\t'.join(out_fields2) + '\n'
                fout.write(lout1)
                fout.write(lout2)

print('Done')


