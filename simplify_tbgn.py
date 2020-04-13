'''
Simplify TbGn file
'''

## PARAMETERS

# initial header fields to keep
field_nos_hit_info = [0, 1, 3, 4, 5, 7]   # Replicon, EffPos, Dir, Locus, Strand, Gene, RelPos
new_h_fields = ['Replicon', 'Locus', 'From', 'To', 'Strand', 'Gene']
# number of first sample column
samp1_field_no = 13     # after hlocs_allrunds fields
# sample set column prefixes
sample_col_prefixes = ['hTot_', 'rTot_', 'h5-90_', 'r5-90_']

## PROGRAM

tbgn_file = input('Input TbGn file name: ')
response = input('Output qn0 columns [y/n]? ')
# make qn0 columns?
make_qn0 = False
if (response[0] == 'y' or response[0] == 'Y'):
    make_qn0 = True
# choose output file name
out_filename = input('Output file name: ')

# Process input file
with open(tbgn_file, 'r') as fin:
    
    # Make new header line
    h = fin.readline()
    hfields = h.strip().split('\t')
    # calcuate field nos for aQ and q0 columns
    no_of_samples = int( ( len(hfields) - samp1_field_no ) / 8 )    # 4 columns per sample per aQ or Q0
    field_nos_aQ = [ i for i in range( samp1_field_no, samp1_field_no + 4 * no_of_samples ) ]
    # choose short sample names
    hTot_aQ_fields = [ hfields[i] for i in range( samp1_field_no, samp1_field_no + 4 * no_of_samples, 4 ) ]
    short_sample_names = [input("Short name for " + sn[ 5 : sn.find('_all')] + ": ") for sn in hTot_aQ_fields ]
    # make new header
    new_h_fields += [ pf + s + '_aQ' for s in short_sample_names for pf in sample_col_prefixes ]
    if (make_qn0):
        new_h_fields += [ pf + s + '_Qn0' for s in short_sample_names for pf in sample_col_prefixes ]
    newh = '\t'.join(new_h_fields) + '\n'

    # Process old file, Make new file
    with open(out_filename, 'w') as fout:
        fout.write(newh)
        for line in fin.readlines():
            infields = line.strip().split('\t')
            outfields = [infields[i] for i in field_nos_hit_info + field_nos_aQ]
            if (make_qn0):
                qn0_vals = [ float(infields[i]) - float(infields[i + 4 * no_of_samples]) for i in field_nos_aQ ]
                qn0_vals = [ int(v) if (int(v) == float(v)) else v for v in qn0_vals ]
                outfields += [ str(v) for v in qn0_vals ]
            outline = '\t'.join(outfields) + '\n'
            fout.write(outline)

print('Done')


