'''
Simplify rcmp file

Provides options to include either q0 or qn0 columns (or anyQ only)
'''

## PARAMETERS

# initial header fields to keep
field_nos_hit_info = [0, 2, 3, 5, 6, 7, 9]   # Replicon, EffPos, Dir, Locus, Strand, Gene, RelPos
new_h_fields = ['Replicon', 'EffPos', 'Dir', 'Locus', 'Strand', 'Gene', 'RelPos']
# number of first sample column
samp1_field_no = 11

## PROGRAM

rcmp_file = input('Input rcmp file name: ')
make_q0 = False
make_qn0 = False
response1 = input('Output q0 columns [y/n]? ')
# make q0 columns?
if (response1[0] == 'y' or response1[0] == 'Y'):
    make_q0 = True
else:
    response2 = input('Output qn0 columns [y/n]? ')
    # make qn0 columns?
    if (response2[0] == 'y' or response2[0] == 'Y'):
        make_qn0 = True
# choose output file name
out_filename = input('Output file name: ')

# Process input file
with open(rcmp_file, 'r') as fin:
    
    # Make new header line
    h = fin.readline()
    hfields = h.strip().split('\t')
    # calcuate field nos for aQ and q0 columns
    no_of_samples = int( ( len(hfields) - samp1_field_no ) / 2 )
    field_nos_aQ = [ i for i in range( samp1_field_no, samp1_field_no + no_of_samples ) ]
    # choose short sample names
    short_sample_names = [input("Short name for " + hfields[sn][:hfields[sn].find('_all')] + ": ") for sn in field_nos_aQ ]
    # make new header
    new_h_fields += [ s + '_aQ' for s in short_sample_names ]
    if (make_q0):
        new_h_fields += [ s + '_q0' for s in short_sample_names ]
    elif (make_qn0):
        new_h_fields += [ s + '_Qn0' for s in short_sample_names ]
    newh = '\t'.join(new_h_fields) + '\n'
    
    # Process old file, Make new file
    with open(out_filename, 'w') as fout:
        fout.write(newh)
        for line in fin.readlines():
            infields = line.strip().split('\t')
            outfields = [infields[i] for i in field_nos_hit_info + field_nos_aQ]
            if (make_q0):
                outfields += [ infields[i + no_of_samples] for i in field_nos_aQ ]
            elif (make_qn0):
                qn0_vals = [ float(infields[i]) - float(infields[i + no_of_samples]) for i in field_nos_aQ ]
                qn0_vals = [ int(v) if (int(v) == float(v)) else v for v in qn0_vals ]
                outfields += [ str(v) for v in qn0_vals ]
            outline = '\t'.join(outfields) + '\n'
            fout.write(outline)

print('Done')


