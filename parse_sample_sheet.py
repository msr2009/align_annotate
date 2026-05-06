"""
parse_sample_sheet.py

small script required for batch_align_annotate.sh
to figure out what the sample names are based on 
HCI genomics core's sample sheet.

(instead of hard-coding the columns, which is how I dealt with this earlier)
"""
import pandas as pd
import sys

def main(sample_sheet_file, idcolumn, namecolumn, outfile):
    dat = pd.read_csv(sample_sheet_file, header=0, sep="\t")
    out_dat_str = dat.to_csv(sep='\t', index=False, columns=[idcolumn, namecolumn], header=True)[:-1]

    if outfile != None:
        f_out = open(outfile, 'w')
        print(out_dat_str, file=f_out)
        f_out.close()
    else:
        print(out_dat_str)

    
if __name__ == '__main__':

    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('sample_sheet', type=str, help="sample sheet file")
    parser.add_argument('--id_column', type=str, dest='idcol', 
                        help="id column header (default='ID')", default='ID')
    parser.add_argument('--name_column', type=str, dest='namecol', 
                        help="name column header (default='Sample Name')", default='Sample Name')
    parser.add_argument('-', '--output_file', type=str, dest='outfilename', 
                        help="output file name (default=STDOUT)", default=None)
    args = parser.parse_args()
    main(args.sample_sheet, args.idcol, args.namecol, args.outfilename) 
