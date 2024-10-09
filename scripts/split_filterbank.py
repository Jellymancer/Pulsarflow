import filtools
import argparse 
import os


def main(args):

    input = filtools.FilterbankIO()
    input.read_header(args.fil)
    
    nsamples = input.total_samples()
    chunksize = nsamples//args.nblocks

    # Split the filterbank file into subbands
    for i in range(args.nblocks):
        start_sample = i*chunksize
        end_sample = (i+1)*chunksize

        input.seek_to_sample(i*chunksize)
        data = input.read_block(chunksize)

        # Determine the output filename
        if args.outdir is None:
            output_filename = f"{i}_" + os.path.splitext(args.fil)[0] + f"_{start_sample}_{end_sample}.fil"
        else:
            output_filename = f"{i}_" + os.path.join(args.outdir, os.path.splitext(args.fil)[0] + f"_{start_sample}_{end_sample}.fil")

        output = filtools.FilterbankIO()
        output.header = input.header
        output.write_header(output_filename)
        output.write_block(data)

parser = argparse.ArgumentParser(description='Split a filterbank file into subbands.')
parser.add_argument('fil', metavar='fil', type=str,
                    help='The filterbank file to split.')
parser.add_argument('--nblocks', type=int,
                    help='The number of subbands to split the filterbank file into.')
parser.add_argument('--outdir', type=str, default=None,
                    help='The directory to save the split filterbank files. If' 
                    'not provided, the files will be saved in the same directory as the input file.')

args = parser.parse_args()

if __name__ == "__main__":
    main(args)








