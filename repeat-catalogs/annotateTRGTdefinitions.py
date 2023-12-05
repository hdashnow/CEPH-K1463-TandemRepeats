
import pandas as pd
import os
import sys
import TRlib

def main(trgt_definitions: str, fasta: str, *, output: str = 'stdout', bed_files: list[str] = None, bed_names: list[str] = None):
        """
        :param trgt_definitions: Path to the TRGT definitions bed file
        :param fasta: Path to the reference FASTA file
        :param bed_files: Paths to one or more bed files
        :param bed_names: Column names for the bed files (must be the same length as bed_files)
        :output: bed files with annotated TRGT definitions
        """
        if output != 'stdout':
            outstream = open(output, 'w')
            sys.stdout = outstream
        else:
            outstream = sys.stdout

        header = ['chrom', 'start', 'end', 'TRid',
                    'longest_homopolymer', 'gc_content', 
                    'n_motifs', 'min_motiflen', 'max_motiflen',
                    'motifs', 'struc',
                ]

        if bed_files:
            if bed_names:
                if len(bed_names) != len(bed_files):
                    raise ValueError('bed_names must be the same length as bed_files')
                header += bed_names
            else:
                header += [os.path.splitext(os.path.basename(bed_file))[0] for bed_file in bed_files]
        
        outstream.write('#' + '\t'.join(header) + '\n')

        # Read the TRGT bed file into a pandas DataFrame
        trgt_df = pd.read_csv(trgt_definitions, sep='\t', header=None)

        # Extract the sequence from the FASTA file using the locus coordinates
        for (chrom, start, end, info), sequence in TRlib.extractfasta(trgt_definitions, fasta, maxout=None, returninfo=True):

            TRid = info.split(';')[0].split('=')[1]
            motifs = info.split(';')[1].split('=')[1]
            motifslist = motifs.split(',')
            struc = info.split(';')[2].split('=')[1]
            longest_homopolymer = TRlib.longesthomopolymer(sequence)
            gc_content = TRlib.GCcontent(sequence)
            n_motifs = len(motifslist)
            min_motiflen = len(min(motifslist, key=len))
            max_motiflen = len(max(motifslist, key=len))


            outstream.write('\t'.join([str(x) for x in [chrom, start, end, TRid,
                                       longest_homopolymer, gc_content, 
                                       n_motifs, min_motiflen, max_motiflen,
                                       motifs, struc,
                                    ]]) + '\n')

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    import defopt
    defopt.run(main)