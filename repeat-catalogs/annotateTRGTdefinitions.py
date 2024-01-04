import pandas as pd
import bioframe as bf
import os
import sys
import TRlib

def annotateFasta(trgt_definitions: str, fasta: str):
        """Extract the sequence from the FASTA file using the locus coordinates
        and annotate the TRGT definitions with the sequence and other information."""
        for (chrom, start, end, info), sequence in TRlib.extractfasta(trgt_definitions, fasta, maxout=None, returninfo=True):

            TRid = info.split(';')[0].split('=')[1]
            motifs = info.split(';')[1].split('=')[1]
            motifslist = motifs.split(',')
            struc = info.split(';')[2].split('=')[1]
            longest_homopolymer = TRlib.longesthomopolymer(sequence)
            if longest_homopolymer is None:
                longest_homopolymer = 'NA'
            gc_content = TRlib.GCcontent(sequence)
            if gc_content is None:
                gc_content = 'NA'
            n_motifs = len(motifslist)
            min_motiflen = len(min(motifslist, key=len))
            max_motiflen = len(max(motifslist, key=len))

            result = [chrom, int(start), int(end), TRid,
                        longest_homopolymer, gc_content, 
                        n_motifs, min_motiflen, max_motiflen,
                        motifs, struc,
                    ]
            yield result

def main(trgt_definitions: str, fasta: str, *, output: str = 'stdout', bed_files: list[str] = None, bed_names: list[str] = None):
        """
        Outputs a bed file with annotated TRGT definitions including length of longest homopolymer, 
        GC content, number of motifs, min and max motif length, and motifs.
        If annotation bed files are provided, number of bp overlap is reported.

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
            else:
                bed_names = [os.path.splitext(os.path.basename(bed_file))[0] for bed_file in bed_files]

        trgt_df = pd.DataFrame(columns=header, data=annotateFasta(trgt_definitions, fasta))

        # Read the bed files into a dictionary of pandas DataFrames
        bed_dfs = {}
        if bed_files:
            for bed_name, bed_file in zip(bed_names,bed_files):
                bed_dfs[bed_name] = pd.read_csv(bed_file, sep='\t', header=None)
                colnames = [str(x) for x in bed_dfs[bed_name].columns]
                colnames[:3] = ['chrom', 'start', 'end']
                bed_dfs[bed_name].columns = colnames

        # Overlap the bed files with the TRGT DataFrame
        for bed_name, bed_df in bed_dfs.items():
            trgt_df = bf.coverage(trgt_df, bed_df)
            trgt_df = trgt_df.rename(columns={'coverage':bed_name})

        # Write the TRGT definitions with annotations to the output file
        trgt_df = trgt_df.rename(columns={'chrom':'#chrom'})
        trgt_df.to_csv(outstream, sep='\t', header=True, index=False)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
    import defopt
    defopt.run(main)