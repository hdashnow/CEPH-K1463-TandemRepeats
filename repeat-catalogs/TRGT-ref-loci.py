# Take a TRF UCSC genome annotation input file and generate a TRGT definition file
# Merge overlapping and nearby loci into compound loci
# Determine locus structure using tr-solve, falling back to TRF if tr-solve finds no motifs

import gzip
import io
import pandas as pd
import bioframe as bf
import matplotlib.pyplot as plt
#import numpy as np
import pysam
import sys
from trsolve import runtrsolve, rmdup

def is_gzip(file_path):
    with open(file_path, 'rb') as file:
        return file.read(2) == b'\x1f\x8b'

def overlap(a, b):
    if min(a[1], b[1]) - max(a[0], b[0]) > 0:
        return True

def fix_overlaps(clusterbed, pct_overlap = 5):
    
    if len(clusterbed) == 1:
        return(clusterbed)
    
    # If multiple, take the largest
    maxlen = (clusterbed['end'] - clusterbed['start']).max()
    best = clusterbed[clusterbed['end'] - clusterbed['start'] == maxlen]
    
    # If tied, take the highest score
    if len(best) > 1:
        maxscore = (best['score']).max()
        best = best[best['score'] == maxscore]
        
    # If still tied, take the smallest period
    if len(best) > 1:
        minperiod = (best['period']).min()
        best = best[best['period'] == minperiod]
        
    # If still tied, take the first
    if len(best) > 1:
        best = best.head(1)
    
    assert len(best) == 1
    
    # Check if any intervals remain that don't overlap the largest, or only overlap by less than x%
    # Not sure if this is needed at all?
#     overlaps = bf.coverage(clusterbed, best)
#     overlaps['pct_overlap'] = overlaps['coverage']/(clusterbed['end'] - clusterbed['start'])*100
#     alsokeep = overlaps[overlaps['pct_overlap'] < pct_overlap]
    
    # Trim overlapping sequence
    trimmed = bf.subtract(clusterbed, best)

    if len(trimmed) > 0:
        return bf.sort_bedframe(pd.concat([best, trimmed], join='inner', ignore_index=True, copy=True))
    else:
        return best.copy()

def merge_loci(clusterbed, id_suffix = '', maxmotifs = 5):
    """Join simple repeats into compound loci
    Non-repetitive sequence between them will be ignored.
    
    Output format example string: 
    chr1	57367043	57367119	ID=chr1_57367043_57367119;MOTIFS=AAAAT,GAAAT;STRUC=(AAAAT)n(GAAAT)n(AAAAT)n
    """

    assert len(set(clusterbed['chrom']))==1 # check all chromosomes are the same
    # Check that none overlap?
    # Check sorted?
    
    if clusterbed.shape[0] > maxmotifs:
        widths = clusterbed['end'] - clusterbed['start']
        clusterbed = clusterbed.loc[widths.nlargest(maxmotifs).index.sort_values()]

    motifs = clusterbed['sequence']
    uniq_motifs = ','.join(list(dict.fromkeys(motifs))) # Requires Python 3.7+ to maintain ordering
    struc = ''.join([f'({x})n'.format() for x in motifs])

    chrom = clusterbed['chrom'].values[0]
    start = clusterbed['cluster_start'].values[0]
    end = clusterbed['cluster_end'].values[0]

    loc_id = f'{chrom}_{start}_{end}'.format()
    outline = f'{chrom}\t{start}\t{end}\tID={loc_id}{id_suffix};MOTIFS={uniq_motifs};STRUC={struc}'.format()
    return outline

def main(repeats: str, fasta: str, output: str, *, replace: str = None, replace_suffix: str = '_pathogenic', trsolve: str = 'tr-solve',
         max_cluster_len: int = 10000, join_dist: int = 50, exclude_ends: int = 250):
    """
    :param repeats: UCSC TRF SimpleRepeats file (tsv)
    :param output: Output file name
    :param fasta: Reference genome fasta file
    :param replace: TRGT file containing manually curated loci to replace overlapping loci with
    :param replace_suffix: Annotate replacement loci IDs with this suffix
    :param trsolve: Path to tr-solve executable
    :param max_cluster_len: Exclude clusters of loci > than this many bp
    :param join_dist: If loci are within this many bp, join them into a compound locus (no interruption kept currently)
    :param exclude_ends: Exclude loci within this many bp of the chromosome ends
    """
    # Clarify variable names
    in_filepath = repeats
    out_filepath = output


    if fasta is None:
        raise ValueError('BED input requires FASTA file')
    else:
        ref = pysam.FastaFile(fasta)

    chrom_lengths = {k: v for k, v in zip(ref.references, ref.lengths)}

    # Read replacement loci
    if replace:
        replacebed = pd.read_table(replace, header=None, 
                                   names=['chrom', 'start', 'end', 'definition'])
        sys.stderr.write(f'Read {len(replacebed)} replacement loci\n')

    # Read input UCSC SimpleRepeats file as a bed
    if is_gzip(in_filepath):
        bed_file = io.TextIOWrapper(gzip.open(in_filepath, 'rb'), encoding='utf-8')
    else:
        bed_file = open(in_filepath)
            
    bed = pd.read_table(bed_file)
    if '#bin' in bed.columns:
        bed.drop(columns='#bin', inplace=True)
    if '#chrom' in bed.columns:
        bed.rename(columns={'#chrom':'chrom'}, inplace=True)
    bed.rename(columns={'chromStart':'start', 'chromEnd': 'end'}, inplace=True)
    assert bf.core.checks.is_bedframe(bed)

    sys.stderr.write(f'Rows in bed: {len(bed)}\n')

    bedclusters = bf.cluster(bed, min_dist=join_dist)
    sys.stderr.write(f'Rows in bedclusters: {len(bedclusters)}\n')
    bedclusters.head()

    n_big = 0
    n_TRF = 0
    n_trsolve = 0

    with open(out_filepath, 'w') as outfile:

        position_last = 0
        for cluster in bedclusters['cluster'].unique():
            keep_cluster = True
            clusterbed = bedclusters.loc[bedclusters['cluster'] == cluster]

            chrom = clusterbed.iloc[0]['chrom']
            start = clusterbed.iloc[0]['cluster_start']
            end = clusterbed.iloc[0]['cluster_end']

            if replace:
                if chrom in replacebed['chrom'].values:
                   
                    for _, row in replacebed[replacebed['chrom'] == chrom].iterrows():
                        # Check for overlap
                        if overlap((row['start'], row['end']), (start, end)):
                            sys.stderr.write(f'Overlapping replacement loci found for {chrom}:{start}-{end}\n')
                            sys.stderr.write(f'Replacing it with {row["definition"]}\n')
                            definition = row['definition'].replace(';MOTIFS=', replace_suffix + ';MOTIFS=')
                            outfile.write(f'{row["chrom"]}\t{row["start"]}\t{row["end"]}\t{definition}\n')
                            # Remove the locus from replacebed so it doesn't get written again
                            replacebed.drop(row.name, inplace=True)
                            keep_cluster = False
                    
            if keep_cluster and chrom in chrom_lengths: # Check if chrom is in reference genome
                if start < exclude_ends or end > chrom_lengths[chrom] - exclude_ends:
                    sys.stderr.write(f'Skipping cluster at {chrom}:{start}-{end} due to proximity to chromosome end\n')
                    continue
                
                if end - start > max_cluster_len:
                    n_big += 1
                    sys.stderr.write(f'Skipping cluster of length {end - start}\n')
                    continue

                try:
                    sequence = ref.fetch(chrom, start, end)
                except KeyError as e:
                    sys.stderr.write(f'Error: {e}\n')
                    continue
        
                # Run tr-solve
                try:
                    motifs, bounds = runtrsolve(sequence, trsolve=trsolve)
                except OSError as e:
                    sys.stderr.write(f'Skipping region {chrom}:{start}-{end} of length {end - start} due to error: {e}\n')
                    n_big += 1
                    continue

                if len(motifs) > 0:
                    unique_motifs = dict.fromkeys(motifs) # can use set(motifs) if order doesn't matter. Assume it does for now
                    motifs_str = ','.join(unique_motifs)
                    struc = ''.join([f'({motif})n' for motif in rmdup(motifs)])
                    trgt_def = f'{chrom}\t{start}\t{end}\tID={chrom}_{start}_{end}_trsolve;MOTIFS={motifs_str};STRUC={struc}'
                    n_trsolve += 1

                else:
                    # No motifs found, use TRF
                    trgt_def = merge_loci(clusterbed, id_suffix = '_TRF')
                    n_TRF += 1
                    
                # Write in TRGT format
                outfile.write(trgt_def + '\n')

            # Update position_last
            position_last = end

        # Write remaining replacement loci, if any
        if replace:
            if len(replacebed) > 0:
                sys.stderr.write(f'{len(replacebed)} replacement loci being added to end of file\n')
                for _, row in replacebed.iterrows():
                    definition = row['definition'].replace(';MOTIFS=', replace_suffix + ';MOTIFS=')
                    outfile.write(f'{row["chrom"]}\t{row["start"]}\t{row["end"]}\t{definition}\n')

    sys.stderr.write(f'Skipped due to size: {n_big}\n')
    sys.stderr.write(f'Solved by tr-solve: {n_trsolve}\n')
    sys.stderr.write(f'Solved by TRF: {n_TRF}\n')

if __name__ == "__main__":
    import defopt
    defopt.run(main)