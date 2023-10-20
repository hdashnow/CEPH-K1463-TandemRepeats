# Take a TRF UCSC genome annotation input file and generate a TRGT definition file
# Merge overlapping and nearby loci into compound loci
# Determine locus structure using tr-solve, falling back to TRF if tr-solve finds no motifs

import gzip
import io
import pandas as pd
import bioframe as bf
import matplotlib.pyplot as plt
import numpy as np
import math
import pysam
import sys
from trsolve import runtrsolve, rmdup

def is_gzip(file_path):
    with open(file_path, 'rb') as file:
        return file.read(2) == b'\x1f\x8b'

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

def merge_loci(clusterbed, id_suffix = ''):
    """Join simple repeats into compound loci
    Non-repetitive sequence between them will be ignored.
    
    Output format example string: 
    chr1	57367043	57367119	ID=chr1_57367043_57367119;MOTIFS=AAAAT,GAAAT;STRUC=(AAAAT)n(GAAAT)n(AAAAT)n
    """

    assert len(set(clusterbed['chrom']))==1 # check all chromosomes are the same
    # Check that none overlap?
    # Check sorted?
    
    motifs = clusterbed['sequence']
    
    uniq_motifs = ','.join(list(dict.fromkeys(motifs))) # Requires Python 3.7+ to maintain ordering
    struc = ''.join([f'({x})n'.format() for x in motifs])

    chrom = clusterbed['chrom'].values[0]
    start = clusterbed['cluster_start'].values[0]
    end = clusterbed['cluster_end'].values[0]

    loc_id = f'{chrom}_{start}_{end}'.format()
    outline = f'{chrom}\t{start}\t{end}\tID={loc_id}{id_suffix};MOTIFS={uniq_motifs};STRUC={struc}'.format()
    return outline



# Settings
in_filepath = 'data/T2T-CHM13v2.0.UCSC.SimpleRepeats.tsv.gz'
out_filepath = 'chm13v2.0_maskedY_rCRS.T2T-CHM13v2.0.UCSC.SimpleRepeats.trsolve.trgt.bed'
fastafile = '/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/ref/T2T/chm13v2.0_maskedY_rCRS.fa.gz'
trsolve = '~/tools/tr-solve-v0.3.0-linux_x86_64'

max_cluster_len = 20000 # Exclude clusters of loci > than this many bp
join_dist = 50 # If loci are within this many bp, join them into a compound locus (no interruption kept currently)

if fastafile is None:
    raise ValueError('BED input requires FASTA file')
else:
    ref = pysam.FastaFile(fastafile)

# Read input UCSC SimpleRepeats file as a bed
if is_gzip(in_filepath):
    bed_file = io.TextIOWrapper(gzip.open(in_filepath, 'rb'), encoding='utf-8')
else:
    bed_file = open(in_filepath)
        
bed = pd.read_table(bed_file)
#bed.drop(columns='#bin', inplace=True)
bed.rename(columns = {'#chrom':'chrom', 'chromStart':'start', 'chromEnd': 'end'}, inplace=True)
assert bf.core.checks.is_bedframe(bed)

print('Rows in bed: ', len(bed))

bedclusters = bf.cluster(bed, min_dist=join_dist)
print('Number of clusters:', len(bedclusters['cluster'].unique()))
bedclusters.head()

n_big = 0
n_TRF = 0
n_trsolve = 0

with open(out_filepath, 'w') as outfile:

    for cluster in bedclusters['cluster'].unique():
        clusterbed = bedclusters.loc[bedclusters['cluster'] == cluster]

        chrom = clusterbed.iloc[0]['chrom']
        start = clusterbed.iloc[0]['cluster_start']
        end = clusterbed.iloc[0]['cluster_end']
        
        if end - start > max_cluster_len:
            n_big += 1
            print('Skipping cluster of length ', end - start)
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
            #if len(clusterbed) == 1:
                # Single motif, no need to fix
            trgt_def = merge_loci(clusterbed, id_suffix = '_TRF')
            n_TRF += 1
            
        # Write in TRGT format
        outfile.write(trgt_def + '\n')

print('skipped due to size: ', n_big)
print('This many loci were solved by tr-solve: ', n_trsolve)
print('This many loci were solved by TRF: ', n_TRF)



#XXX TODO
# If any cluster overlap a known pathogenic locus, replace them with the manually curated pathogenic locus


