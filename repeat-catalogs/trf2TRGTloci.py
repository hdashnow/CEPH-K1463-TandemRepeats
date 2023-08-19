import gzip
import io
import pandas as pd
import bioframe as bf
import matplotlib.pyplot as plt
import numpy as np
import sys
import time
start_time = time.time()

def is_gzip(file_path):
    with open(file_path, 'rb') as file:
        return file.read(2) == b'\x1f\x8b'
    
def plot_intervals(bed):
    """Assumes sorted intervals - need to test this?
    
    bed is a pandas datafram with at least the columns chrom, start and end
    """
    assert len(set(bed['chrom']))==1 # check all chromosomes are the same
    
    intervals = list(zip(bed['start'], bed['end']))
    layers = [[intervals[0]]]
    minstart, maxend = intervals[0]
    
    for interval in intervals[1:]:
        
        if interval[0] < minstart: minstart = interval[0]
        if interval[1] > maxend: maxend = interval[1]
        for layer in layers:
            if layer[-1][-1] < interval[0]:
                layer.append(interval)
                break
        else:
            layers.append([interval])

    plt.figure()
    for i, lay in enumerate(layers):
        x1, x2 = zip(*lay)
        plt.hlines([i + 1] * len(x1), x1, x2, lw=30, color='lightblue')

    #To add text. This isn't working yet...

    # for bar in zip(x1, x2):
    #     if text is not None:
    #         plt.text(np.mean(bar), i + 1, 'test', fontsize=12,
    #                  horizontalalignment='center', verticalalignment='center', )

    plt.xlim(minstart - 2, maxend + 1)
    plt.ylim(0, len(layers) + 1)
    return plt

def fix_overlaps(clusterbed, pct_overlap = 5):
    
    if len(clusterbed) == 1:
        return(clusterbed)
    
    # If multiple, take the largest
    maxlen = (clusterbed['end'] - clusterbed['start']).max()
    best = clusterbed[clusterbed['end'] - clusterbed['start'] == maxlen]
    
    # If tied, take the highest percent match
    if len(best) > 1:
        max_match = (best['perMatch']).max()
        best = best[best['perMatch'] == max_match]

    # If still tied, take the smallest period
    if len(best) > 1:
        minperiod = (best['period']).min()
        best = best[best['period'] == minperiod]

    # If stll tied, take the highest score
    if len(best) > 1:
        maxscore = (best['score']).max()
        best = best[best['score'] == maxscore]
        
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

def merge_loci(clusterbed):
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
    outline = f'{chrom}\t{start}\t{end}\tID={loc_id};MOTIFS={uniq_motifs};STRUC={struc}'.format()
    return outline


def main(tsv: str, out: str, max_cluster_len: int = 5000, join_dist: int = 250):
    """
    :param tsv: UCSC table browser TRF SimpleRepeats in TSV format.
    :param out: Output file for TRGT locus definitions in bed-like format. 
    :param max_cluster_len: Exclude clusters of loci > than this many bp
    :join_dist: If loci are within this many bp, join them into a compound locus (no interruption kept currently)
    """

    if is_gzip(tsv):
        bed_file = io.TextIOWrapper(gzip.open(tsv, 'rb'), encoding='utf-8')
    else:
        bed_file = open(tsv)
            
    bed = pd.read_table(bed_file)
    bed.drop(columns='#bin', inplace=True)
    bed.rename(columns = {'chromStart':'start', 'chromEnd': 'end'}, inplace=True)
    assert bf.core.checks.is_bedframe(bed)

    print('Rows in bed: ', len(bed), file=sys.stderr)
    #Parse and cluster bed file
    bedclusters = bf.cluster(bed, min_dist=None)
    print('Number of clusters:', len(bedclusters['cluster'].unique()), file=sys.stderr)
    bedclusters.head()

    ## Plot some metrics about genome annotation

    ### Gaps between loci
    bed['direction'] = '.'
    # # Calculate distance to closest repeat
    # closest_intervals = bf.closest(bed, df2=None, return_distance=True, direction_col = 'direction')
    # display(closest_intervals)

    n_big = 0

    # Fix overlapping intervals by choosing the largest then trimming the others
    all_list = []
    for cluster in bedclusters['cluster'].unique():
        clusterbed = bedclusters.loc[bedclusters['cluster'] == cluster]
        
        if clusterbed.iloc[0]['cluster_end'] - clusterbed.iloc[0]['cluster_start'] > max_cluster_len:
            n_big += 1
            continue

        if len(clusterbed) == 1:
            all_list.append(clusterbed)
            
        if len(clusterbed) > 1:
            non_overlapping = fix_overlaps(clusterbed)
            all_list.append(non_overlapping)
        
    all_non_overlappping = pd.concat(all_list, ignore_index=True, copy=True)
        
    # Cluster loci within x bp of each other

    # Drop old cluster columns
    all_non_overlappping.drop(columns=['cluster', 'cluster_start', 'cluster_end'], inplace = True)
    bedclusters2 = bf.cluster(bed, min_dist=join_dist)
    
    # Join nearby loci into compound ones
    n_big = 0
    merged = 0
    max_merge = 0
    biggest_cluster = None

    metrics = open('metrics.txt', 'w')
    metrics.write('chrom\tstart\tlength\tmotifs\n')
    with open(out, 'w') as outfile:
        for cluster in bedclusters2['cluster'].unique():
            clusterbed = bedclusters2.loc[bedclusters2['cluster'] == cluster]

            chrom = clusterbed['chrom'].values[0]
            start = clusterbed['cluster_start'].values[0]
            length = clusterbed['cluster_end'].values[0] - clusterbed['cluster_start'].values[0]
            motifs = len(clusterbed)
            metrics.write(f'{chrom}\t{start}\t{length}\t{motifs}\n')
            
            if clusterbed['cluster_end'].values[0] - clusterbed['cluster_start'].values[0] > max_cluster_len:
                n_big += 1
                continue
            
            if len(clusterbed) > 1:
                merged += 1
                if len(clusterbed) > max_merge:
                    max_merge = len(clusterbed)
                    biggest_cluster = clusterbed

            outfile.write(merge_loci(clusterbed) + '\n')
            
            
                #all_merged = pd.concat([all_merged, merged_loci], ignore_index=True, copy=True)
        
    print(f'skipped due to size > {max_cluster_len}: {n_big}', file=sys.stderr)
    print('This many loci are compound: ', merged, file=sys.stderr)
    #print('Biggest cluster:\n', biggest_cluster, file=sys.stderr)

    # If any cluster overlap a known pathogenic locus, replace them with the manually curated pathogenic locus


if __name__ == "__main__":
    import doctest
    doctest.testmod()
    import defopt
    defopt.run(main)

    elapsed = time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))
    print(f'Elapsed time: {elapsed}', file=sys.stderr)

