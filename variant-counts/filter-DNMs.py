import pandas as pd
import glob
import sys

# # Find files ending with .gz using glob
# sample_files = glob.glob('data/TRGTdn/*.gz')#[0:1]
# print(sample_files)

usecols = ['trid', 'denovo_coverage','per_allele_reads_father', 'per_allele_reads_mother',
           'per_allele_reads_child', 'index', 'father_MC', 'mother_MC', 'child_MC',
           'father_overlap_coverage', 'mother_overlap_coverage', 'allele_ratio', 'child_ratio']
colformats = {'trid': str, 'denovo_coverage': 'UInt16', 'per_allele_reads_father': str,
                'per_allele_reads_mother': str, 'per_allele_reads_child': str,
                'index': 'UInt8', 'father_MC': str, 'mother_MC': str, 'child_MC': str,
                'father_overlap_coverage': str, 'mother_overlap_coverage': str,
                'allele_ratio': 'float32', 'child_ratio': 'float32'}
# Not currently used
# 'allele_origin', 'mean_diff_father', 'mean_diff_mother',
# 'genotype', 'allele_coverage', 'child_coverage', 'father_dropout_prob',
# 'mother_dropout_prob',  'denovo_status', 'father_dropout', 'mother_dropout',
# 'child_dropout', 'maxlh']

def get_chunks(files):
    """Read the next chunk for all of the files in the list.
    """
    all_sample_readers = []
    all_sample_ids = []
    # Open a reader for each file and store it in a list
    for file in files:
        reader = pd.read_csv(file, chunksize=1000000, sep='\t', usecols=usecols, dtype=colformats)
        all_sample_readers.append(reader)
        all_sample_ids.append(file.split('/')[-1].split('.')[0].split('_')[0])
    # Keep reading the next chunk for all the files until there are no more chunks
    while True:
        try:
            all_sample_chunks = []
            for sample, reader in zip(all_sample_ids, all_sample_readers):
                this_chunk = next(reader)
                all_sample_chunks.append(this_chunk)
                yield sample, this_chunk
            # Merge the chunks into a single DataFrame
            #merged = pd.concat(all_sample_chunks)
            #yield merged
        except StopIteration: #XXX Am I missing the end of the file here?
            break

def sum_depth(depths):
    """Sum the depth for each allele. Or any other string of comma separated numbers."""
    return [sum(pd.to_numeric(x.split(','))) for x in depths]

def min_depth(depths):
    """Sum the depth for each allele"""
    return [min(pd.to_numeric(x.split(','))) for x in depths]

def min_depth_trio(chunk):
    """Find the minimum depth accross trio for each locus in the chunk"""

    # Sum the depth for each individual then take the minimum accross the trio
    return chunk[['per_allele_reads_father',
          'per_allele_reads_mother',
          'per_allele_reads_child']].apply(sum_depth, axis=0).min(axis=1)

def min_allelic_depth(chunk):
    """Find the minimum allelic depth accross trio for each locus in the chunk"""

    # Take the minimum accross the trio
    return chunk[['per_allele_reads_father',
          'per_allele_reads_mother',
          'per_allele_reads_child']].apply(min_depth, axis=0).min(axis=1)

def parent_overlap_coverage(chunk):
    """Count the proportion of reads in the parents that support the denovo allele"""
    parent_overlap_coverage = chunk[['father_overlap_coverage', 'mother_overlap_coverage']].apply(sum_depth, axis=0).sum(axis=1)
    parent_coverage = chunk[['per_allele_reads_father', 'per_allele_reads_mother']].apply(sum_depth, axis=0).sum(axis=1)
    return parent_overlap_coverage / parent_coverage

min_depth_filtered = 0
complex_filtered = 0
def callable_loci(chunk, min_depth=10, min_allelic_depth=15):
    """Filter the chunk to keep only the loci that pass the filters"""
    global min_depth_filtered
    global complex_filtered
    
    min_depth_filtered += sum(chunk['min_allelic_depth'] < min_allelic_depth)
    complex_filtered += sum(chunk['min_motiflen'] != chunk['max_motiflen'])

    #print(chunk[chunk['child_MC'].str.contains('_') == True])
    #print(chunk[chunk['child_MC'].str.contains('_') == False])
    #print(chunk[chunk['min_motiflen'] == chunk['max_motiflen']])

    #print(f'Filtered out {min_depth_filtered} loci with depth < 10')
    #print(f'Filtered out {complex_filtered} loci with complex motifs')

    keep_loci = (
        #(chunk['min_depth'] >= min_depth) &
        (chunk['min_allelic_depth'] >= min_allelic_depth) &
        #(chunk['mother_dropout_prob'] < 0.005) &
        #(chunk['father_dropout_prob'] < 0.005) &
        #(chunk['mother_dropout'] == "N") &
        #(chunk['father_dropout'] == "N") &
        (chunk['min_motiflen'] == chunk['max_motiflen'])
        #(chunk['child_MC'].str.contains('_') == False)
        )
    return chunk[keep_loci]

def novel_allele(row):
    """Check if allele is absent from both parents
    
    Once we have phasing, this could be adapted.
    """
    kid_mc = row['child_MC'].split(',')[row['index']]
    if kid_mc in row['father_MC'].split(',') or kid_mc in row['mother_MC'].split(','):
        return False
    else:
        return True

def identify_dnms(chunk):
    """Identify the de novo mutations"""
    # Find the alleles that are not present in the parents
    pd.options.mode.chained_assignment = None
    chunk['candidate_dnm'] = chunk.apply(novel_allele, axis=1)
    pd.options.mode.chained_assignment = 'warn'
    is_dnm = (
        (chunk['candidate_dnm'] == True) &
        (chunk['denovo_coverage'] >= 2) #&
        # (chunk['father_overlap_coverage'] == '0,0') &
        # (chunk['mother_overlap_coverage'] == '0,0') &
        (chunk['parent_overlap_proportion'] <= 0.05) #& # proportion of reads supporting the denovo across both parents
        # (chunk['child_ratio'] >= 0.25) &
        # (chunk['child_ratio'] <= 0.75) &
        #(chunk['allele_ratio'] >= 0.5) # proportion of reads supporting that allele that support the denovo allele
        )

    return is_dnm

def count_dnms(chunk):
    """Count the number of de novo mutations in chunk"""

    n_dnms = identify_dnms(chunk).sum()
    n_callable = chunk.shape[0]

    return n_dnms, n_callable

def main(outfile: str, outdnms: str, samples: str,
         TRids: str = '',
         annotations: str = '../repeat-catalogs/human_GRCh38_no_alt_analysis_set.palladium-v1.0.trgt.annotations.bed',
         ):
    """
    Count the number of de novo mutations in each sample and write to a file.

    :param annotations: Path to the repeat annotation file.
    :param outfile: Path to the output file.
    :param outdnms: Path to the output file containing all DNMs.
    :param samples: Path to the TRGTdn files (comma separated).
    :param TRids: file of TRids to include (one per line).
    """
    #T2T:
    #chm13v2.0_maskedY_rCRS.palladium-v1.0.trgt.annotations.bed

    annot_df = pd.read_csv(annotations, sep='\t',
                            usecols=['TRid', 'longest_homopolymer', 'gc_content', 'n_motifs',
                                     'min_motiflen', 'max_motiflen', 'start', 'end'],
                            dtype={'TRid': str, 'longest_homopolymer': 'UInt16',
                                   'gc_content': 'float32', 'n_motifs': 'UInt16',
                                   'min_motiflen': 'UInt16', 'max_motiflen': 'UInt16',
                                   'start': 'UInt32', 'end': 'UInt32'})
    annot_df.rename(columns={'TRid': 'trid'}, inplace=True)

    if TRids:
        keepTRids = pd.read_csv(TRids, sep='\t', header=None, names=['trid'])['trid']

    # Create a generator to read the next chunk
    sample_files = samples.split(',')
    chunk_generator = get_chunks(sample_files)

    all_dnms = pd.DataFrame(columns=['Sample', 'DNMs', 'Callable', 'Group', 'Value'])
    #with open('out.txt', 'w') as outfile:
    #print('Sample\tDNMs\tCallable\tMotiflen\tDNMrate')

    all_dnms_list = []
    all_dnms_loci_list = []
    for sample, chunk in chunk_generator:

        # Filter to only the TRids of interest
        if TRids:
            #print(chunk['trid'], TRids['trid'])
            chunk = chunk[chunk['trid'].isin(keepTRids)]

        # Annotate the loci with the repeat annotations
        chunk = pd.merge(chunk, annot_df, on='trid', how='left')

        chunk['min_depth'] = min_depth_trio(chunk)
        chunk['min_allelic_depth'] = min_allelic_depth(chunk)
        chunk['parent_overlap_proportion'] = parent_overlap_coverage(chunk)
        
        chunk = callable_loci(chunk)

        chunk['reflen'] = chunk['end'] - chunk['start']

        # Save all the DNM rows
        all_dnms_loci_list.append(chunk[identify_dnms(chunk)])

        for motiflen, group in chunk.groupby('min_motiflen'):
            dnms, callable = count_dnms(group)
            all_dnms = all_dnms_list.append(pd.Series({'Sample': sample,
                                                    'DNMs': dnms, 'Callable': callable,
                                                    'Group': 'motiflen', 'Value': motiflen,
                                                    #'DNMrate': dnms/callable
                                                    }))
        
        for longest_homopolymer, group in chunk.groupby('longest_homopolymer'):
            dnms, callable = count_dnms(group)
            all_dnms = all_dnms_list.append(pd.Series({'Sample': sample,
                                                    'DNMs': dnms, 'Callable': callable,
                                                    'Group': 'longest_homopolymer',
                                                    'Value': longest_homopolymer,
                                                    #'DNMrate': dnms/callable
                                                    }))
        
        try:
            for reflen, group in chunk.groupby('reflen'):
                dnms, callable = count_dnms(group)
                all_dnms = all_dnms_list.append(pd.Series({'Sample': sample,
                                                        'DNMs': dnms, 'Callable': callable,
                                                        'Group': 'reflen', 'Value': reflen,
                                                        #'DNMrate': dnms/callable
                                                        }))
        except KeyError:
            sys.stderr.write(f'reflen not found in group code\n')
        
        try:
            for gc_content, group in chunk.groupby('gc_content'):
                dnms, callable = count_dnms(group)
                all_dnms = all_dnms_list.append(pd.Series({'Sample': sample,
                                                        'DNMs': dnms, 'Callable': callable,
                                                        'Group': 'gc_content', 'Value': gc_content,
                                                        #'DNMrate': dnms/callable
                                                        }))
        except KeyError:
            sys.stderr.write(f'gc_content not found in group code\n')
        

    all_dnms = pd.DataFrame(all_dnms_list)
    # merge by sample and group
    all_dnms = all_dnms.groupby(['Sample', 'Group', 'Value']).sum().reset_index()
    all_dnms['DNMrate'] = all_dnms['DNMs']/all_dnms['Callable']

    all_dnms.to_csv(outfile, sep='\t', index=False)

    all_dnms_loci = pd.concat(all_dnms_loci_list)
    all_dnms_loci.to_csv(outdnms, sep='\t', index=False)

    global min_depth_filtered
    global complex_filtered
    sys.stderr.write(f'Filtered out {min_depth_filtered} loci with insufficient depth\n')
    sys.stderr.write(f'Filtered out {complex_filtered} loci with complex motifs\n')

if __name__ == "__main__":
    # import doctest
    # doctest.testmod()
    import defopt
    defopt.run(main)
    