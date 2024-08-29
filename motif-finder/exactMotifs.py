import pathlib
import cyvcf2
import re

def main(vcf: pathlib.Path, output: str):
    """
    Takes a TRGT VCF file and counts exact motifs
    on all REF and ALT sequences. For each repeat found it reports if it matches
    the expected TRGT motif.

    :param vcf: Path to the input VCF file.
    :param output: Prefix for output files, including temporary fasta. Make sure it's unique.
    :param trfmod: Path to the tr-solve executableble if it's not in the PATH. 
    """

    vcfreader = cyvcf2.VCF(vcf)
    samples = vcfreader.samples

    # Detect motif structure in VCF alleles and write results to a tsv file
    tsvfile = f'{output}.tsv'
    with open(tsvfile, 'w') as tsvout:
        header = ['locid', 'TRID', 'alleleID', 'sample', 'motif', 'count', 'longest', 'sequence']
        tsvout.write('\t'.join(header) + '\n')
        for vcf_line in vcfreader:
            locid = '-'.join(str(x) for x in [vcf_line.CHROM, vcf_line.POS])
            alleles = [] # format is a list of alleles where each is a list of dicts with motif, count, longest
            for alleleID, seq in zip(['REF'] + [f'ALT.{i+1}' for i in range(len(vcf_line.ALT))], [vcf_line.REF] + vcf_line.ALT):
                seq = seq.upper()
                this_allele = []
                
                for motif in vcf_line.INFO['MOTIFS'].split(','):
                    motif = motif.upper()
                    count = count_exact_motifs(seq, motif)
                    longest = count_longest_exact_motif(seq, motif)
                    this_allele.append({'count': count, 'longest': longest, 'seq': seq, 'alleleID': alleleID, 'motif': motif})
                alleles.append(this_allele)
            for sample, sample_GT in zip(samples, vcf_line.genotypes):
                if len(sample_GT) != 3:
                    continue
                for GT in sample_GT[:2]:
                    allele = alleles[GT] # Get the allele corresponding to the genotype
                    for motif_dict in allele:
                        tsvout.write(f'{locid}\t{vcf_line.INFO['TRID']}\t{motif_dict['alleleID']}\t{sample}\t{motif_dict['motif']}\t{motif_dict['count']}\t{motif_dict['longest']}\t{motif_dict['seq']}\n')

def count_exact_motifs(sequence: str, motif: str):
    """
    Count the number of exact motifs in a sequence.

    :param sequence: The sequence to search for motifs.
    :param motif: The motif to search for.
    :return: The number of exact motifs found.

    >>> count_exact_motifs('CAGCAGCAG', 'CAG')
    3
    >>> count_exact_motifs('CAGCAGCAGCGGCAG', 'CAG')
    4
    >>> count_exact_motifs('GGCGGCGGCGGCGGCAGCGGCGGCTGCGGCGGCGGCGGCGGCAGCCATATATGCCA', 'GCN')
    16
    """
    if 'N' in motif:
        return len(re.findall(motif.replace('N', '.'), sequence))
    else:
        return sequence.count(motif)

def count_longest_exact_motif(sequence: str, motif: str):
    """
    Count the number of uninterrupted exact motifs in a sequence.

    :param sequence: The sequence to search for motifs.
    :param motif: The motif to search for.
    :return: The length of the longest exact motif found.

    >>> count_longest_exact_motif('CAGCAGCAGCGGCAG', 'CAG')
    3
    >>> count_longest_exact_motif('AAAATAAAATAAAATAAAATAAAATAAAATAAAATAAATAAA', 'GAAAT')
    0
    >>> count_longest_exact_motif('AAAATAAAATAAAATAAAATAAAATAAAATAAAATAAATAAA', 'AAAAT')
    7
    >>> count_longest_exact_motif('AAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAATAAAAAAATAAAATAAAATAAAAT', 'AAAAT')
    7
    >>> count_longest_exact_motif('GGCGGCGGCGGCGGCAGCGGCGGCTGCGGCGGCGGCGGCGGCAGCCATATATGCC', 'GCN')
    15
    """
    count = 0
    longest = 0

    if 'N' in motif:
        # replace N with . in motif
        pattern = re.compile(motif.replace('N', '.'))
        splits = re.split(pattern, ('^' + sequence + '$'))
    else:
        splits = ('^' + sequence + '$').split(motif)
    if len(splits) == 1:
        return 0
    for i in splits:
        if i == '':
            count += 1
            longest = max(longest, count)
        else:
            count = 0
    return longest + 1

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    import defopt
    defopt.run(main)