#!/usr/bin/env python

import sys
import cyvcf2
import pathlib
import subprocess
import re

def extractvcf(infile: pathlib.Path, minsize: int = 8):
    cyvcf2_vcf = cyvcf2.VCF(infile)
    for locus in cyvcf2_vcf:
        if locus.INFO['SVTYPE'] == 'INS' and locus.INFO['SVLEN'] >= minsize:
            yield (locus.CHROM, locus.POS, locus.POS), locus.ALT

def extractfasta(bed: pathlib.Path, fasta: pathlib.Path):
    """
    Extract sequences from a FASTA file using a BED file.
    :param bed: Path to the BED file.
    :param fasta: Path to the FASTA file.
    :return: DNA sequence as a string.
    """
    for locus in bed:
        yield (locus.chrom, locus.start, locus.end), [fasta.fetch(locus.chrom, locus.start, locus.end)]

def runtrsolve(sequence: str, trsolve: pathlib.Path):
    """
    Run tr-solve on a sequence by calling the binary.
    Example command:
    echo input | tr-solve > output

    :param sequence: DNA sequence as a string.
    :return: List of TRGT motifs.


    >>> runtrsolve('ATATATATATATATATATATATAT', trsolve='~/tools/tr-solve-v0.2.0-linux_x86_64')
    ['AT']
    >>> runtrsolve('CAGCAGCAGCAGCAGCAGAAAAAAAAAAAAAAAAAAAA', trsolve='~/tools/tr-solve-v0.2.0-linux_x86_64')
    ['CAG', 'A']
    """
    command = f'echo {sequence} | {trsolve}'
    result = subprocess.run(command, shell=True, check=True, 
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    motifs = result.stdout.decode('utf-8').strip().split('\t')
    
    try:
        return [re.sub("[\(].*?[\)]", "", x) for x in motifs[2].split(',')]
    except IndexError:
        return []

def main(infile: pathlib.Path, outfile: str = 'stdout', *, 
         fasta: pathlib.Path = None, minsize: int = 8,
         trtools: pathlib.Path = pathlib.Path('tr-solve')):
    """
    Convert a VCF or BED file to TRGT repeat definitions.
    :param infile: Path to the input VCF or BED file.
    :param outfile: Path to the output TRGT repeat definitions file.
    :param fasta: Path to the reference FASTA file (required for BED inputs).
    :param minsize: Minimum insertion size to use from VCF.
    :param trtools: Path to the tr-solve binary.
    """
    if infile.suffix == '.vcf' or infile.suffixes[-2:] == ['.vcf', '.gz']:
        sequences = extractvcf(infile, minsize=minsize)
    elif infile.suffix == '.bed':
        sequences = extractfasta(infile, fasta)
    else:
        raise ValueError(f'Unknown input file type: {infile.suffix}')
    
    if outfile == 'stdout':
        f = sys.stdout
    else:
        f = open(outfile, 'w')

    for (chrom, start, end), alts in sequences:
        for alt in alts:
            motifs = runtrsolve(alt, trsolve=trtools)
            
            if len(motifs) == 0:
                f.write(f'{alt}\tNone\n')
                continue
            unique_motifs = dict.fromkeys(motifs) # can use set(motifs) if order doesn't matter
            motifs_str = ','.join(unique_motifs)
            struc = ''.join([f'({motif})n' for motif in motifs])
            # XXX Delete the input sequence printing once tested - remove {alt}\t below and f.write(f'{alt}\tNone\n') above
            f.write(f'{alt}\t{chrom}\t{start}\t{end}\tID={chrom}_{start}_{end};MOTIFS={motifs_str};STRUC={struc}\n')

    f.flush()

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    import defopt
    defopt.run(main)