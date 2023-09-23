#!/usr/bin/env python

import sys
import cyvcf2
import pathlib
import subprocess
import re
import pysam

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
    # read fasta file
    ref = pysam.FastaFile(fasta)

    with open(bed) as f:
        for line in f:
            locus = line.strip().split()
            sequence = ref.fetch(locus[0], int(locus[1]), int(locus[2]))
            yield (locus[0], locus[1], locus[2]), [sequence]

def runtrsolve(sequence: str, trsolve: pathlib.Path):
    """
    Run tr-solve on a sequence by calling the binary.
    Example command:
    echo input | tr-solve > output
    Report the bases of the TRGT motifs.

    DNA sequence as a string.
    Return a list of TRGT motifs and the start/end of the left and rightmost motifs.


    >>> runtrsolve('ATATATATATATATATATATATAT', trsolve='~/tools/tr-solve-v0.2.0-linux_x86_64')
    (['AT'], [0, 24])
    >>> runtrsolve('CAGCAGCAGCAGCAGCAGAAAAAAAAAAAAAAAAAAAA', trsolve='~/tools/tr-solve-v0.2.0-linux_x86_64')
    (['CAG', 'A'], [0, 38])
    >>> runtrsolve('GGCACGGCATATATATATATATATATATATAT', trsolve='~/tools/tr-solve-v0.2.0-linux_x86_64')
    (['AT'], [8, 32])
    """
    command = f'echo {sequence} | {trsolve}'
    result = subprocess.run(command, shell=True, check=True, 
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    motif_string = result.stdout.decode('utf-8').strip().split('\t')

    if len(motif_string) == 1:
        return [], [None, None]
    
    # Parse the motifs and their start/end positions
    # e.g. 'CAG(0-18),A(18-38)' -> ['CAG', 'A'], [0, 38]
    re_motifs = re.compile(r'([A-Z]+)\(([0-9]+)-([0-9]+)\)').findall(motif_string[2])
    minpos = None
    maxpos = None
    motifs = []
    for motif, start, end in re_motifs:
        motifs.append(motif)
        if minpos is None:
            minpos = int(start)
        if maxpos is None:
            maxpos = int(end)
        if int(start) < minpos:
            minpos = int(start)
        if int(end) > maxpos:
            maxpos = int(end)

    return motifs, [minpos, maxpos]

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
        if fasta is None:
            raise ValueError('BED input requires FASTA file')
        sequences = extractfasta(infile, fasta)
    else:
        raise ValueError(f'Unknown input file type: {infile.suffix}')
    
    if outfile == 'stdout':
        f = sys.stdout
    else:
        f = open(outfile, 'w')

    for (chrom, start, end), alts in sequences:
        for alt in alts:
            motifs, bounds = runtrsolve(alt, trsolve=trtools)
            
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