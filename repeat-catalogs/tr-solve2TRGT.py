#!/usr/bin/env python

import sys
import cyvcf2
import pathlib
import subprocess
import re
import pysam
from itertools import groupby

def extractvcf(infile: pathlib.Path, minins: int = 8):
    cyvcf2_vcf = cyvcf2.VCF(infile)
    for locus in cyvcf2_vcf:
        if locus.INFO['SVTYPE'] == 'INS' and locus.INFO['SVLEN'] >= minins:
            yield (locus.CHROM, locus.POS, locus.POS), locus.ALT

def extractfasta(bed: pathlib.Path, fasta: pathlib.Path, maxout: int = 10000):
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
            if int(locus[2]) - int(locus[1]) > maxout:
                continue
            try:
                sequence = ref.fetch(locus[0], int(locus[1]), int(locus[2]))
            except ValueError as e:
                sys.stderr.write(f'Error: {e}\n')
                continue
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
    try:
        result = subprocess.run(command, shell=True, check=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f'Error: {command}\n')
        sys.stderr.write(f'{e.stderr.decode("utf-8")}\n')
        return [], [None, None] #XXX: Should this be an error?
    motif_string = result.stdout.decode('utf-8').strip().split('\t')

    if len(motif_string) < 3: # No motif structure reported
        return [], [None, None]
    
    # Parse the motifs and their start/end positions
    # e.g. 'CAG(0-18),A(18-38)' -> ['CAG', 'A'], [0, 38]
    try:
        re_motifs = re.compile(r'([A-Z]+)\(([0-9]+)-([0-9]+)\)').findall(motif_string[2])
    except IndexError as e:
        sys.stderr.write(f'Error: Cannot parse {motif_string} - {e}\n')
        raise
    minpos = None
    maxpos = None
    motifs = []
    for motif, start, end in re_motifs:
        assert 'N' not in motif, f'N found in motif: {motif}'

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

def rmdup(seq):
    """remove adjacent duplicates from a list

    >>> rmdup([1, 1, 2, 3, 3, 3, 4, 5, 5])
    [1, 2, 3, 4, 5]
    >>> rmdup(['a', 'a', 'b', 'c', 'c', 'c', 'a', 'd', 'e', 'e'])
    ['a', 'b', 'c', 'a', 'd', 'e']
    """
    return [x[0] for x in groupby(seq)]

def main(infile: pathlib.Path, outfile: str = 'stdout', *, 
         fasta: pathlib.Path = None, minins: int = 8, maxout: int = 10000,
         trtools: pathlib.Path = pathlib.Path('tr-solve'),
         seqout: bool = False):
    """
    Convert a VCF or BED file to TRGT repeat definitions.
    :param infile: Path to the input VCF or BED file.
    :param outfile: Path to the output TRGT repeat definitions file.
    :param fasta: Path to the reference FASTA file (required for BED inputs).
    :param minins: Minimum insertion size to use from VCF.
    :param maxout: Maximum region size to output (smaller ones skipped).
    :param trtools: Path to the tr-solve binary.
    :param seqout: Output all input sequences. When false, sequences with no motif found are skipped.
    """
    if infile.suffix == '.vcf' or infile.suffixes[-2:] == ['.vcf', '.gz']:
        sequences = extractvcf(infile, minins=minins)
    elif infile.suffix == '.bed':
        if fasta is None:
            raise ValueError('BED input requires FASTA file')
        sequences = extractfasta(infile, fasta, maxout=maxout)
    else:
        raise ValueError(f'Unknown input file type: {infile.suffix}')
    
    if outfile == 'stdout':
        f = sys.stdout
    else:
        f = open(outfile, 'w')

    for (chrom, start_orig, end_orig), alts in sequences:
        for alt in alts:
            motifs, bounds = runtrsolve(alt, trsolve=trtools)
            
            if len(motifs) == 0:
                if seqout:
                    f.write(f'{alt}\tNone\n')
                continue
            start = int(start_orig)
            end = int(end_orig)

            # Trim the bounds to the bases covered by the motifs
            if bounds[0] is not None:
                start += bounds[0]
            if bounds[1] is not None:
                end = start + bounds[1] - bounds[0]

            if start > int(end_orig):
                start = int(end_orig)
            if end < int(start_orig):
                end = int(start_orig)

            unique_motifs = dict.fromkeys(motifs) # can use set(motifs) if order doesn't matter. Assume it does for now
            motifs_str = ','.join(unique_motifs)
            struc = ''.join([f'({motif})n' for motif in rmdup(motifs)])
            if seqout:
                outstring = f'{alt}\t'
            else:
                outstring = ''
            outstring += f'{chrom}\t{start}\t{end}\tID={chrom}_{start}_{end};MOTIFS={motifs_str};STRUC={struc}\n'
            f.write(outstring)

    f.flush()

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    import defopt
    defopt.run(main)