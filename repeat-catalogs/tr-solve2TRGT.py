#!/usr/bin/env python

import sys
import cyvcf2
import pathlib
import subprocess
import re
import pysam
from itertools import groupby
from TRlib import extractfasta

def extractvcf(infile: pathlib.Path, minins: int = 8):
    cyvcf2_vcf = cyvcf2.VCF(infile)
    for locus in cyvcf2_vcf:
        if locus.INFO['SVTYPE'] == 'INS' and locus.INFO['SVLEN'] >= minins:
            yield (locus.CHROM, locus.POS, locus.POS), locus.ALT

def parsetrsolve(motif_string: str):
    """
    Parse the output of tr-solve version 0.2.1+
    Parse the motifs and their start/end positions
    e.g. 'CAG_0_18,A_18_38' -> ['CAG', 'A'], [0, 38]
    
    >>> parsetrsolve('AGAGGCGCGGCGCGCCGGCGCAGGCGCAG_0_597')
    ('AGAGGCGCGGCGCGCCGGCGCAGGCGCAG', 0, 597)
    >>> parsetrsolve('CAG_0_18')
    ('CAG', 0, 18)
    >>> parsetrsolve('A_18_38')
    ('A', 18, 38)
    """

    motif, start, end = motif_string.split('_')
        # except ValueError as e:
        #     sys.stderr.write(f'Error: Cannot parse {motif_list} - {e}\n')
        #     raise
    assert 'N' not in motif, f'N found in motif: {motif}'
    start = int(start)
    end = int(end)

    return motif, start, end

def runtrsolve(sequence: str, trsolve: pathlib.Path):
    """
    Run tr-solve on a sequence by calling the binary.
    Example command:
    echo input | tr-solve > output
    Report the bases of the TRGT motifs.

    DNA sequence as a string.
    Return a list of TRGT motifs and the start/end of the left and rightmost motifs.


    >>> runtrsolve('ATATATATATATATATATATATAT', trsolve='~/tools/tr-solve-v0.2.1-linux_x86_64')
    (['AT'], [0, 24])
    >>> runtrsolve('CAGCAGCAGCAGCAGCAGAAAAAAAAAAAAAAAAAAAA', trsolve='~/tools/tr-solve-v0.2.1-linux_x86_64')
    (['CAG', 'A'], [0, 38])
    >>> runtrsolve('GGCACGGCATATATATATATATATATATATAT', trsolve='~/tools/tr-solve-v0.2.1-linux_x86_64')
    (['AT'], [8, 32])
    """
    sequence = sequence.upper()
    if 'N' in sequence:
        sys.stderr.write(f'N found in sequence: {sequence}\n')
        # Remove Ns from the sequence
        sequence = re.sub('N', '', sequence)
        sys.stderr.write(f'Removed Ns from sequence to produce: {sequence}\n')
    command = f'echo {sequence} | {trsolve}'
    #sys.stderr.write(f'Running: {command}\n')
    try:
        result = subprocess.run(command, shell=True, check=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f'Error: {command}\n')
        sys.stderr.write(f'{e.stderr.decode("utf-8")}\n')
        return [], [None, None] #XXX: Should this be an error?
    
    motif_string_list = result.stdout.decode('utf-8').strip().split('\t')
    if len(motif_string_list) < 3: # No motif structure reported
        return [], [None, None]
    
    # Trim the positions to the bases covered by the motifs
    minpos = None
    maxpos = None
    motifs = []

    for motif_string in motif_string_list[2].split(','):
        motif, start, end = parsetrsolve(motif_string)

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
         seqout: bool = False, trim: bool = False):
    """
    Convert a VCF or BED file to TRGT repeat definitions.
    :param infile: Path to the input VCF or BED file.
    :param outfile: Path to the output TRGT repeat definitions file.
    :param fasta: Path to the reference FASTA file (required for BED inputs).
    :param minins: Minimum insertion size to use from VCF.
    :param maxout: Maximum region size to output (smaller ones skipped).
    :param trtools: Path to the tr-solve binary.
    :param seqout: Output all input sequences. When false, sequences with no motif found are skipped.
    :param trim: Trim the output coordinates to the bases covered by the motifs.
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
            #sys.stderr.write(f'{motifs}\t{bounds}\n')
            
            if len(motifs) == 0:
                if seqout:
                    f.write(f'{alt}\tNone\n')
                continue
            start = int(start_orig)
            end = int(end_orig)

            if trim:

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