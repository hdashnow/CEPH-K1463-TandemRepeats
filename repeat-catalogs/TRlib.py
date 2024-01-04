import pysam
import pathlib
import sys
import re

def extractfasta(bed: pathlib.Path, fasta: pathlib.Path, maxout: int = 10000, returninfo: bool = False):
    """
    Extract sequences from a FASTA file using a BED file.
    :param bed: Path to the BED file.
    :param fasta: Path to the FASTA file.
    :param maxout: Maximum region size to output (smaller ones skipped).
    :param returninfo: Return the locus information as well as the sequence.
    :return: DNA sequence as a string.
    """
    # read fasta file
    ref = pysam.FastaFile(fasta)

    with open(bed) as f:
        for line in f:
            locus = line.strip().split()
            if maxout is not None and int(locus[2]) - int(locus[1]) > maxout:
                continue
            try:
                sequence = ref.fetch(locus[0], int(locus[1]), int(locus[2]))
            except KeyError as e:
                sys.stderr.write(f'Error: {e}\n')
                continue
            if returninfo:
                yield (locus[0], locus[1], locus[2], locus[3]), sequence
            else:
                yield (locus[0], locus[1], locus[2]), [sequence]

def longesthomopolymer(sequence):
    """
    Find the longest homopolymer in a sequence. 
    If there are multiple longest homopolymers, the last one is returned.
    Split on Ns
    :param sequence: DNA sequence as a string.
    :return: length of the longest homopolymer.

    >>> longesthomopolymer('ATTTTGGGGGCCCCC')
    5
    >>> longesthomopolymer('ATTTTGGGGGNNNNNNNCCCCC')
    5
    >>> longesthomopolymer('ATTNTTGGGNGGNNNNNNNCCCNCC')
    3
    >>> longesthomopolymer('') is None
    True
    >>> longesthomopolymer('NNNNNNNNNNNNNN') is None
    True
    """
    # if sequence is all Ns, return None
    if len(sequence) == sequence.count('N'):
        return None
    maximum = count = 0
    for s in filter(None, re.split('N', sequence, flags=re.IGNORECASE)):
        current = ''
        for c in s:
            if c == current:
                count += 1
            else:
                count = 1
                current = c
            maximum = max(count,maximum)
    return maximum

def GCcontent(sequence):
    """
    Calculate the GC content of a sequence.
    :param sequence: DNA sequence as a string.
    :return: GC content as a float.

    >>> GCcontent('ATTTTGGGGGCCCCCATTTT')
    0.5
    >>> GCcontent('ATTT')
    0.0
    >>> GCcontent('ATTTTGGGGGNNNNNNNNNNNCCCCCATTTT')
    0.5
    >>> GCcontent('GGGGCCC')
    1.0
    >>> GCcontent('') is None
    True
    >>> GCcontent('NNNNNNNNNNN') is None
    True
    """
    if len(sequence) == sequence.count('N'):
        return None
    sequence = sequence.upper()
    return (sequence.count('G') + sequence.count('C')) / (len(sequence) - sequence.count('N'))

if __name__ == "__main__":
    import doctest
    doctest.testmod()