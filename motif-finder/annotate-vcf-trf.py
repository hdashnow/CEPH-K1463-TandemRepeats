# Annotate the repeat motif structure of ref and/or alts in a VCF using Tandem Repeats Finder
# Calls out to https://github.com/lh3/TRF-mod/tree/dev

import subprocess
import sys
import cyvcf2
import pathlib
import typing as ty
import time
import numpy as np
start_time = time.time()

def main(vcf: pathlib.Path, output: str, trfmod: pathlib.Path = 'trf-mod',
         min:int = 25, max:int = 5000, ref: bool = False):
    """
    Takes a VCF file and runs Tandem Repeats Finder (via TRF-mod fork @lh3)
    on the REF and/or ALT sequences for each locus.
    For each repeat found it reports key statistics and the repeat motif.
    Note, with default settings minimum repeat called by TRF is 25 bp

    :param vcf: Path to the input VCF file.
    :param output: Prefix for output files, including temporary fasta. Make sure it's unique.
    :param trfmod: Path to the trf-mod executableble if it's not in the PATH.
    :param min: Minimum indel size to report (positive for insertions, negative for deletions)
    :param max: Maximum indel size to report (positive for insertions, negative for deletions)
    :param ref: Report the repeat structure of the reference allele.
    """
    vcfreader = cyvcf2.VCF(vcf)

    # Write alleles in fasta format for TRF-mod to read
    fastafile = f'{output}.fasta'
    with open(fastafile, 'w') as fasta:
        for vcf_line in vcfreader:

            info_string = format_vcf_info(vcf_line)
            
            if ref and len(vcf_line.REF) >= min and len(vcf_line.REF) <= max:
                header = f'{info_string};allele=REF'
                fasta.write(f'>{header}\n{vcf_line.REF}\n')
            for label, seq in zip([f'ALT.{i}' for i in range(len(vcf_line.ALT))], vcf_line.ALT):
                if len(seq) >= min and len(seq) <= max:
                    header = f'{info_string};allele={label}'
                    fasta.write(f'>{header}\n{seq}\n')

    # Run TRF-mod on the fasta file and save in an iterator
    trfmod_lines = run_trf(trfmod, fastafile)

    # Note, positions are relative to the start of the locus, not the chromosome
    trfmod_header = 'VCFid TRFstart TRFend TRFperiod TRFcopyNum TRFfracMatch TRFfracGap TRFscore TRFentroy TRFmotif'.split()

    # Iterate through the TRF-mod output
    trf_cols = ['TRFstart', 'TRFend', 'TRFperiod', 'TRFcopyNum', 'TRFfracMatch', 'VCFid', 'TRFmotif']
    vcf_cols1 = ['CHROM', 'POS', 'reflen', 'altlen']
    header =  vcf_cols1 + [ 'allele_type'] + trf_cols 

    with open(f'{output}.tsv', 'w') as outfh:
        outfh.write('#' + '\t'.join(header) + '\n')

        for trf_line in trfmod_lines:
            # Parse the TRF-mod line
            # if len(trf_line) == 0:
            #     continue
            trfmod_fields = trf_line.split()
            info_dict = parse_vcf_info(trfmod_fields[0].strip('>'))

            trfmod_dict = dict(zip(trfmod_header, trfmod_fields))
            trfmod_dict['VCFid'] = info_dict['ID']

            # What type of allele is this?
            trfmod_allele = info_dict['allele']

            if trfmod_allele.startswith('ALT'):
                allele_type = 'ALT'
                allele_pos = int(trfmod_allele.lstrip('ALT.'))
                altlen = info_dict['altlens'][allele_pos]

            elif trfmod_allele == 'REF':
                allele_type = 'REF'
                altlen = None

            else:
                raise ValueError(f'Unknown allele type: {trfmod_allele}')
            
            try:
                vcf_chrom = info_dict['CHROM']
                vcf_pos = info_dict['POS']
                vcf_reflen = info_dict['reflen']
            except KeyError:
                print(info_dict)
                raise

            # Write output as a tab-separated file
            outfh.write('\t'.join([ str(x) for x in
                [vcf_chrom, vcf_pos, vcf_reflen] +
                [altlen, allele_type] +
                [trfmod_dict[field] for field in trf_cols]
            ]) + '\n')

def format_vcf_info(vcf_line: cyvcf2.Variant) -> str:
    """
    Format the VCF info field as a string to be written to the fasta header
    format: key1=value1;key2=value2;...
    e.g. CHROM=chr1;POS=10033;reflen=1;altlen=2;ID=chr1-10034-INS-1
    """
    chrom = vcf_line.CHROM
    pos = vcf_line.POS
    if vcf_line.ID is None:
        locid = '-'.join(str(x) for x in [vcf_line.CHROM, vcf_line.POS]) # might need to be more specific
    else:
        locid = vcf_line.ID
    if ';' in locid:
        locid = locid.replace(';', ',')
    if '=' in locid:
        locid = locid.replace('=', ':')
    reflen = len(vcf_line.REF)
    altlens = ','.join([str(len(alt)) for alt in vcf_line.ALT])
    outstring = f'CHROM={chrom};POS={pos};reflen={reflen};altlens={altlens};ID={locid}'
    return outstring

def parse_vcf_info(info_string: str) -> dict:
    """
    >>> parse_vcf_info('CHROM=chr1;POS=10033;reflen=1;altlens=2;ID=chr1-10034-INS-1')
    {'CHROM': 'chr1', 'POS': 10033, 'reflen': 1, 'altlens': [2], 'ID': 'chr1-10034-INS-1'}
    """
    info_dict = dict()
    for keyval in info_string.split(';'):
        try:
            key, val = keyval.split('=')
        except ValueError:
            print(info_string)
            raise
        if key == 'altlens':
            val = [int(x) for x in val.split(',')]
        info_dict[key] = val
        if key == 'POS':
            info_dict[key] = int(val)
        if key == 'reflen':
            info_dict[key] = int(val)
    return info_dict

def norm_chrom(chrom: str) -> str:
    """
    Normalize the chromosome name. This may fail for species with more than 90 autosomes or Z/W.
    >>> norm_chrom('chr1')
    '01'
    >>> norm_chrom('chrX')
    '91'
    >>> norm_chrom('chr2')
    '02'
    """
    chrom = chrom.lstrip('chr')
    if chrom == 'X':
        chrom = '91'
    elif chrom == 'Y':
        chrom = '92'
    elif chrom == 'M':
        chrom = '93'
    chrom = chrom.zfill(2)
    return chrom


# Genomic coordinates object
class GenomicPos(ty.NamedTuple):
    chrom: str
    pos: int

    # Initialize from a cyvcf2 record
    @classmethod
    def from_cyvcf(cls, record: cyvcf2.Variant) -> 'GenomicPos':
        return cls(record.CHROM, record.POS)
    # Initialize from chrom and pos
    @classmethod
    def from_chrom_pos(cls, chrom: str, pos: int) -> 'GenomicPos':
        return cls(chrom, pos)
        
    # Print as a string
    def __str__(self):
        return f'{self.chrom}:{self.pos}'

    # Compare two genomic positions
    def __lt__(self, other):
        """
        >>> GenomicPos('chr22', 200) < GenomicPos('chrX', 100)
        True
        >>> GenomicPos('chr1', 100) < GenomicPos('chr1', 1000)
        True
        >>> GenomicPos('chr9', 1) < GenomicPos('chr10', 1)
        True
        """

        if self.chrom == other.chrom:
            return self.pos < other.pos
        else:
            # Add leading zeros to the chromosome name
            # To make sure that chr10 comes after chr9 and not before chr1

            # Remove leading 'chr' if present
            self_chrom_num = norm_chrom(self.chrom)
            other_chrom_num = norm_chrom(other.chrom)

            return self_chrom_num < other_chrom_num

    # Test if equal
    def __eq__(self, other):
        return self.chrom == other.chrom and self.pos == other.pos

def run_trf(trfmod, fastafile):
    trfmod = pathlib.Path(trfmod).resolve()
    fasta = pathlib.Path(fastafile).resolve()
    trfmod_command = f'{trfmod} {fasta}'.split()

    proc = subprocess.Popen(trfmod_command, stdout=subprocess.PIPE)

    while True:
        output = proc.stdout.readline()

        if output == b'' and proc.poll() is not None:
            break
        if output:
            if len(output) > 0:
                yield output.decode('utf-8').strip()

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    import defopt
    defopt.run(main)

    elapsed = time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))
    print(f'Elapsed time: {elapsed}', file=sys.stderr)
