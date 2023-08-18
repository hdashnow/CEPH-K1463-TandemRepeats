# Annotate the repeat motif structure of ref and alts in a VCF using Tandem Repeats Finder
# Calls out to https://github.com/lh3/TRF-mod/tree/dev

import subprocess
import sys
import cyvcf2
import pathlib
import typing as ty
import time
import numpy as np
start_time = time.time()

def main(vcf: pathlib.Path, output: str, trfmod: pathlib.Path = 'trf-mod'):
    """
    Takes a TRGT VCF file and runs Tandem Repeats Finder (via TRF-mod fork @lh3)
    on all REF and ALT sequences. For each repeat found it reports if it matches
    the expected TRGT motif.

    :param vcf: Path to the input VCF file.
    :param output: Prefix for output files, including temporary fasta. Make sure it's unique.
    :param trfmod: Path to the trf-mod executableble if it's not in the PATH.
    """
    vcfreader = cyvcf2.VCF(vcf)

    # Write alleles in fasta format for TRF-mod to read
    fastafile = f'{output}.fasta'
    with open(fastafile, 'w') as fasta:
        for vcf_line in vcfreader:
            locid = '-'.join(str(x) for x in [vcf_line.CHROM, vcf_line.POS])
            for label, seq in zip(['REF'] + [f'ALT.{i}' for i in range(len(vcf_line.ALT))], [vcf_line.REF] + vcf_line.ALT):
                fasta.write(f'>{locid}-{label}\n{seq}\n')

    # Run TRF-mod on the fasta file and save in an iterator
    trfmod_lines = run_trf(trfmod, fastafile)

    # Note, positions are relative to the start of the locus, not the chromosome
    trfmod_header = 'VCFid TRFstart TRFend TRFperiod TRFcopyNum TRFfracMatch TRFfracGap TRFscore TRFentroy TRFmotif'.split()


    # Match up the TRF-mod output with the VCF
    # Iterate through the VCF and the TRF-mod output in parallel. Need to keep track of 
    # the current locus as there are multiple lines per locus in the TRF-mod output
    trf_cols = ['TRFstart', 'TRFend', 'TRFperiod', 'TRFcopyNum', 'TRFfracMatch', 'VCFid', 'TRFmotif']
    vcf_cols1 = ['CHROM', 'POS', 'TRGTmc']
    vcf_cols2 = ['TRGTmotif', 'TRGTstruc']
    header =  vcf_cols1 + [ 'allele_type', 'match'] + trf_cols + vcf_cols2
    with open(f'{output}.tsv', 'w') as outfh:
        outfh.write('#' + '\t'.join(header) + '\n')

        vcfreader = cyvcf2.VCF(vcf)

        # create a new vcf Writer using the input vcf as a template.
        vcfreader.add_info_to_header({'ID': 'REFmc',
                                      'Description': 'Reference allele motif count calculated by TandemRepeatsFinder', 
                                      'Type': 'Character', 'Number': '1'})
        vcfreader.add_format_to_header({'ID': 'ALTmc',
                                        'Description': 'Alternate allele motif count calculated by TandemRepeatsFinder', 
                                        'Type': 'Character', 'Number': '1'})
        outvcf = cyvcf2.Writer(f'{output}.vcf', vcfreader)

        vcf_line = next(vcfreader) # Read the first line of the VCF
        trf_line = next(trfmod_lines) # Read the first line of the TRF-mod output

        while True:
            # Parse the VCF line
            vcf_locus = GenomicPos.from_cyvcf(vcf_line)
            # Parse the TRF-mod line
            # if len(trf_line) == 0:
            #     continue
            trfmod_fields = trf_line.split()
            trfmod_dict = dict(zip(trfmod_header, trfmod_fields))
            trfmod_locus = GenomicPos.from_chrom_pos(trfmod_dict['VCFid'].split('-')[0],
                                                     int(trfmod_dict['VCFid'].split('-')[1]))
            
            # If they match, write the output
            if trfmod_locus == vcf_locus:
                # What type of allele is this?
                trfmod_allele = trfmod_dict['VCFid'].split('-')[2]
                vcf_motifs = vcf_line.INFO.get('MOTIFS')

                if vcf_motifs == trfmod_dict['TRFmotif']:
                    match = True
                else:
                    match = False

                if trfmod_allele.startswith('ALT'):
                    allele_type = 'ALT'
                    allele_pos = int(trfmod_allele.lstrip('ALT.'))
                    mc = int(vcf_line.format('MC')[0].split(',')[allele_pos])
                    if match:
                        vcf_line.set_format('ALTmc', str(np.array([trfmod_dict['TRFcopyNum']])))
                    else:
                        vcf_line.set_format('ALTmc', '0')

                elif trfmod_allele == 'REF':
                    allele_type = 'REF'
                    mc = None
                    match = None
                    if match:
                        vcf_line.INFO['REFmc'] = trfmod_dict['TRFcopyNum']
                    else:
                        vcf_line.INFO['REFmc'] = '0'


                # Write output as a tab-separated file
                outfh.write('\t'.join([ str(x) for x in
                    [vcf_line.CHROM, vcf_line.POS, mc] +
                    [allele_type, match] +
                    [trfmod_dict[field] for field in trf_cols] +
                    [vcf_motifs, vcf_line.INFO.get('STRUC')]
                ]) + '\n')

                try:
                    trf_line = next(trfmod_lines)
                    continue

                except StopIteration:
                    break

            # Choose which file to read next
            # If the TRF-mod locus is before the VCF locus, read the next TRF-mod line
            if trfmod_locus < vcf_locus:
                try:
                    trf_line = next(trfmod_lines)
                except StopIteration:
                    break

            else:
                # Write the VCF line to the output VCF (since we're finished with this VCF line at this point)
                outvcf.write_record(vcf_line)

                # If the TRF-mod locus is after the VCF locus, read the next VCF line
                if trfmod_locus > vcf_locus:
                    try:
                        vcf_line = next(vcfreader)
                    except StopIteration:
                        break
                # If the TRF-mod locus is the same as the VCF locus, read both
                elif trfmod_locus == vcf_locus:
                    try:
                        trf_line = next(trfmod_lines)
                        vcf_line = next(vcfreader)
                    except StopIteration:
                        break

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
    print(f'Elapsed time: {elapsed}')
