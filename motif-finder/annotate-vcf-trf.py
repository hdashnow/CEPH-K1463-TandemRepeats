# Annotate the repeat motif structure of ref and alts in a VCF using Tandem Repeats Finder
# Calls out to https://github.com/lh3/TRF-mod/tree/dev

import subprocess
import sys
import cyvcf2
import pathlib
import typing as ty

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
        for v in vcfreader:
            locid = '-'.join(str(x) for x in [v.CHROM, v.POS])
            for label, seq in zip(['REF'] + [f'ALT.{i}' for i in range(len(v.ALT))], [v.REF] + v.ALT):
                fasta.write(f'>{locid}-{label}\n{seq}\n')

    # Run TRF-mod on the fasta file and write output to trfmod_out
    trfmod = pathlib.Path(trfmod).resolve()
    fasta = pathlib.Path(fastafile).resolve()
    trfmod_command = f'{trfmod} {fasta}'.split()        
    #print(''.join(trfmod_command))
    result = subprocess.run(trfmod_command, capture_output=True, text=True)

    # Note, positions are relative to the start of the locus, not the chromosome
    trfmod_header = 'VCFid start end period copyNum fracMatch fracGap score entroy TRFmotif'.split()
    # with open(f'{output}.trfmod', 'w') as trfmod_out:
    #     trfmod_out.write('\t'.join(trfmod_header) + '\n')
    #     trfmod_out.write(result.stdout)

    # Match up the TRF-mod output with the VCF
    # Iterate through the VCF and the TRF-mod output in parallel. Need to keep track of 
    # the current locus as there are multiple lines per locus in the TRF-mod output
    trfmod_cols = ['VCFid', 'start', 'end', 'period', 'copyNum',  'TRFmotif']
    header = ['CHROM', 'POS', 'TRGTmotifs', 'TRGTstruc'] + trfmod_cols + ['match']
    with open(f'{output}.tsv', 'w') as outfh:
        outfh.write('#' + '\t'.join(header) + '\n')

        vcfreader = cyvcf2.VCF(vcf)
        v = next(vcfreader)
        vcf_locus = v.CHROM + '-' + str(v.POS)

        trfmod_lines = iter(result.stdout.split('\n'))
        for trfmod_line in trfmod_lines:
            # Parse the TRF-mod line
            if len(trfmod_line) == 0:
                continue
            trfmod_fields = trfmod_line.split()
            trfmod_dict = dict(zip(trfmod_header, trfmod_fields))
            trfmod_locus = '-'.join(trfmod_dict['VCFid'].split('-')[0:2])
            
            # If we've moved on to a new locus, read the next line from TRF-mod
            while trfmod_locus == vcf_locus:
                vcf_motifs = v.INFO.get('MOTIFS')
                if vcf_motifs == trfmod_dict['TRFmotif']:
                    match = True
                else:
                    match = False
                # Write output as a tab-separated file
                outfh.write('\t'.join([ str(x) for x in
                    [v.CHROM, v.POS, 
                    vcf_motifs, v.INFO.get('STRUC')] +
                    [trfmod_dict[field] for field in trfmod_cols] + [match]
                ]) + '\n')
                
                try:
                    v = next(vcfreader)
                    vcf_locus = v.CHROM + '-' + str(v.POS)
                except StopIteration:
                    break

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    import defopt
    defopt.run(main)
