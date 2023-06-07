import sys
import cyvcf2
import pathlib
import typing as ty
import numpy as np
import plotly.express as px


def distance(a: ty.List[int], b: ty.List[int], pow=1) -> int:
    """
    Calculate the distance between two alleles.
    pow of 1 is manhattan distance, pow of 2 is euclidean distance.
    """
    return (abs(a[0] - b[0])**pow + abs(a[1] - b[1])**pow)**(1/pow)

def mc(v: cyvcf2.Variant, tag: str) -> ty.List[int]:
    mc = v.format(tag)[0]
    if isinstance(mc, np.ndarray):
        return mc.tolist()
    if mc == '.':
        return [-1]
    return [int(x) for x in mc.split(',')]

def find_distance(mom: cyvcf2.Variant, dad: cyvcf2.Variant, kid: cyvcf2.Variant,
                  pow: int, tag: str) -> ty.Tuple[int, ty.List[int]]:
    mom_mc = mc(mom, tag)
    dad_mc = mc(dad, tag)
    kid_mc = mc(kid, tag)
    parent_mcs = [[mom_mc[0], dad_mc[0]], [mom_mc[0], dad_mc[1]], [mom_mc[1], dad_mc[0]], [mom_mc[1], dad_mc[1]]]
    dists = [(distance(parent_mc, kid_mc, pow=pow), parent_mc) for parent_mc in parent_mcs]
    return min(dists, key=lambda x: x[0])

def vmc_fmt(variant: cyvcf2.Variant, tag: str) -> str:
    """
    Return the tag format field as a string.
    """
    return f'{",".join(str(x) for x in mc(variant, tag))}'

def mc_fmt(mc: ty.List[int]) -> str:
    """
    Return the tag format field as a string.
    """
    return f'{",".join(str(x) for x in mc)}'

def main(mom_vcf: pathlib.Path, dad_vcf: pathlib.Path, kid_vcfs:
         ty.List[pathlib.Path], *, pow: int = 1, output_prefix: str = "trgt-mc-",
         tag: str = "MC", exclude_chroms: list[str] = ['chrX', 'chrY']):
    """
    :param mom_vcf: Path to the maternal TRGT VCF file.
    :param dad_vcf: Path to the paternal TRGT VCF file.
    :param kid_vcfs: List of paths to one or more child TRGT VCF files.
    :param pow: Power to use for distance calculation. 1 is manhattan distance, 2 is euclidean distance.
    :param output_prefix: prefix for output files.
    :param tag: VCF format field used for lengths.
    :param exclude_chroms: chromosomes to exclude.
    """
    vcfs = [cyvcf2.VCF(mom_vcf), cyvcf2.VCF(dad_vcf)]
    vcfs.extend([cyvcf2.VCF(kid_vcf) for kid_vcf in kid_vcfs])
    kid_ids = [vcf.samples[0] for vcf in vcfs[2:]]
    fh = open(output_prefix + 'dists.txt', 'w')
    print(f"#variant\tkid_id\tmom_{tag}\tdad_{tag}\tkid_{tag}\t{tag}_dist{pow}\tparent_ht",
          file=fh)
    dists = []
    for vs in zip(*vcfs):
        mom, dad = vs[:2]
        if mom.CHROM in exclude_chroms: continue
        assert mom.POS == dad.POS and mom.CHROM == dad.CHROM \
                and mom.REF == dad.REF
        for i, kid in enumerate(vs[2:]):
            assert mom.POS == kid.POS and mom.CHROM == kid.CHROM \
                    and mom.REF == kid.REF
            try:
                d, parental_ht = find_distance(mom, dad, kid, pow, tag)
            except IndexError:
                d = -1
                parental_ht = [-1, -1]
            dists.append(d)
            print(f'{mom.CHROM}:{mom.POS}:{mom.REF}:{"".join(mom.ALT) or "."}\t{kid_ids[i]}\t{vmc_fmt(mom, tag)}\t{vmc_fmt(dad, tag)}\t{vmc_fmt(kid, tag)}\t{d}\t{mc_fmt(parental_ht)}', file=fh)

    fig = px.histogram(x=dists)
    m = ["manhattan", "euclidean"][pow - 1]
    fig.update_layout(xaxis_range=[-1, 25], xaxis_title=f"{m}-distance between child and parent alleles ({tag})")
    fig.write_html(output_prefix + 'dists.html')
    print(f'wrote {output_prefix}dists.html and {output_prefix}dists.txt',
          file=sys.stderr)


if __name__ == "__main__":
    import defopt
    defopt.run(main)
