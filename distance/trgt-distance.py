import sys
import cyvcf2
import pathlib
import typing as ty
import numpy as np
import plotly.express as px
from itertools import cycle


def distance(a: ty.List[int], b: ty.List[int], power=1) -> int:
    """
    Calculate the distance between two alleles.
    power of 1 is manhattan distance, power of 2 is euclidean distance.
    """
    return sum(abs(a[i] - b[i])**power for i in range(len(a)))**(1/power)

def get_tag(v: cyvcf2.Variant, tag: str) -> ty.List[int]:
    if tag == "GT":
        return "/".join(map(str, v.genotypes[0][:-1]))

    val = v.format(tag)[0]
    if isinstance(val, np.ndarray):
        return val.tolist()
    if val == '.':
        return [np.nan]
    return [int(x) for x in val.split(',')]

def get_norm_tag(v: cyvcf2.Variant, tag: str) -> ty.List[int]:
    factors = [1] # Normalization factor (to get values into bp space)
    if tag == 'MC':
        motifs = v.INFO.get('MOTIFS').split(',')
        factors = [len(x) for x in motifs]
    if tag == 'AP':
        factors = get_tag(v, 'AL')

    vals = get_tag(v, tag)

    for i in vals:
        if type(i) != str:
            if i < 0:
                raise ValueError(f'{v.POS} {tag} is negative: {i}') 
            if np.isnan(i):
                raise ValueError # These are expected missing values

    return [x*f for x, f in zip(vals, cycle(factors))]

def generate_combinations(mom: ty.List[ty.List[int]], dad: ty.List[ty.List[int]]) -> ty.List[ty.List[int]]:
    """
    Generate all possible combinations of alleles from mom and dad.
    >>> generate_combinations([[1, 0], [2, 3]], [[0, 1], [3, 4]])
    [[1, 0, 2, 3], [1, 1, 2, 4], [0, 0, 3, 3], [0, 1, 3, 4]]

    Note that when we take the 0th item, we must take the 0th from all lists for that parent.
    """
    result = []
    for i in [0, 1]:
        for j in [0, 1]:
            c = []
            for k in range(len(mom)):
                c.extend([mom[k][i], dad[k][j]])
            result.append(c)
    return result


def flip(arr):
    """
    >>> flip([0, 2, 0, 2])
    [2, 0, 2, 0]
    """
    result = arr[:]
    result[::2] = arr[1::2]
    result[1::2] = arr[::2]
    return result


def min_distance(mom_mc: ty.List[ty.List[int]], dad_mc: ty.List[ty.List[int]],
                kid_mc: ty.List[int], power: int) -> ty.Tuple[int, ty.List[int]]:
    """
    >>> min_distance([[6021, 15812], [6032, 15805], [5828.002823531628, 15366.006873607635]], [[3874, 4398], [3886, 4408], [3730.9991244077682, 4270.000632762909]], [4398, 15808, 4408, 15805, 4270.000632762909, 15317.999767303467], power = 1)
    (52.0071063041687, [15812, 4398, 15805, 4408, 15366.006873607635, 4270.000632762909])
    """
    parent_mcs = generate_combinations(mom_mc, dad_mc)
    dists = [(distance(parent_mc, kid_mc, power=power), parent_mc) for parent_mc in parent_mcs]
    result_md = min(dists, key=lambda x: x[0])

    parent_mcs = generate_combinations(dad_mc, mom_mc)
    dists = [(distance(parent_mc, kid_mc, power=power), parent_mc) for parent_mc in parent_mcs]
    result_dm = min(dists, key=lambda x: x[0])
    # flip order of dm (dad/mom) to match md (mom/dad)
    result_dm = (result_dm[0], flip(result_dm[1]))

    result = min([result_md, result_dm], key=lambda x: x[0])
    return result

def find_distance(mom: cyvcf2.Variant, dad: cyvcf2.Variant, kid: cyvcf2.Variant,
                  power: int, tags: ty.List[str]) -> ty.Tuple[int, ty.List[int]]: #, ty.List[int]]:
    mom_mc = []
    dad_mc = []
    kid_mc = []
    for tag in tags:
        try:
            mom_mc.append(get_norm_tag(mom, tag))
            dad_mc.append(get_norm_tag(dad, tag))
            kid_mc.extend(get_norm_tag(kid, tag))
        except ValueError:
            raise ValueError

    return min_distance(mom_mc, dad_mc, kid_mc, power = power)

def vmc_fmt(variant: cyvcf2.Variant, tag: str) -> str:
    """
    Return the tag format field as a string.
    """
    if tag == 'ALLR' or tag == 'MS':
        return variant.format(tag)[0]
    try:
        t = get_norm_tag(variant, tag)
    except ValueError:
        t = [np.nan]
    if isinstance(t, str):
        return t
    if tag == 'GT':
        return ''.join(t)
    return f'{",".join(str(x) for x in t)}'

def mc_fmt(mc: ty.List[int]) -> str:
    """
    Return the tag format field as a string.
    """
    return f'{",".join(str(x) for x in mc)}'

def main(mom_vcf: pathlib.Path, dad_vcf: pathlib.Path, kid_vcfs:
         ty.List[pathlib.Path], *, power: int = 1, output_prefix: str = "trgt-mc-",
         dist_tags: ty.List[str] = ["MC", "AL", "AP"], extra_tags: ty.List[str] = ["SD"],
         exclude_chroms: ty.List[str] = ['chrX', 'chrY']):
    """
    :param mom_vcf: Path to the maternal TRGT VCF file.
    :param dad_vcf: Path to the paternal TRGT VCF file.
    :param kid_vcfs: List of paths to one or more child TRGT VCF files.
    :param power: Power to use for distance calculation. 1 is manhattan distance, 2 is euclidean distance.
    :param output_prefix: prefix for output files.
    :param dist_tags: VCF format field used for lengths.
    :param exclude_chroms: chromosomes to exclude.
    """
    vcfs = [cyvcf2.VCF(mom_vcf), cyvcf2.VCF(dad_vcf)]
    vcfs.extend([cyvcf2.VCF(kid_vcf) for kid_vcf in kid_vcfs])
    kid_ids = [vcf.samples[0] for vcf in vcfs[2:]]
    fh = open(output_prefix + 'dists.txt', 'w')
    header = "#variant\tkid_id\t" \
        + "\t".join(f"{x}_{tag}" for x in ["mom", "dad", "kid"] for tag in dist_tags)  \
        + "\t" \
        + "\t".join(f"{x}_{tag}" for x in ["mom", "dad", "kid"] for tag in extra_tags) \
        + "\tmotifs" \
        + "\tdist\tparent_ht"
    print(header, file=fh)
    dists = []
    for vs in zip(*vcfs):
        mom, dad = vs[:2]
        if mom.CHROM in exclude_chroms: continue
        if not (mom.INFO.get('TRID') == dad.INFO.get('TRID')):
            sys.stderr.write(f'Warning: Skipped locus due to parent TRID mismatch. Mom: {mom.INFO.get("TRID")}, Dad: {dad.INFO.get("TRID")}\n'.format())
            continue
        if not (mom.POS == dad.POS and mom.CHROM == dad.CHROM and mom.REF == dad.REF):
            sys.stderr.write(f'Warning: Parent postition mismatch. Mom: {mom.CHROM} {mom.POS}, Dad: {dad.CHROM} {dad.POS}\n'.format())
        for i, kid in enumerate(vs[2:]):
            if not (mom.INFO.get('TRID') == kid.INFO.get('TRID')):
                sys.stderr.write(f'Warning: Skipped kid {kid_ids[i]} at this locus due to kid TRID mismatch. Parents: {mom.INFO.get("TRID")}, Kid: {kid.INFO.get("TRID")}\n'.format())
                continue
            if not (mom.POS == kid.POS and mom.CHROM == kid.CHROM and kid.REF == kid.REF):
                sys.stderr.write(f'Warning: Kid postition mismatch. Parents: {mom.CHROM} {mom.POS}, Kid: {kid.CHROM} {kid.POS}\n'.format())
            try:
                d, parental_ht = find_distance(mom, dad, kid, power, dist_tags)
            except ValueError:
                d = -1
                parental_ht = [-1, -1]
            dists.append(d)
            motifs = kid.INFO.get('MOTIFS')

            line = f'{mom.INFO.get("TRID")}\t{kid_ids[i]}\t'
            line += "\t".join(f"{vmc_fmt(x, tag)}" for x in [mom, dad, kid] for tag in dist_tags)
            line += "\t"
            line += "\t".join(f"{vmc_fmt(x, tag)}" for x in [mom, dad, kid] for tag in extra_tags)
            line += f"\t{motifs}"
            line += f"\t{d}\t{','.join(str(p) for p in parental_ht)}"
            print(line, file=fh)

    fig = px.histogram(x=dists)
    m = ["manhattan", "euclidean"][power - 1]
    fig.update_layout(xaxis_range=[-1, 25], xaxis_title=f"{m}-distance between child and parent alleles ({','.join(dist_tags)})")
    fig.write_html(output_prefix + 'dists.html')
    print(f'wrote {output_prefix}dists.html and {output_prefix}dists.txt',
          file=sys.stderr)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    import defopt
    defopt.run(main)
