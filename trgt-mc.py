import cyvcf2
import pathlib
import typing as ty

def distance(a: ty.List[int], b: ty.List[int], pow=1) -> int:
    """
    Calculate the distance between two alleles.
    pow of 1 is manhattan distance, pow of 2 is euclidean distance.
    """
    return (abs(a[0] - b[0])**pow + abs(a[1] - b[1])**pow)**(1/pow)

def mc(v: cyvcf2.Variant) -> ty.List[int]:
    m = [int(x) for x in v.format("MC")[0].split(',')]
    return m

def find_distance(mom: cyvcf2.Variant, dad: cyvcf2.Variant, kid: cyvcf2.Variant, pow: int) -> ty.Tuple[int, ty.List[int]]:
    mom_mc = mc(mom)
    dad_mc = mc(dad)
    kid_mc = mc(kid)
    parent_mcs = [[mom_mc[0], dad_mc[0]], [mom_mc[0], dad_mc[1]], [mom_mc[1], dad_mc[0]], [mom_mc[1], dad_mc[1]]]
    dists = [(distance(parent_mc, kid_mc, pow=pow), parent_mc) for parent_mc in parent_mcs]
    return min(dists, key=lambda x: x[0])

def vmc_fmt(variant: cyvcf2.Variant) -> str:
    """
    Return the MC format field as a string.
    """
    return f'{",".join(str(x) for x in mc(variant))}'

def mc_fmt(mc: ty.List[int]) -> str:
    """
    Return the MC format field as a string.
    """
    return f'{",".join(str(x) for x in mc)}'

def main(mom_vcf: pathlib.Path, dad_vcf: pathlib.Path, kid_vcfs: ty.List[pathlib.Path], *, pow: int = 1):
    """
    :param mom_vcf: Path to the maternal TRGT VCF file.
    :param dad_vcf: Path to the paternal TRGT VCF file.
    :param kid_vcfs: List of paths to one or more child TRGT VCF files.
    :param pow: Power to use for distance calculation. 1 is manhattan distance, 2 is euclidean distance.
    """
    vcfs = [cyvcf2.VCF(mom_vcf), cyvcf2.VCF(dad_vcf)]
    vcfs.extend([cyvcf2.VCF(kid_vcf) for kid_vcf in kid_vcfs])
    kid_ids = [vcf.samples[0] for vcf in vcfs[2:]]
    print(f"#variant\tkid_id\tmom_mc\tdad_mc\tkid_mc\tmc_dist{pow}\tparent_ht")
    for vs in zip(*vcfs):
        mom, dad = vs[:2]
        assert mom.POS == dad.POS and mom.CHROM == dad.CHROM \
                and mom.REF == dad.REF
        for i, kid in enumerate(vs[2:]):
            assert mom.POS == kid.POS and mom.CHROM == kid.CHROM \
                    and mom.REF == kid.REF
            d, parental_ht = find_distance(mom, dad, kid, pow)
            print(f'{mom.CHROM}:{mom.POS}:{mom.REF}:{"".join(mom.ALT) or "."}\t{kid_ids[i]}\t{vmc_fmt(mom)}\t{vmc_fmt(dad)}\t{vmc_fmt(kid)}\t{d}\t{mc_fmt(parental_ht)}')






if __name__ == "__main__":
    import defopt
    defopt.run(main)
