# distances of allele lengths in trios for STRs from TRGT

```
D=/path/to/vcfs/TRGT/

mom=200100_M_adotto_v02.sorted.vcf.gz
dad=2189_SF_adotto_v02.sorted.vcf.gz
kid1=200103
kid2=200104

python trgt-distance.py \
	$D/${mom} \
	$D/${dad} \
	-k $D/${kid1}*.vcf.gz $D/${kid2}*.vcf.gz \
	--dist-tags AL MC AP \
	--extra-tags SD \
	--output-prefix trgt.
```

will write: `trgt.dists.html` and `trgt.dists.txt`

