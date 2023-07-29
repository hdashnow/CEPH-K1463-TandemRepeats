
mom=mom.mini.vcf
dad=dad.mini.vcf
kid=kid.mini.vcf


python ../trgt-distance.py \
	${mom} \
	${dad} \
	-k ${kid} \
	--dist-tags AL MC AP \
	--extra-tags SD GT \
	--output-prefix trgt.mini.
