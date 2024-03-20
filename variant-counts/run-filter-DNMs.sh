# Run multiple instances of filter-DNMs.sh in parallel in the background


# Loci overlapping gymrek study

TRids="data/TRIDs_gymrek_matches_0.9_GRCh38.txt"

for f in data/TRGTdn/*.gz; do
    sample=$(basename $f .gz)
    python filter-DNMs.py $sample.gymrek.dnms.tsv $f $TRids &
done

# All loci

for f in data/TRGTdn/*.gz; do
    sample=$(basename $f .gz)
    python filter-DNMs.py $sample.dnms.tsv $f &
done
