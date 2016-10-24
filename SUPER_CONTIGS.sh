

python3 igor_mandric.py assembly_contigs.fa reference.fa IGOR

/home/igorm/MUMmer3.23/nucmer -p nucmerhit reference.fa IGOR.fa --maxmatch
/home/igorm/MUMmer3.23/delta-filter -i 97 -l 180 -r nucmerhit.delta > nucmerhit.filter
/home/igorm/MUMmer3.23/show-coords -dTlro nucmerhit.filter > nucmerhit.coords

python build_scaffolds.py > perfect.repeat.scaf


