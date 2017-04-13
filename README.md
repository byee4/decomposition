installs may require some stuff:

conda install mkl mkl-service numpy pandas sklearn bokeh matplotlib

# Examples!

### usage (to plot log2(rpkm) from featurecounts.txt, given a conditions file, saving principle components to file) :
```
python runner.py -l2 \
-f \
-rpkm \
-i featureCounts.txt \
-o pca.png \
-c conditions.tsv \
-cc treatment \
-pc
```

### usage (to plot a generic matrix) given a conditions file :
```
python runner.py -i examples/iris.txt \
-o examples/iris.png \
-c examples/iris.names \
-cc names
```

### usage (to plot a generic matrix) given a conditions file (keeping intermediates):
```
python runner.py -i examples/alll_clip_alone_kmer_statistics.txt \
-o examples/kmer.png \
-c examples/alll_clip_alone_kmer_statistics_conditions.txt \
-cc condition \
-k
```