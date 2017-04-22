installs may require some stuff:

conda install mkl mkl-service numpy pandas sckit-learn bokeh matplotlib

# Examples!

### usage (to plot log2(rpkm) from featurecounts.txt, given a conditions file, saving principle components to file) :
```
python runner.py -l2 \
-f \
-rpkm \
-i examples/data/counts.txt \
-o examples/data/pca.png \
-c examples/data/conditions.txt \
-cc response \
-k
```

### usage (to plot a generic matrix) given a conditions file :
```
python runner.py -i examples/data/iris.txt \
-o examples/data/iris.png \
-c examples/data/iris.names \
-cc names
```

### usage (to plot a generic matrix) given a conditions file (keeping intermediates):
```
python runner.py -i examples/data/alll_clip_alone_kmer_statistics.txt \
-o examples/data/kmer.png \
-c examples/data/alll_clip_alone_kmer_statistics_conditions.txt \
-cc condition \
-k
```