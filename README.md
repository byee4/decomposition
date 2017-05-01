installs may require some stuff:

```
conda create -n decomposition python=2.7
source activate decomposition
conda install mkl mkl-service numpy pandas scikit-learn bokeh matplotlib
```

```
git clone https://github.com/byee4/decomposition
cd decomposition
python setup.py build
python setup.py install
```
# Examples!

### usage (to plot log2(rpkm) from [featurecounts](http://bioinf.wehi.edu.au/featureCounts/) counts.txt, given a conditions file (--conditions), saving principle components to file) :

Example [conditions.txt](https://github.com/byee4/decomposition/blob/master/examples/data/conditions.txt)

```bash
decompose \
--input examples/data/counts.txt \
--output examples/data/pca.png \
--conditions examples/data/conditions.txt \
--algorithm PCA \
-cc response \
-f \
-l2 \
-rpkm \
-k
```

### usage (to plot a generic matrix) given a conditions file :
```bash
decompose -i examples/data/iris.txt \
-o examples/data/iris.png \
-c examples/data/iris.names \
-cc names
```

### usage (to plot a generic matrix) given a conditions file (keeping intermediates):
```bash
decompose -i examples/data/alll_clip_alone_kmer_statistics.txt \
-o examples/data/kmer.png \
-c examples/data/alll_clip_alone_kmer_statistics_conditions.txt \
-cc condition \
-k
```
