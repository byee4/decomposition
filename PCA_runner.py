import matplotlib

matplotlib.use('Agg')
import sys
import os
import numpy as np
import seaborn as sns
import pandas as pd
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import PCAPlotter

DEBUG = 0
TESTRUN = 0
PROFILE = 0

__all__ = []
__version__ = 0.1
__date__ = '2015-12-19'
__updated__ = '2015-12-19'

def counts_to_rpkm(df):
    """
    calculates and returns counts normalized by gene length and seq depth
    :param df: featureCounts counts.txt file
    :return: matrix of rpkms
    """
    counts = df.ix[:, 5:]
    lengths = df['Length']
    mapped_reads = counts.sum()
    return (counts * pow(10, 9)).div(mapped_reads, axis=1).div(lengths, axis=0)


def rgb_to_hex(rgb):
    """
    Returns list of html hex codes for each rgb code in rgb list
    :param rgb: list of rgb
    :return: list of hex
    """
    hexcolors = ['#%02x%02x%02x' % (c[0] * 255, c[1] * 255, c[2] * 255) for c in rgb]  # rgb to hex
    return hexcolors


def expr_to_hex(expr, cmap='Purples', is_norm=True):
    """
    From a list of expression values, return a list of html hex values
    :param expr: list of expression values
    :param cmap: colormapping
    :param is_norm: normalize expression 0-1
    :return: list of hex
    """
    cmap = plt.get_cmap(cmap)

    if is_norm:
        norm = [float(i) / max(expr) for i in expr]  # normalize expr val [0-1]
    else:
        norm = expr

    rgbs = [cmap(color) for color in norm]  # add color
    return rgb_to_hex(rgbs)

def expr_to_rgb(expr, cmap='Purples', is_norm=True):
    cmap = plt.get_cmap(cmap)

    if is_norm:
        norm = [float(i) / max(expr) for i in expr]  # normalize expr val [0-1]
    else:
        norm = expr

    rgbs = [cmap(color) for color in norm]  # add color
    return rgbs

def color_by_condition(df, col_string):
    """
    Takes df of conditions (rows of samples, cols of conditions)
    and returns a list of colors matching each distinct condition.

    :param df: pandas DataFrame
    :param col_string: name of column specifying the condition
    :return: list of colors associated with each distinct condition
    """
    max_conditions = set(df[col_string])
    # colors = sns.color_palette("hls", (len(max_conditions)))
    colormap = {}
    c = 1
    for condition in max_conditions:
        colormap[condition] = {
            'color':c, # colors[c],
            'condition':condition
        }
        c+=1
    colormapped = df[col_string].apply(lambda x: colormap[x])
    return dict(colormapped)

def subset(df, subset_file):
    """
    returns a dataframe subset given a file containing a list of indices
    :param df: dataframe of gene expression values
    :param subset_file: file of genes to subset, line delimited
    :return: dataframe
    """
    genes = [gene.strip() for gene in open(subset_file, 'r')]
    return df.loc[genes, :].dropna(axis=0)

def main(argv=None):  # IGNORE:C0111
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)

    # Setup argument parser
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--output", dest="output", required=True,
                        help='output pdf/svg/png (will also output PC components as [prefix].txt')
    parser.add_argument("-i", "--input", dest="input", required=True,
                        help="input matrix as featureCounts or matrix")
    parser.add_argument("-rpkm", "--rpkm", dest="rpkm", default=False, action='store_true',
                        help="flag converts expression counts to rpkm (requires featureCounts)")
    parser.add_argument("-l2", "--log2", dest="log2", default=False, action='store_true',
                        help="flag log2-transforms expression counts")
    parser.add_argument("-f", "--featureCounts", dest="featureCounts", default=False, action='store_true',
                        help="True if featurecounts-format")
    parser.add_argument("-V", "--version", action='version', version=program_version_message)
    parser.add_argument("-sc", "--sum_cutoff", dest="cutoff", type=int, default=0)
    parser.add_argument("-g", "--gene", dest="gene", default=None, type=str)
    parser.add_argument("-s", "--subset", dest="subset", default=None)
    parser.add_argument("-c", "--conditions", dest="conditions", default=None)
    parser.add_argument("-cc", "--conditions-col", dest="conditions_col", default=None)
    parser.add_argument("-pc", "--save-pc", dest="savepc", action='store_true', default=True)
    # Process arguments
    args = parser.parse_args()

    counts_file = args.input
    output_file = args.output
    subset_file = args.subset
    conditions_file = args.conditions
    conditions_col = args.conditions_col

    is_featureCounts = args.featureCounts
    is_rpkm = args.rpkm
    is_log2 = args.log2

    sum_cutoff = args.cutoff
    gene_id = args.gene

    """ read in counts file """
    counts = pd.read_table(counts_file, index_col=0, comment='#')


    """ featureCounts """
    if is_featureCounts:
        if is_rpkm:
            counts = counts_to_rpkm(counts)
        else:
            counts = counts.ix[:, 5:]

    """ colors """
    colors = {}
    cmap = None

    if gene_id in counts.index: # color by single gene expression
        cmap = 'Purples'
        expr = counts.loc[gene_id]
        if is_log2:
            expr = np.log2(expr+1)
        for key, value in expr.iteritems():
            colors[key] = {'color': value / max(expr), 'condition': 'expression'}

    elif conditions_file and os.path.exists(conditions_file): # color by condition
        flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e",
                  "#2ecc71"] # hardcoded for now..
        sns.set_palette(flatui)
        cmap = ListedColormap(flatui)
        conditions_df = pd.read_table(
            conditions_file,
            index_col=0
        )
        colors = color_by_condition(conditions_df, conditions_col)
    else: #
        for col in counts.columns:
            colors[col] = {'color': 'blue', 'condition': 'condition'}

    colors = pd.DataFrame(colors).T
    colors.to_csv('examples/colors.txt',sep='\t')
    """ do pca on select genes only """
    if subset_file and os.path.exists(subset_file):
        counts = subset(counts, subset_file)


    """ removes rows whos sum (reads or rpkm) < cutoff """
    if sum_cutoff > 0:
        counts = counts[counts.sum(axis=1) > sum_cutoff]

    """ log2 """
    if is_log2:
        counts = np.log2(counts+1)

    fig, ax = plt.subplots()
    plotter = PCAPlotter.pcaplot(counts, cmap, colors=colors, ax=ax, bokeh=False)
    leg = plt.legend(loc='best', shadow=False, frameon = 1)
    leg.get_frame().set_edgecolor('b')
    leg.get_frame().set_facecolor('w')
    fig.savefig(output_file)

    plotter.prcomp.to_csv(os.path.splitext(output_file)[0] + '_.txt', sep='\t')

if __name__ == "__main__":
    if DEBUG:
        sys.argv.append("-h")
        sys.argv.append("-v")
        sys.argv.append("-r")
    if TESTRUN:
        import doctest

        doctest.testmod()
    if PROFILE:
        import cProfile
        import pstats

        profile_filename = 'profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())