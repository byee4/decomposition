import matplotlib

matplotlib.use('Agg')
import sys
import os
import logging

import matplotlib.pyplot as plt

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

import PCAPlotter
import TSNEPlotter
import ICAPlotter
import KMeansPlotter
import color_helpers as ch
from ClusterExperiment import ClusterExperiment

DEBUG = 0
TESTRUN = 0
PROFILE = 0

SEP = "\t"

__all__ = []
__version__ = 0.2
__date__ = '2015-12-19'
__updated__ = '2017-2-13'

def main(argv=None):  # IGNORE:C0111
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (
        program_version,
        program_build_date
    )

    # Setup argument parser
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--output",
                        dest="output",
                        required=True,
                        help='output pdf/svg/png (will also output PC " + \
                        "components as [prefix].txt')
    parser.add_argument("-i", "--input",
                        dest="input",
                        required=True,
                        help="input matrix as featureCounts or matrix")
    parser.add_argument("-rpkm", "--rpkm",
                        dest="rpkm",
                        default=False,
                        action='store_true',
                        help="flag converts expression counts to rpkm " + \
                             "(requires featureCounts)")
    parser.add_argument("-l2", "--log2",
                        dest="log2",
                        default=False,
                        action='store_true',
                        help="flag log2-transforms expression counts")
    parser.add_argument("-f", "--featureCounts",
                        dest="featureCounts",
                        default=False,
                        action='store_true',
                        help="True if featurecounts-format")
    parser.add_argument("-V", "--version",
                        action='version',
                        version=program_version_message)
    parser.add_argument("-sc", "--sum_cutoff",
                        dest="cutoff",
                        type=int,
                        default=0,
                        help="any more than this level (reads) will not " + \
                             "be counted.")
    parser.add_argument("-s", "--subset",
                        dest="subset",
                        default=None,
                        help="line-delimited file containing the genes of " + \
                             "interest")
    parser.add_argument("-g", "--gene",
                        dest="gene_id",
                        default=None,
                        type=str,
                        help="gene id of expression values by which to " + \
                             "color the PCA.")
    parser.add_argument("-c", "--conditions",
                        dest="conditions",
                        default=None,
                        help="table of (row:samples,col:conditions) that " + \
                             "define sample conditions")
    parser.add_argument("-cc", "--conditions-col",
                        dest="conditions_col",
                        default=None,
                        help="column within the conditions file " + \
                             "(-cc flag must be on)" + \
                             "by which to categorize the plot legend.")
    parser.add_argument("-k", "--keep-intermediates",
                        dest="keep",
                        default=False,
                        action='store_true',
                        help="True if we want to keep all intermediates")
    parser.add_argument("-a", "--algorithm",
                        dest="algorithm",
                        default='PCA',
                        type=str,
                        help="Algorithm ([PCA] by default, or 'tSNE')")

    # Process arguments
    args = parser.parse_args()

    # io
    counts_file = args.input
    output_file = args.output
    subset_file = args.subset
    conditions_file = args.conditions
    conditions_col = args.conditions_col
    algorithm = args.algorithm.upper()

    is_featurecounts = args.featureCounts
    is_rpkm = args.rpkm
    is_log2 = args.log2
    keep_intermediates = args.keep
    sum_cutoff = args.cutoff
    gene_id = args.gene_id

    # prefix
    prefix = os.path.splitext(output_file)[0]

    # Process logging info
    logger = logging.getLogger('PCA_runner')
    logger.setLevel(logging.INFO)
    ih = logging.FileHandler(prefix + ".log")
    eh = logging.FileHandler(prefix + ".err")
    ih.setLevel(logging.INFO)
    eh.setLevel(logging.ERROR)
    logger.addHandler(ih)
    logger.addHandler(eh)
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ih.setFormatter(formatter)
    eh.setFormatter(formatter)
    logger.info("starting program")

    """ read in counts file """
    logger.info(sys.argv)
    experiment = ClusterExperiment(
        counts_file=counts_file,
        conditions_file=conditions_file,
        conditions_col=conditions_col,
        gene_id=gene_id,
        is_featurecounts=is_featurecounts,
    )

    """ do pca on select genes only """
    if subset_file and os.path.exists(subset_file):
        logger.info("SUBSET on: {}".format(subset_file))
        logger.info(
            "SUBSET SIZE (before): {}".format(
                experiment.counts.data.shape[0]
            )
        )
        experiment.counts.subset(subset_file)
        if keep_intermediates:
            experiment.counts.data.to_csv(prefix + ".subset.txt",sep=SEP)
        logger.info(
            "SUBSET SIZE (after): {}".format(
                experiment.counts.data.shape[0]
            )
        )

    """ rpkm """
    if is_rpkm:
        logger.info("RPKM FLAG ON")
        experiment.counts.as_rpkm()
        if keep_intermediates:
            experiment.counts.data.to_csv(prefix + ".rpkm.txt", sep=SEP)

    """ removes rows whos sum (reads) < cutoff """
    if sum_cutoff > 0:
        logger.info("CUTOFF AT {} READ ROW SUMS".format(sum_cutoff))
        logger.info(
            "CUTOFF SIZE (before): {}".format(
                experiment.counts.data.shape[0]
            )
        )
        experiment.counts.min_row_sum_cutoff(sum_cutoff)
        if keep_intermediates:
            experiment.counts.data.to_csv(prefix + ".cutoff.txt", sep=SEP)
        logger.info(
            "CUTOFF SIZE (before): {}".format(
                experiment.counts.data.shape[0]
            )
        )

    """ log2 """
    if is_log2:
        logger.info("LOG2 FLAG ON")
        experiment.counts.as_log2(1)
        if keep_intermediates:
            experiment.counts.data.to_csv(prefix + ".log2.txt", sep=SEP)

    """ save metadata """
    if keep_intermediates:
        experiment.metadata.to_csv(prefix + ".metadata.txt", sep=SEP)

    """ get appropriate cmap """
    if conditions_file is not None and conditions_col is not None:
        cmap = ch.hex_to_cmap(experiment.metadata.shape[1])
    else:
        cmap = 'Purples'

    """ plot stuff """
    fig, ax = plt.subplots()

    if algorithm == 'PCA':
        plotter = PCAPlotter.pcaplot(
            experiment,
            cmap,
            ax=ax, bokeh=False)
        plotter.prcomp.to_csv(prefix + '.pcacomp.txt', sep=SEP)
    elif algorithm == 'TSNE':
        plotter = TSNEPlotter.tsneplot(
            experiment,
            cmap,
            ax=ax, bokeh=False)
        plotter.tcomp.to_csv(prefix + '.tsnecomp.txt', sep=SEP)
    elif algorithm == 'ICA':
        plotter = ICAPlotter.icaplot(
            experiment,
            cmap,
            ax=ax, bokeh=False)
        plotter.icacomp.to_csv(prefix + '.icacomp.txt', sep=SEP)
    elif algorithm == 'KMEANS':
        plotter = KMeansPlotter.kmeansplot(
            experiment,
            cmap,
            ax=ax, bokeh=False)
        plotter.kclusters.to_csv(prefix + '.kcomp.txt', sep=SEP)
    else:
        print("invalid algorithm. Exiting..")
        sys.exit(1)
    leg = plt.legend(loc='best', shadow=False, frameon = 1)

    leg.get_frame().set_edgecolor('b')
    leg.get_frame().set_facecolor('w')
    fig.savefig(output_file)



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