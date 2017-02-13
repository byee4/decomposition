import matplotlib

matplotlib.use('Agg')
import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

import PCAPlotter
from ClusterExperiment import ClusterExperiment

DEBUG = 0
TESTRUN = 0
PROFILE = 0

__all__ = []
__version__ = 0.1
__date__ = '2015-12-19'
__updated__ = '2015-12-19'

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

    parser.add_argument("-s", "--subset", dest="subset", default=None)

    parser.add_argument("-g", "--gene", dest="gene_id", default=None, type=str)
    parser.add_argument("-c", "--conditions", dest="conditions", default=None)
    parser.add_argument("-cc", "--conditions-col", dest="conditions_col", default=None)
    # Process arguments
    args = parser.parse_args()

    counts_file = args.input
    output_file = args.output
    subset_file = args.subset
    conditions_file = args.conditions
    conditions_col = args.conditions_col

    is_featurecounts = args.featureCounts
    is_rpkm = args.rpkm
    is_log2 = args.log2

    sum_cutoff = args.cutoff
    gene_id = args.gene_id

    """ read in counts file """

    counts = pd.read_table(counts_file, index_col=0, comment='#')

    """ trim & transform featureCounts table """
    """
    if is_featureCounts:
        counts = FeatureCountsTable(counts)
        if is_rpkm:
            counts.as_rpkm()
    else:
        counts = CountsTable(counts)
    """
    """ colors """


    # colors = Metadata(conditions_file, conditions_col)





    experiment = ClusterExperiment(
        counts_file=counts_file,
        conditions_file=conditions_file,
        conditions_col=conditions_col,
        gene_id=gene_id,
        is_featurecounts=is_featurecounts,
    )
    """ do pca on select genes only """
    if subset_file and os.path.exists(subset_file):
        experiment.counts.subset(subset_file)
    """ removes rows whos sum (reads or rpkm) < cutoff """
    if sum_cutoff > 0:
        experiment.counts.min_row_sum_cutoff(sum_cutoff)
    """ log2 """
    if is_log2:
        experiment.counts.as_log2(1)
    """ rpkm """
    if is_rpkm:
        experiment.counts.as_rpkm()
    fig, ax = plt.subplots()
    plotter = PCAPlotter.pcaplot(
        experiment.counts.data,
        experiment.cmap,
        colors=experiment.metadata,
        ax=ax, bokeh=False)
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