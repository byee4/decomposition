import numpy as np
import pandas as pd

class ExpressionTable():

    def __init__(self, data_file, lengths_file = 0):
        data = pd.read_table(data_file, index_col=0)
        self.data = data
        self.length = self._set_lengths(lengths_file)
        self._pseudocount = 0
        self._is_log2 = False
        self._is_rpkm = False

    def as_log2(self, pseudocount=0):
        """
        log2 transforms self.data

        Parameters
        ----------
        pseudocount : int
            number of fake counts to add to prevent log2(0) inf problem

        Returns
        -------

        """
        self._is_log2 = True
        self._pseudocount = pseudocount
        self.data = np.log2(self.data + pseudocount)


    def subset(self, subset_file):
        """
        removes any row (gene) that is NOT contained in the subset file

        Parameters
        ----------
        subset_file : basestring
            line-delimited file name of gene ids to subset the raw matrix on.

        Returns
        -------

        """
        attributes = [attr.strip() for attr in open(subset_file, 'r')]
        self.data = self.data.loc[attributes, :].dropna(axis=0)
    def min_row_sum_cutoff(self, min_expr_sum=0):
        """
        removes any row (gene) that does NOT meet the minimum row sum
        of read counts (or RPKM if flagged) requirements.

        Parameters
        ----------
        min_expr_sum : int
            sum across all samples to filter

        Returns
        -------

        """
        self.data = self.data[self.data.sum(axis=1) >= min_expr_sum]
    def as_rpkm(self):
        pass
    def _set_lengths(self, lengths_file):
        pass

class FeatureCountsTable(ExpressionTable):

    def __init__(self, counts_file):
        counts = pd.read_table(counts_file, index_col=0, skiprows=1)
        self.data = counts.ix[:,5:]
        self.length = counts['Length']

        self._pseudocount = 0
        self._is_rpkm = False

    def as_rpkm(self):
        """

        Parameters
        ----------
        df : pandas.DataFrame
            featureCounts table formatted by default as having the following
            columns preceding expression data:

            Geneid	Chr	Start	End	Strand	Length  expression**

        Returns
        -------
        pandas.DataFrame containing rpkm
        """

        counts = self.data
        lengths = self.length
        mapped_reads = counts.sum()
        self.is_rpkm = True
        self.data = (
            counts * pow(10, 9)
        ).div(
            mapped_reads,
            axis=1
        ).div(
            lengths,
            axis=0
        )

    def min_row_sum_cutoff(self, min_expr_sum=0):
        """
        removes any row (gene) that does NOT meet the minimum row sum
        of read counts (or RPKM if flagged) requirements.

        Parameters
        ----------
        min_expr_sum : int
            sum across all samples to filter

        Returns
        -------

        """
        self.data = self.data[self.data.sum(axis=1) >= min_expr_sum]
        self.length = self.length.ix[self.data.index]