import matplotlib

matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from bokeh.models import ColumnDataSource
import seaborn as sns
import color_helpers as ch
import numpy as np
import pandas as pd

__all__ = []
__version__ = 0.1
__date__ = '2015-12-19'
__updated__ = '2015-12-19'


class _PCAPlotter():

    def __init__(self, expt, cmap = 'Purples'):
        """

        Parameters
        ----------
        data : pandas.DataFrame
            A table of gene expression in the format (genes, samples)
        cmap : matplotlib.colors.Colormap
            colormap instance corresponding to the name

        Attributes
        ----------
        self.cmap : matplotlib.colors.Colormap
            colormap instance corresponding to the name
        self.expt : Experiment.Experiment
            object containing info about the decomposition expt. This
            includes: self.expt.counts (counts table or matrix to be fitted)
            and self.expt.metadata (labels specifying the grouping of each
            column).
        self.pca : sklearn.decomposition.PCA
        self.source : bokeh.models.ColumnDataSource

        """
        self.cmap = plt.get_cmap(cmap)
        self.expt = expt
        self.pca, self.prcomp = self._fit_transform()
        self.source = self._columnsource()

    def get_pc_components(self):
        """
        Returns a DataFrame of how much each feature contributes to the PC.

        Returns
        -------
        pc_components : pandas.DataFrame

        """
        pc_cols = range(len(self.pca.components_))
        pc_components = pd.DataFrame(
            index=self.expt.counts.data.index,
            columns=pc_cols
        )

        for n in range(0, len(self.pca.components_)):
            for i, j in zip(
                    self.expt.counts.data.index,
                    np.abs(self.pca.components_[n])
            ):
                pc_components.ix[i, n] = j
        return pc_components

    def _fit_transform(self):
        """
        Transforms the expression data to principal component space.

        Returns
        -------
        prcomp : pandas.DataFrame
            table containing principle components ordered by variance
        """

        smusher = PCA()
        prcomp = smusher.fit_transform(self.expt.counts.data.T)
        prcomp = pd.DataFrame(prcomp, index=self.expt.counts.data.columns)

        return smusher, prcomp

    def _columnsource(self):
        """
        Creates and returns the ColumnDataSource object needed by Bokeh plots.

        Returns
        -------
        ColumnDataSource : bokeh.models.ColumnDataSource
            Object which allows set_color() method to
            interactively update colors in bokeh.
        """
        self.expt.metadata['hex'] = ch.expr_series_to_hex(
            self.expt.metadata['color'],
            self.cmap,
            is_norm=True
        )

        return ColumnDataSource(
            data=dict(
                x=self.prcomp[0],
                y=self.prcomp[1],
                idx=self.prcomp.index,
                fill_color=self.expt.metadata['hex'],
            )
        )

    def _matplotlib(self, ax=None):
        """

        Parameters
        ----------
        ax : matplotlib.axes._subplots.AxesSubplot
            subplot axes

        Returns
        -------

        """
        if ax is None:
            ax = plt.gca()

        colors = sns.color_palette(
            "hls", len(set(self.expt.metadata['color']))
        )
        i = 0
        for c in set(self.expt.metadata['condition']):

            indices = self.prcomp.ix[self.expt.metadata[self.expt.metadata['condition'] == c].index]
            ax.scatter(indices[0], indices[1], label=c, color=colors[i])
            i += 1

    def _bokeh(self, ax):
        """

        Parameters
        ----------
        ax : bokeh.plotting.figure.Figure

        Returns
        -------

        """
        print(ax)
        ax.scatter('x', 'y', radius=0.2,
                   fill_color='fill_color', fill_alpha=0.6,
                   line_color=None, source=self.source)

    def set_color(self, gene_id):
        """
        Updates self.ColumnDataSource 'fill_color' column to interactively
        change point colors.

        Parameters
        ----------
        colors : pandas.DataFrame
            Table describing samples as row indices, 'condition' and
            corresponding 'color' as columns

        Returns
        -------

        """
        print('setting color')
        self.expt.recolor(gene_id)
        self.expt.metadata['hex'] = ch.expr_series_to_hex(
            self.expt.metadata['color'],
            self.cmap,
            is_norm=True
        )
        self.source.data['fill_color'] = self.expt.metadata['hex']

    def plot(self, bokeh=False, ax=None):
        """

        Parameters
        ----------
        bokeh : Boolean
            True if plotting bokeh figure, else matplotlib axes
        ax : matplotlib.axes._subplots.AxesSubplot or bokeh.plotting.figure.Figure

        Returns
        -------

        """
        if bokeh:
            self._bokeh(ax)
        else:
            self._matplotlib(ax)

    def update_cmap(self):
        pass


def pcaplot(expt, cmap, ax=None, bokeh=False):
    """

    Parameters
    ----------
    expt : Experiment
        Object defining the expression data and conditions for samples.
    cmap : basestring
        colormap string
    ax : matplotlib.axes._subplots.AxesSubplot or bokeh.plotting.figure.Figure
    bokeh : Boolean
        True if plotting bokeh figure, else matplotlib axes

    Returns
    -------
    _PCAPlotter object

    """
    plotter = _PCAPlotter(expt, cmap)
    plotter.plot(bokeh=bokeh, ax=ax)
    return plotter