import matplotlib

matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from bokeh.models import ColumnDataSource

import color_helpers as ch

__all__ = []
__version__ = 0.1
__date__ = '2017-2-13'
__updated__ = '2017-2-13'


class _TSNEPlotter():

    def __init__(self, expt, cmap = 'Purples',
                 method = 'exact', random_state = 1):
        """

        Parameters
        ----------
        expt : ClusterExperiment
        cmap : basestring
        method : basestring
            gradient calculation algorithm
            @see TSNE(method)
        random_state : int
            tsne random state for tsne seed generator
            @see TSNE(random_state)
        """
        self.method = method
        self.random_state = random_state
        self.cmap = plt.get_cmap(cmap)
        self.expt = expt
        self.tcomp = self._fit_transform()
        self.source = self._columnsource()



    def _fit_transform(self):
        """
        Transforms the expression data to principal component space.

        Returns
        -------
        prcomp : pandas.DataFrame
            table containing principle components ordered by variance
        """

        manifolder = TSNE(method=self.method, random_state = self.random_state)
        tcomp = manifolder.fit_transform(self.expt.counts.data.T)
        tcomp = pd.DataFrame(tcomp, index=self.expt.counts.data.columns)
        return tcomp

    def _columnsource(self):
        """

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
                x=self.tcomp[0],
                y=self.tcomp[1],
                idx=self.tcomp.index,
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

        for c in set(self.expt.metadata['condition']):
            indices = self.tcomp.ix[self.expt.metadata[self.expt.metadata['condition'] == c].index]

            color = self.expt.metadata[self.expt.metadata['condition'] == c]['color']
            rgbs = [self.cmap(n / self.expt.metadata['color'].max()) for n in color]

            ax.scatter(indices[0], indices[1], label=c, color=rgbs)

    def _bokeh(self, ax):
        """

        Parameters
        ----------
        ax : bokeh.plotting.figure.Figure

        Returns
        -------

        """
        ax.scatter('x', 'y', radius=0.1,
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

def tsneplot(expt, cmap, ax=None, bokeh=False):
    """

    Parameters
    ----------
    expt : ClusterExperiment
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
    plotter = _TSNEPlotter(expt, cmap)
    plotter.plot(bokeh=bokeh, ax=ax)
    return plotter