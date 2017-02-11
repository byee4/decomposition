import matplotlib

matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from bokeh.models import ColumnDataSource

__all__ = []
__version__ = 0.1
__date__ = '2015-12-19'
__updated__ = '2015-12-19'


class _PCAPlotter():

    def __init__(self, data, colors):
        """

        Parameters
        ----------
        data : pandas.DataFrame
            A table of gene expression in the format (genes, samples)
        colors : pandas.Dataframe
            A table defining the color and condition for each sample name
        prcomp : pandas.DataFrame
            A table of pc components
        source : bokeh.models.ColumnDataSource
            Object which allows set_color() method to interactively update
            colors in bokeh.
        ax : matplotlib.axes._subplots.AxesSubplot

        """
        self.data = data
        self.colors = colors
        self.prcomp = self._fit_transform()
        self.source = self._columnsource()

    def _fit_transform(self):
        """
        Transforms the expression data to principal component space.

        Returns
        -------
        prcomp : pandas.DataFrame
            table containing principle components ordered by variance
        """

        smusher = PCA()
        prcomp = smusher.fit_transform(self.data.T)
        prcomp = pd.DataFrame(prcomp, index=self.data.columns)
        return prcomp

    def _columnsource(self):
        """

        Returns
        -------
        ColumnDataSource : bokeh.models.ColumnDataSource
            Object which allows set_color() method to
            interactively update colors in bokeh.
        """

        return ColumnDataSource(
            data=dict(
                x=self.prcomp[0],
                y=self.prcomp[1],
                idx=self.prcomp.index,
                fill_color=self.colors['color'],
            )
        )

    def _matplotlib(self, cmap='Purples', ax=None):
        """

        Parameters
        ----------
        cmap : basestring
            colormap string
        ax : matplotlib.axes._subplots.AxesSubplot
            subplot axes

        Returns
        -------

        """
        if ax is None:
            ax = plt.gca()
        cmap = plt.get_cmap(cmap)
        for c in set(self.colors['condition']):
            indices = self.prcomp.ix[self.colors[self.colors['condition'] == c].index]
            color = self.colors[self.colors['condition'] == c]['color']
            rgbs = [cmap(n) for n in color]
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
    def set_color(self, colors):
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
        self.source.data['fill_color'] = colors['color']

    def plot(self, bokeh=False, cmap=None, ax=None):
        """

        Parameters
        ----------
        bokeh : Boolean
            True if plotting bokeh figure, else matplotlib axes
        cmap : colormap string
        ax : matplotlib.axes._subplots.AxesSubplot or bokeh.plotting.figure.Figure

        Returns
        -------

        """
        if bokeh:
            self._bokeh(ax)
        else:
            self._matplotlib(cmap, ax)

def pcaplot(data, cmap, colors=None, ax=None, bokeh=False):
    """

    Parameters
    ----------
    data : pandas.DataFrame
        table containing elements as rows, pc as cols
    cmap : basestring
        colormap string
    colors : pandas.DataFrame
        A table defining the color and condition for each sample name
    ax : matplotlib.axes._subplots.AxesSubplot or bokeh.plotting.figure.Figure
    bokeh : Boolean
        True if plotting bokeh figure, else matplotlib axes

    Returns
    -------
    _PCAPlotter object

    """
    plotter = _PCAPlotter(data, colors)
    plotter.plot(bokeh=bokeh, cmap=cmap, ax=ax)
    return plotter