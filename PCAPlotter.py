import matplotlib

matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
# import seaborn as sns
from sklearn.decomposition import PCA

# from bokeh.io import push_notebook, show, output_notebook # for plotting inline
# from bokeh.plotting import figure, output_file
# from bokeh.plotting import Figure, show, output_file # for saving to file
# from ipywidgets import interact
# from bokeh.models import HoverTool
from bokeh.models import ColumnDataSource

__all__ = []
__version__ = 0.1
__date__ = '2015-12-19'
__updated__ = '2015-12-19'


class _PCAPlotter():

    def __init__(self, data, colors, ax):
        """

        Parameters
        ----------
        data : pandas.DataFrame
            A table of gene expression in the format (genes, samples)
        colors : pandas.Dataframe
            A table defining the color and condition for each sample name
        prcomp : pandas.DataFrame
        source : pandas.DataFrame

        """
        self.data = data
        self.colors = colors
        self.prcomp = self._fit_transform()
        self.source = None
        self.ax = ax

    def _fit_transform(self):
        """Transforms the expression data to principal component space"""
        smusher = PCA()
        prcomp = smusher.fit_transform(self.data.T)
        prcomp = pd.DataFrame(prcomp, index=self.data.columns)
        return prcomp

    def _matplotlib(self, cmap='Purples'):
        """
        plots 2d PCA

        :param df: dataframe[gene,sample]
        :param output_file: output .png/.jpg/.svg
        :param colors: list of floats corresponding to how we want to color points
        :param cmap: colormap
        :return: dataframe of prcomps
        """
        if self.ax is None:
            self.ax = plt.gca()
        cmap = plt.get_cmap(cmap)
        for c in set(self.colors['condition']):
            indices = self.prcomp.ix[self.colors[self.colors['condition'] == c].index]
            color = self.colors[self.colors['condition'] == c]['color']
            rgbs = [cmap(n) for n in color]
            self.ax.scatter(indices[0], indices[1], label=c, color=rgbs)
    def _bokeh(self):
        self.source = ColumnDataSource(
            data=dict(
                x=self.prcomp[0],
                y=self.prcomp[1],
                idx=self.prcomp.index,
                fill_color=self.colors['color'],
            )
        )

        self.ax.scatter('x', 'y', radius=0.1,
                   fill_color='fill_color', fill_alpha=0.6,
                   line_color=None, source=self.source)
    def set_color(self, colors):
        self.source.data['fill_color'] = colors['color']

    def plot(self, bokeh=False, cmap=None):
        if bokeh:
            self._bokeh()
        else:
            self._matplotlib(cmap)

def pcaplot(data, cmap, colors=None, ax=None, bokeh=False):
    plotter = _PCAPlotter(data, colors, ax)
    plotter.plot(bokeh=bokeh, cmap=cmap)
    return plotter