import numpy as np
from scipy.spatial import distance
import matplotlib.pyplot as plt

class ArrayComp:
    
    "quantitatively compare similarity metrics between two 2D arrays (model and observation)"

    def __init__(self, model, observation):
        """
        
        :param: model: array-like data
        :param: observation: array-like data

        --------
        Example:
        --------
        
        from cgeniepy.skill import ModelSkill

        ## generate 2d data
        import numpy as np
        x = np.random.rand(36, 36)
        y = np.random.rand(36, 36)

        ## calculate skill score
        ms = ModelSkill(x, y)
        ms.mscore()
        """
        self.model = model
        self.data = observation

        ## check data type, if not float, convert to float
        if self.model.dtype != np.float64:
            self.model = self.model.astype(np.float64)
        if self.data.dtype != np.float64:
            self.data = self.data.astype(np.float64)
        

    def safe_unveil(self, array):
        "get pure array from a numpy masked array object"
        if array.__class__ != np.ma.core.MaskedArray:
            return array
        else:
            return array.filled(np.nan)

    def intersect_index(self, verbose=False):
        """ Return the index where corresponding values are not
        nan in both input arrays. One then can filter the array
        by the output boolean array.
        """

        # If both are not NaN, return it
        array1 = self.safe_unveil(self.model)
        array2 = self.safe_unveil(self.data)

        indx_array = np.logical_and(~np.isnan(array1), ~np.isnan(array2))

        if verbose is True:
            num = indx_array.flatten()[indx_array.flatten() == True].shape[0]
            print("Summary: {} elements simultaneously exist.".format(num))

        return indx_array

    def mscore(self):
        """
        Calculate skill metric M-score. See more in the paper Watterson, I. G. (1996)
        """

        # get common set
        indx = self.intersect_index()
        sub_data1 = self.model[indx]
        sub_data2 = self.data[indx]

        # calculate M-score
        mse = np.square(np.subtract(sub_data1, sub_data2)).mean()
        v1 = sub_data1.var()
        v2 = sub_data2.var()
        g1 = sub_data1.mean()
        g2 = sub_data2.mean()

        mscore = (2 / np.pi) * np.arcsin(1 - mse / (v1 + v2 + np.square(g1 - g2)))

        return mscore

    def pearson_r(self):
        """
        calculate pearson correlation coefficient for 2D array
        """
        # get common set
        indx = self.intersect_index()
        sub_data1 = self.model[indx]
        sub_data2 = self.data[indx]

        # correlation matrix
        corr_mat = np.corrcoef(sub_data1, sub_data2)
        corr = corr_mat[0, 1]

        return corr

    def cos_similarity(self):
        """
        Cosine similarity that ignores the magnitude of data, but focuses on the direction
        """
        # get common set
        indx = self.intersect_index()
        sub_data1 = self.model[indx].ravel()
        sub_data2 = self.data[indx].ravel()

        if sub_data1.mean() != 0 and sub_data2.mean() != 0:
            cos_sim = 1 - distance.cosine(sub_data1.flatten(), sub_data2.flatten())
        else:
            cos_sim = np.nan

        return cos_sim

    def rmse(self):
        """
        Root Mean Sqaure Error (rmse, or rmsd)
        """
        data1 = self.model
        data2 = self.data

        error_2d = data1 - data2
        error_1d = error_2d.ravel()[~np.isnan(error_2d.ravel())]
        rmse = np.sqrt(np.square(error_1d).mean())

        return rmse

    def nrmse(self):
        "Normalised RMSE facilitating the comparison between datasets or models with different scales"

        norm_rmse = self.rmse()/np.nanmean(self.model)
        return norm_rmse

    def crmse(self):
        """
        centred Root Mean Sqaure Error (rmse, or rmsd). See Talor, K. E. (2001) JGR
        """
        # select data
        indx = self.intersect_index()
        sub_data1 = self.model[indx].ravel()
        sub_data2 = self.data[indx].ravel()

        # calculate std
        sigma1 = np.std(sub_data1)
        sigma2 = np.std(sub_data2)

        # pearson correlation
        corr_mat = np.corrcoef(sub_data1, sub_data2)
        corr = corr_mat[0, 1]

        # central rmse
        crmse = sigma1**2 + sigma2**2 - 2 * sigma1 * sigma2 * corr

        return crmse

    def taylor_diagram(self):
        pass
 
    
class DataFrameComp(ArrayComp):

    def __init__(self, df, model_col, observation_col):
        self.model = df[model_col].to_numpy()
        self.data = df[observation_col].to_numpy()
        super().__init__(self.model, self.data)

    def plot(self):
        ## a x-y plot
        fig, ax = plt.subplots()
        ax.scatter(self.model, self.data)
        ax.set_xlabel("Model")
        ax.set_ylabel("Observation")
        ax.set_title("Model vs Observation")

 
# class TaylorDiagram(object):
#     """
#     Taylor diagram.
#     Plot model standard deviation and correlation to reference data in a single-quadrant polar plot,
#     with r=stddev and theta=arccos(correlation).

#     modified from Yannick Copin, https://gist.github.com/ycopin/3342888
#     reference: https://matplotlib.org/stable/gallery/axisartist/demo_floating_axes.html
#     """

#     def __init__(
#             self,
#             fig=None,
#             figscale=1,
#             subplot=111,
#             xmax=None,
#             tmax=np.pi / 2,
#             ylabel="Standard Deviation",
#             rotation=None,
#     ):

#         """
#         Set up Taylor diagram axes, i.e. single quadrant polar
#         plot, using `mpl_toolkits.axisartist.floating_axes`.

#         Parameters:

#         * fig: input Figure or None
#         * subplot: subplot definition
#         * xmax: the length of radius, xmax can be 1.5* reference std
#         """

#         # --------------- tickers --------------------------
#         # Correlation labels (if half round)
#         cor_label = np.array([0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99, 1])

#         # add the negative ticks if more than half round
#         excess_theta = tmax - np.pi / 2
#         if excess_theta > 0:
#             cor_label = np.concatenate((-cor_label[:0:-1], cor_label))

#         # convert to radian
#         rad = np.arccos(cor_label)
#         # tick location
#         gl = gf.FixedLocator(rad)
#         # tick formatting: bind radian and correlation coefficient
#         tf = gf.DictFormatter(dict(zip(rad, map(str, cor_label))))

#         # --------------- coordinate -----------------------
#         # Standard deviation axis extent (in units of reference stddev)
#         # xmin must be 0, which is the centre of round

#         self.xmin = 0
#         self.xmax = xmax
#         self.tmax = tmax

#         # ------- curvilinear coordinate definition -------
#         # use built-in polar transformation (i.e., from theta and r to x and y)
#         tr = PolarAxes.PolarTransform()
#         ghelper = fa.GridHelperCurveLinear(
#             tr,
#             extremes=(0, self.tmax, self.xmin, self.xmax),
#             grid_locator1=gl,
#             tick_formatter1=tf,
#         )

#         # ------- create floating axis -------
#         if fig is None:
#             fig_height = 4.5 * figscale
#             fig_width = fig_height * (1 + np.sin(excess_theta))
#             fig = plt.figure(figsize=(fig_width, fig_height), dpi=100)

#         ax = fa.FloatingSubplot(fig, subplot, grid_helper=ghelper)
#         fig.add_subplot(ax)

#         # Adjust axes
#         # Angle axis
#         ax.axis["top"].label.set_text("Correlation")
#         ax.axis["top"].toggle(ticklabels=True, label=True)
#         # inverse the direction
#         ax.axis["top"].set_axis_direction("bottom")
#         ax.axis["top"].major_ticklabels.set_axis_direction("top")
#         ax.axis["top"].label.set_axis_direction("top")

#         # X axis
#         ax.axis["left"].set_axis_direction("bottom")

#         # Y axis direction & label
#         ax.axis["right"].toggle(all=True)
#         ax.axis["right"].label.set_text(ylabel)
#         ax.axis["right"].set_axis_direction("top")
#         # ticklabel direction
#         ax.axis["right"].major_ticklabels.set_axis_direction("left")

#         ax.axis["bottom"].set_visible(False)

#         # ------- Set instance attribute ----------
#         self.fig = fig
#         # Graphical axes
#         self._ax = ax
#         # grid line
#         self._ax.grid(True, zorder=0, linestyle="--")
#         # aspect ratio
#         self._ax.set_aspect(1)
#         # A parasite axes for further plotting data
#         self.ax = ax.get_aux_axes(tr)
#         # Collect sample points for latter use (e.g. legend)
#         self.samplePoints = []

#     def add_ref(self, refstd, reflabel="Observation", linestyle="-", color="k"):
#         """add a reference point"""
#         self.refstd = refstd
#         # Add reference point
#         # slightly higher than 0 so star can be fully seen
#         l = self.ax.plot(0.01, self.refstd, "k*", ls="", ms=10)
#         # xy for the point, xytext for the text (the coordinates are
#         # defined in xycoords and textcoords, respectively)
#         self.ax.annotate(
#             reflabel,
#             xy=(0.01, self.refstd),
#             xycoords="data",
#             xytext=(-25, -30),
#             textcoords="offset points",
#         )
#         # add stddev contour
#         t = np.linspace(0, self.tmax)
#         r = np.zeros_like(t) + self.refstd
#         self.ax.plot(t, r, linestyle=linestyle, color=color)
#         self.samplePoints.append(l)

#     def add_scatter(self, stddev, corrcoef, *args, **kwargs):
#         """
#         Add sample (*stddev*, *corrcoeff*) to the Taylor
#         diagram. *args* and *kwargs* are directly propagated to the
#         `Figure.plot` command.
#         """

#         l = self.ax.scatter(
#             np.arccos(corrcoef), stddev, *args, **kwargs
#         )  # (theta, radius)
#         self.samplePoints.append(l)

#         return l

#     def add_contours(self, levels=5, **kwargs):
#         """
#         Add constant centered RMS difference contours, defined by *levels*.
#         """

#         rs, ts = np.meshgrid(
#             np.linspace(self.xmin, self.xmax), np.linspace(0, self.tmax)
#         )
#         # Compute centered RMS difference
#         crmse = np.sqrt(self.refstd**2 + rs**2 - 2 * self.refstd * rs * np.cos(ts))
#         contours = self.ax.contour(ts, rs, crmse, levels, linestyles="--", **kwargs)
#         self.ax.clabel(contours, contours.levels[::1], inline=False)

#         return contours

#     def add_legend(self, *args, **kwargs):
#         return self.ax.legend(*args, **kwargs)

#     def add_annotation(self, *args, **kwargs):
#         return self.ax.annotation(*args, **kwargs)

#     def savefig(self, *args, **kwargs):
#         self.fig.savefig(*args, **kwargs)
    
