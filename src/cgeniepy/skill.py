from typing import Union
import numpy as np
from scipy.spatial import distance
from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt

from matplotlib.projections import PolarAxes
import mpl_toolkits.axisartist.floating_axes as fa
import mpl_toolkits.axisartist.grid_finder as gf


class ArrComparison:
    
    "quantitatively compare similarity metrics between two 2D arrays (model and observation)"

    def __init__(self, model, observation, model_name="Model", obs_name="Observation", label=None):
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

        self.model_name = model_name
        self.obs_name = obs_name
        self.label = label

        ## if is not numeric, convert to numeric
        if not np.issubdtype(self.model.dtype, np.number):
            self.model = self.model.astype(float)
        if not np.issubdtype(self.data.dtype, np.number):
            self.data = self.data.astype(float)
            
        

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

    def spearman_r(self):
        """
        calculate spearman correlation
        """
        indx = self.intersect_index()
        sub_data1 = self.model[indx]
        sub_data2 = self.data[indx]

        corr, _ = spearmanr(sub_data1, sub_data2)

        return corr

    def mae(self):
        """
        Mean Absolute Error (MAE)
        """
        # get common set
        indx = self.intersect_index()
        
        sub_data1 = self.model[indx]
        sub_data2 = self.data[indx]

        mae = np.abs(sub_data1 - sub_data2).mean()

        return mae

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

    def mse(self):
        """
        Mean Squared Error (MSE)
        """
        data1 = self.model
        data2 = self.data

        error_2d = data1 - data2
        error_1d = error_2d.ravel()[~np.isnan(error_2d.ravel())]
        mse = np.square(error_1d).mean()

        return mse

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
    
    def lm_pvalue(self, method='pearson'):
        """
        get p-value from linear regression
        """
        
        indx = self.intersect_index()
        sub_data1 = self.model[indx].ravel()
        sub_data2 = self.data[indx].ravel()        

        if method == 'pearson':
            _, p = pearsonr(sub_data1, sub_data2)
        elif method == 'spearman':
            _, p = spearmanr(sub_data1, sub_data2)
        else:
            raise ValueError("method should be either 'pearson' or 'spearman'")            
        
        
        return p

    def plot(self, ax=None, diagonal=True, log_scale=False,
             title="Model vs Observation",             
             savefig_name=None, metric=['rmse', 'r2', 'm', 'p'], *args, **kwargs):
        ## a x-y plot
        if not ax:
            fig, ax = plt.subplots()
        plt.rcParams['font.family'] = 'sans-serif'

        if log_scale:
            ax.set_xscale('log')
            ax.set_yscale('log')
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.scatter(self.model, self.data, zorder=2,*args, **kwargs)
        ax.minorticks_on()
        ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
        ax.set_xlabel(self.model_name)
        ax.set_ylabel(self.obs_name)
        ## add 1:1 line
        if diagonal:            
            ## add real 1:1 line
            ax.plot(ax.get_xlim(), ax.get_xlim(), ls="-", c=".3")            
            
        
        ## add metrics
        if 'rmse' in metric:
            ax.text(0.05, 0.8, "RMSE: {:.3f}".format(self.rmse()), transform = ax.transAxes)
            
        if 'r2' in metric:                 
             ax.text(0.05, 0.85, "Pearson R: {:.3f}".format(self.pearson_r()),transform = ax.transAxes)
                
        if 'm' in metric:
            ax.text(0.05, 0.75, "M-score: {:.3f}".format(self.mscore()),transform = ax.transAxes)
            
        if 'p' in metric:
            ax.text(0.05, 0.9, "p-value: {:.4f}".format(self.lm_pvalue()),transform = ax.transAxes)
            
            
        
        if self.label:
            title = title + " ({})".format(self.label)
            
            

        ax.set_title(title, loc='left', fontweight='bold', fontsize=12)

        if savefig_name:
            fig.savefig(savefig_name)

    
class DFComparison(ArrComparison):

    """
    Quantitatively compare similarity metrics between two columns in a dataframe-like objects
    """

    def __init__(self, df, model_col, observation_col, *args, **kwargs):
        self.model = df[model_col].to_numpy()
        self.data = df[observation_col].to_numpy()
        super().__init__(self.model, self.data, *args, **kwargs)

 
class TaylorDiagram(object):
    """
    Taylor diagram.

    A visualisation of model-data comparison considering (1) RMSE; (2) Correlation; (3) Standard Deviation.
    
    In practical it is a polar axis with point with `std` as radiance and `corr` as theta.
    The cRMSE then measures the distance from reference point, which is at (ref_std, 0)

    see https://en.wikipedia.org/wiki/Taylor_diagram#/media/File:Taylor_diagram_fig2.png
    and https://bookdown.org/david_carslaw/openair/sections/model-evaluation/taylor-diagram.html (Figure 20.2)

    This class is modified from Yannick Copin, https://gist.github.com/ycopin/3342888
    """

    def __init__(self, ac: Union[ArrComparison,list]):
        if isinstance(ac, list):
            self.mult_comp = True
        else:
            self.mult_comp = False
            
        self.ac = ac
        self.extract_data()

    def extract_data(self):

        if not self.mult_comp:        
            self.corr = self.ac.pearson_r()
            self.model_std = np.nanstd(self.ac.model)
            self.obs_std = np.nanstd(self.ac.data)
            self.label = self.ac.label
        else:
            corr, model_std, obs_std, label = [],[],[],[]
            for ac_i in self.ac:
                corr.append(ac_i.pearson_r())
                model_std.append(np.nanstd(ac_i.model))
                obs_std.append(np.nanstd(ac_i.data))
                label.append(ac_i.label)                

            self.corr = corr
            self.model_std = model_std
            self.obs_std = obs_std
            self.label = label
                

    def setup_ax(self, fig=None, positive_only=True, crmse_contour=False):

        tr = PolarAxes.PolarTransform()

        ## 1.creating extrems for gridhelper
          ## ymax = 1/2 pi if
        if positive_only:
            ymax= 1/2 * np.pi
        else:
            ymax = np.pi

        x0 = 0
        if self.mult_comp:
            self.ref_std = 1
        else:
            self.ref_std = self.obs_std

        x1 = self.ref_std * 1.6
        extremes=(0, ymax, x0, x1)

        # 2. grid locator     
        rlocs = np.array([0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99, 1])
        # a concatenate with ticks' reverse
        rlocs = np.concatenate((-rlocs[:0:-1], rlocs))
        ## arc-cos -> convert to polar angles
        tlocs = np.arccos(rlocs)
        grid_locator1 = gf.FixedLocator(tlocs)    # Positions

        # 3. tick formatters for each axis
        tf1 = gf.DictFormatter(dict(zip(tlocs, map(str, rlocs))))

        ## gridhelper to create the special axis
        gridhelper = fa.GridHelperCurveLinear(
                tr, extremes=extremes,
                grid_locator1=grid_locator1,
                tick_formatter1=tf1)

        ## create figure
        if not fig:
            self.fig = plt.figure(dpi=300)
        else:
            self.fig = fig
            
        ax = self.fig.add_axes(111, axes_class=fa.FloatingAxes, grid_helper=gridhelper)

        # Adjust axes
        plt.rcParams['font.family'] = 'sans-serif'
        ax.axis["top"].set_axis_direction("bottom")   # "Angle axis"
        ax.axis["top"].toggle(ticklabels=True, label=True)
        ax.axis["top"].major_ticklabels.set_axis_direction("top")
        ax.axis["top"].label.set_axis_direction("top")
        ax.axis["top"].label.set_text("Correlation")

        ax.axis["left"].set_axis_direction("bottom")  # "X axis"
        ax.axis["left"].label.set_text("Standard deviation")

        ax.axis["right"].set_axis_direction("top")    # "Y-axis"
        ax.axis["right"].toggle(ticklabels=True)

        ax.axis["bottom"].toggle(ticklabels=False, label=False)

        ## add grid lines
        ax.grid(linestyle='dashed', axis='x')


        self._ax = ax # Graphical axes
        self.ax = ax.get_aux_axes(tr) # Polar coordinates

        if crmse_contour:
            rs, ts = np.meshgrid(np.linspace(x0, x1), np.linspace(0, ymax))
            RMSE=np.sqrt(np.power(self.ref_std, 2) + np.power(rs, 2) - (2.0 * self.ref_std * rs  *np.cos(ts)))
            contours = self.ax.contour(ts, rs, RMSE, 5, linestyles='dashed', colors='#d4af37', linewidths=1.5,
                                       zorder=1)
            ## add label
            self.ax.clabel(contours, inline=True, fontsize=10)

            ## add reference std line
            theta_ref = np.linspace(0, 2*np.pi, 100)
            radiance_ref = theta_ref*0 + self.ref_std
            self.ax.plot(theta_ref, radiance_ref, color='lightblue', linestyle='dashed',
                         linewidth=1.5)

            ## add reference point (because crmse is the distance to reference point)
            self.add_point(1, self.ref_std, marker='*',label='ref', edgecolor='k')            


    def add_point(self, correlation, std, *args, **kwargs):
        ## x->theta, y->radiance        
        self.ax.scatter(np.arccos(correlation), std,*args, **kwargs)        

    def plot(self, cmap=None,add_legend=True, *args, **kwargs):
        if not self.mult_comp:
            ## add model point
            self.add_point(self.corr, self.model_std, label=self.label, *args, **kwargs)
        else:
            if not cmap:
                color_list = plt.cm.get_cmap('tab10', len(self.ac)).colors
            else:
                color_list = cmap.colors            
            
            for i in range(len(self.ac)):
                ## normalised std
                self.add_point(self.corr[i], self.model_std[i]/self.obs_std[i], label=self.label[i],
                              c=color_list[i],  *args, **kwargs)

        if add_legend:
            # outside the box
            self.ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.0))


    def savefig(self, *args, **kwargs):
        self.fig.savefig(*args, **kwargs)

    def add_line(self,theta1, radius1,theta2, radius2, *args, **kwargs):
        """
        Add a line in the polar axis

        If you need to convert the correlation into theta/radian, use np.arccos(correlation)
        """
        line_2d = plt.Line2D([theta1,theta2], [radius1,radius2])
        self.ax.add_line(line_2d, *args, **kwargs)
