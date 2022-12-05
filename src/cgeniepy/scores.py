import numpy as np
from scipy.spatial import distance
from .grid import mask_Arctic_Med


class ModelSkill:
    "quantitatively compare similarity metrics between two 2D arrays"

    def __init__(self, model, observation, mask_MedArc=True):

        if mask_MedArc:
            self.model = mask_Arctic_Med(model, policy="na")
            self.data  = mask_Arctic_Med(observation, policy="na")
        else:
            self.model = model
            self.data = observation

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
