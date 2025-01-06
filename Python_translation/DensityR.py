import numpy as np
import sys
import pandas as pd
from scipy import integrate, stats

from statsmodels.sandbox.nonparametric import kernels
from statsmodels.tools.decorators import cache_readonly

from statsmodels.tools.validation import array_like, float_like
import statsmodels.nonparametric.bandwidths
from statsmodels.nonparametric.kdetools import silverman_transform, forrt, revrt
#from statsmodels.nonparametric.linbin import fast_linbin
from statsmodels.sandbox.nonparametric import kernels

import pyximport
pyximport.install(setup_args={"script_args" : ["--verbose"]})
try:
    from Federated_Differential_Methylation_Analysis.Python_translation.linbinR import fast_linbin
except ModuleNotFoundError:
    sys.path.append('/home/silke/Documents/FedEWAS')
    from Federated_Differential_Methylation_Analysis.Python_translation.linbinR import fast_linbin
# Kernels Switch for estimators

kernel_switch = dict(
    gau=kernels.Gaussian,
    epa=kernels.Epanechnikov,
    uni=kernels.Uniform,
    tri=kernels.Triangular,
    biw=kernels.Biweight,
    triw=kernels.Triweight,
    cos=kernels.Cosine,
    cos2=kernels.Cosine2,
    tric=kernels.Tricube
)

class KDEUnivariate_rDensity:
    def __init__(self, endog):
        self.endog = array_like(endog, "endog", ndim=1, contiguous=True)
    
    def fit(self, kernel="gau", bw="silverman_r", fft=True, weights=None, gridsize=None, adjust=1, low = None, high = None, cut=3, clip=(-np.inf, np.inf)):
        self.bw_method = bw
        endog = self.endog
        density, grid, bw = kdensityfft_rDensity(endog, kernel=kernel, bw=bw, adjust=adjust, weights=weights, gridsize=gridsize, clip=clip, low = low, high = high, cut = cut)
        self.density = density
        self.support = grid
        self.bw = bw
        self.kernel = kernel_switch[kernel](h=bw)
        return self
    
def kdensityfft_rDensity(x, kernel="gau", bw="silverman_r", weights=None, gridsize=None, adjust=1, clip=(-np.inf, np.inf), low = None, high = None, cut=3, retgrid=True):
    """
        Rosenblatt-Parzen univariate kernel density estimator - with input parameters that match those in the Density() function in the stats package in r.

        Parameters
        ----------
        x : array_like
            The variable for which the density estimate is desired.
        kernel : str
            ONLY GAUSSIAN IS CURRENTLY IMPLEMENTED.
            "gau" for Gaussian.
        bw : str, float, callable
            The bandwidth to use. Choices are:

            - "scott" - 1.059 * A * nobs ** (-1/5.), where A is
            `min(std(x),IQR/1.34)`
            - "silverman" - .9 * A * nobs ** (-1/5.), where A is
            `min(std(x),IQR/1.34)`
            - "normal_reference" - C * A * nobs ** (-1/5.), where C is
            calculated from the kernel. Equivalent (up to 2 dp) to the
            "scott" bandwidth for gaussian kernels. See bandwidths.py
            - If a float is given, its value is used as the bandwidth.
            - If a callable is given, it's return value is used.
            The callable should take exactly two parameters, i.e.,
            fn(x, kern), and return a float, where:

            * x - the clipped input data
            * kern - the kernel instance used

        weights : array or None
            WEIGHTS ARE NOT CURRENTLY IMPLEMENTED.
            Optional  weights. If the x value is clipped, then this weight is
            also dropped.
        gridsize : int
            If gridsize is None, min(len(x), 512) is used. Note that the provided
            number is rounded up to the next highest power of 2.
        adjust : float
            An adjustment factor for the bw. Bandwidth becomes bw * adjust.
            clip : tuple
            Observations in x that are outside of the range given by clip are
            dropped. The number of observations in x is then shortened.
        low: int
            value to be used as the lowest number on the grid
        high: int
            value to be used as the highest number on the grid
        cut : float
            Defines the length of the grid past the lowest and highest values of x
            so that the kernel goes to zero. The end points are
            -/+ cut*bw*{x.min() or x.max()} 
            Will be ignored if low and high are provided
        retgrid : bool
            Whether or not to return the grid over which the density is estimated.

        Returns
        -------
        density : ndarray
            The densities estimated at the grid points.
        grid : ndarray, optional
            The grid points at which the density is estimated.

        Notes
        -----
        Only the default kernel is implemented. Weights are not implemented yet.
        This follows Silverman (1982) with changes suggested by Jones and Lotwick
        (1984). However, the discretization step is replaced by linear binning
        of Fan and Marron (1994). This should be extended to accept the parts
        that are dependent only on the data to speed things up for
        cross-validation.

        References
        ----------
        Fan, J. and J.S. Marron. (1994) `Fast implementations of nonparametric
            curve estimators`. Journal of Computational and Graphical Statistics.
            3.1, 35-56.
        Jones, M.C. and H.W. Lotwick. (1984) `Remark AS R50: A Remark on Algorithm
            AS 176. Kernal Density Estimation Using the Fast Fourier Transform`.
            Journal of the Royal Statistical Society. Series C. 33.1, 120-2.
        Silverman, B.W. (1982) `Algorithm AS 176. Kernel density estimation using
            the Fast Fourier Transform. Journal of the Royal Statistical Society.
            Series C. 31.2, 93-9.
        """   
    x = np.asarray(x)
    # will not work for two columns.
    x = x[np.logical_and(x > clip[0], x < clip[1])]

    # Get kernel object corresponding to selection
    kern = kernel_switch[kernel]()

    if callable(bw):
        bw = float(bw(x, kern))
        # user passed a callable custom bandwidth function
    elif isinstance(bw, str):
        if bw == "silverman_r":
            from scipy.stats import scoreatpercentile
            normalize = 1.34
            IQR = (scoreatpercentile(x, 75) - scoreatpercentile(x, 25)) / normalize
            std_dev = np.std(x, axis=0, ddof=1)
            if IQR > 0:
                signma = np.minimum(std_dev, IQR)
            else:
                signma = std_dev

            bw = 0.9*signma*(len(x)**(-0.2))
        else:
            bw = statsmodels.nonparametric.bandwidths.select_bandwidth(x, bw, kern)
        
    else:
        bw = float_like(bw, "bw")

    bw *= adjust

    nobs = len(x)  # after trim

        # 1 Make grid and discretize the data
    if gridsize is None:
        gridsize = np.max((nobs, 512.0))
    gridsize = 2 ** np.ceil(np.log2(gridsize))  # round to next power of 2
    if low is None:
        a = np.min(x) - cut * bw
    a = low - 4 * bw
    if high is None:
        b = np.max(x) + cut * bw
    b = high + 4 * bw
    grid, delta = np.linspace(a, b, int(gridsize), retstep=True)
    RANGE = b - a

    weights = np.repeat((1/len(x)), len(x))
   
    binned = fast_linbin(x, a, b, gridsize, weights) 

    grid_cords = np.linspace(0, 2*(b-a), num=int(2*gridsize))
    grid_cords_rev = grid_cords[::-1].copy()
    set_with = np.negative(grid_cords_rev[int(gridsize):-1]) # get the negative of the 2 to the nth 
    #element of grid_cords in the reverse order
    grid_cords[(int(gridsize)+1):(2*(int(gridsize+1)))+1] = set_with
    from scipy.stats import norm
    grid_cords_dist = norm.pdf(grid_cords, scale = bw)
    dens_fft = np.fft.fft(binned)
    kords_fft = np.fft.fft(grid_cords_dist)
    kords_fft_con = np.conj(kords_fft)
    kords_fft = np.fft.fft(grid_cords_dist)
    kords_fft_con = np.conj(kords_fft)

    combined = np.fft.ifft(dens_fft*kords_fft_con) * len(binned) 
    
    final_kord = np.maximum(0, np.real(combined)[:int(gridsize)])
    xcords = np.linspace(0,5000, num=int(gridsize))
    out = np.interp(xcords, grid, final_kord)
    
    if retgrid:
        return out, xcords, bw
    else:
        return out, bw

