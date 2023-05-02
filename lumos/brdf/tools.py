"""
Tools for working with BRDF models and data
"""

import numpy as np
import scipy.optimize
import lumos.conversions

def fit(
        data_file,
        model_func,
        bounds,
        p0,
        log_space = True,
        clip = 0
        ):
    
    """
    Fits a model to experimental data.

    :param data_file: A path to a .csv file containing BRDF data. DESCRIBE FORMAT.
    :type data_file: str
    :param model_func: A BRDF model to fit
    :type model_func: function
    :param bounds: Bounds passed to :function:`scipy.optimize.curve_fit`
    :type bounds: tuple
    :param p0: Initial guess for parameters passed to :function:`scipy.optimize.curve_fit`
    :type p0: tuple
    :param log_space: Whether or not to fit BRDF in log_space
    :type log_space: bool
    :param clip: Removes BRDF data below this cutoff from fitting
    :type clip: float
    """

    data = np.loadtxt(data_file, skiprows = 1)

    phi_in = np.deg2rad(data[:, 0])
    theta_in = np.deg2rad(data[:, 1])
    phi_out = np.deg2rad(data[:, 2])
    theta_out = np.deg2rad(data[:, 3])
    brdf = data[:, 4]

    mask = brdf > clip

    phi_in = phi_in[mask]
    theta_in = theta_in[mask]
    phi_out = phi_out[mask]
    theta_out = theta_out[mask]
    brdf = brdf[mask]

    indexes = np.arange(brdf.size)

    ix, iy, iz = lumos.conversions.spherical_to_unit(phi_in, theta_in)
    ox, oy, oz = lumos.conversions.spherical_to_unit(phi_out, theta_out)

    def fit_function(idx, *params):
        model_brdf = model_func(*params)

        idx = idx.astype(int)

        f = model_brdf(
            (ix[idx], iy[idx], iz[idx]),
            (0, 0, 1),
            (ox[idx], oy[idx], oz[idx])
            )
        
        return np.log10(f) if log_space else f

    popt, _ = scipy.optimize.curve_fit(fit_function, 
                                       indexes, 
                                       np.log10(brdf) if log_space else brdf, 
                                       bounds = bounds,
                                       p0 = p0)

    return popt

def pack_binomial_parameters(self, n, m, l1, l2, *params):
    """
    Convert list into B & C matrices, which can then be passed to the Binomial Model

    :param n: n
    :type n: int
    :param m: m
    :type m: int
    :param l1: l1
    :type l1: int
    :param l2: l2, ensure l2 > l1
    :type l2: int
    :param params: list of values
    :type params: list[float]
    :return: B, C, d
    :rtype: :class:`np.ndarray`, :class:`np.ndarray`, float
    """
    params = np.array(params)
    B = np.reshape( params[:n * m], (n, m) )
    C = np.reshape( params[n * m:], (n, l2 - l1))
    return B, C