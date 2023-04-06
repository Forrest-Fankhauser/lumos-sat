"""
Tools for fitting experimental BRDF data to models
"""

import numpy as np
import scipy.optimize
import lumos.conversions

def fit_model(
        data_file : str, 
        model_func : callable, 
        p0 : tuple[float, ...], 
        bounds : tuple[float, ...],
        log_space : bool = True,
        clip : float = 1e-4
        ) -> tuple[float, ...]:
    
    """
    Fits a model to experimental data.

    Parameters:
        data_file (str) : File containing experimental data.
        model_func (callable) : Given parameters *params, returns BRDF callable
        p0 (tuple) : Initial guess for fitting parameters, passed to scipy.optimize.curve_fit
        bounds (tuple) : Bounds for fitting parameters, passed to scipy.optimize.curve_fit
        clip (float) : Clips BRDF data below this value
    
    Returns:
        popt (tuple) : Optimal parameters returned from scipy.optimize.curve_fit
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