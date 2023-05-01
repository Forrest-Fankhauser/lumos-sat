""" A catalog of useful and common BRDF models """
import numpy as np

def LAMBERTIAN(albedo):
    """
    Lambertian scattering model.

    :math:`BRDF = \\frac{\\rho}{\\pi}`

    :param albedo: Albedo, should be between 0.0 and 1.0
    :type albedo: float
    :return: BRDF function, where f = BRDF(incident vector, normal_vector, outgoing_vector)
    :rtype: function
    """
    def BRDF(incident_vector, normal_vector, outgoing_vector):
        ix, iy, iz = incident_vector
        nx, ny, nz = normal_vector
        ox, oy, oz = outgoing_vector

        brdf = np.where(ix*nx + iy*ny + iz*nz < 0, 0, albedo / np.pi)
        brdf = np.where(ox*nx + oy*ny + oz*nz < 0, 0, albedo / np.pi)

        return brdf
    
    return BRDF

def ABG(A, B, g):
    """
    ABG scattering model. 
    Described here: `Using the ABG model <https://support.zemax.com/hc/en-us/articles/4408106806163-Fitting-ABg-scattering-coefficients-from-raw-measured-data>`_

    :math:`BRDF = \\frac{A}{B + x^g}`

    :param A: Amplitude parameter, A/B is the specular peak
    :type A: float
    :param B: Knee parameter, B determines where the BRDF transitions to exponential decay
    :type B: float
    :param g: Slope parameter, g determines the slope of the BRDF in log-log space
    :type g: float
    :return: BRDF function, where f = BRDF(incident vector, normal_vector, outgoing_vector)
    :rtype: function
    """
    def BRDF(incident_vector, normal_vector, outgoing_vector):
        ix, iy, iz = incident_vector
        nx, ny, nz = normal_vector
        ox, oy, oz = outgoing_vector
        
        dot = ix*nx + iy*ny + iz*nz
        rx = 2 * dot * nx - ix
        ry = 2 * dot * ny - iy
        rz = 2 * dot * nz - iz

        dot = ox*nx + oy*ny + oz*nz
        beta_x = ox - dot * nx
        beta_y = oy - dot * ny
        beta_z = oz - dot * nz

        dot = rx*nx + ry*ny + rz*nz
        beta0_x = rx - dot * nx
        beta0_y = ry - dot * ny
        beta0_z = rz - dot * nz

        x = np.sqrt( (beta_x - beta0_x)**2 + (beta_y - beta0_y)**2 + (beta_z - beta0_z)**2 )

        brdf = A / (B + x**g)

        brdf = np.where(ix*nx + iy*ny + iz*nz < 0, 0, brdf)
        brdf = np.where(ox*nx + oy*ny + oz*nz < 0, 0, brdf)
        
        return brdf
    
    return BRDF

def GAUSSIAN(A, sigma):
    """
    Gaussian scattering model. 

    :math:`BRDF = A exp(\\frac{\\hat{w}_r\\cdot\\hat{w}_o - 1}{\\sigma})`

    :param A: Height of specular peak of BRDF
    :type A: float
    :param sigma: Width of specular peak
    :type sigma: float
    :return: BRDF function, where f = BRDF(incident vector, normal_vector, outgoing_vector)
    :rtype: function
    """

    def BRDF(incident_vector, normal_vector, outgoing_vector):
        ix, iy, iz = incident_vector
        nx, ny, nz = normal_vector
        ox, oy, oz = outgoing_vector
        
        dot = ix*nx + iy*ny + iz*nz
        rx = 2 * dot * nx - ix
        ry = 2 * dot * ny - iy
        rz = 2 * dot * nz - iz

        dot = rx * ox + ry * oy + rz * oz
        brdf = A * np.exp( (dot - 1) / sigma )

        brdf = np.where(ix*nx + iy*ny + iz*nz < 0, 0, brdf)
        brdf = np.where(ox*nx + oy*ny + oz*nz < 0, 0, brdf)
        
        return brdf
    
    return BRDF

def PHONG(Kd, Ks, n):
    """
    Phong scattering model. 

    :math:`BRDF = \\frac{K_d}{\\pi} + K_s\\frac{n + 2}{2 \\pi} (\\hat{w}_r \\cdot \\hat{w}_o)^n`

    :param Kd: Diffuse coefficient. Ensure :math:`Kd + Ks \\leq 1`
    :type Kd: float
    :param Ks: Specular coefficient. Ensure :math:`Kd + Ks \\leq 1`
    :type Ks: float
    :param n: Controls width of specular peak.
    :type n: float
    :return: BRDF function, where f = BRDF(incident vector, normal_vector, outgoing_vector)
    :rtype: function
    """

    def BRDF(incident_vector, normal_vector, outgoing_vector):
        ix, iy, iz = incident_vector
        nx, ny, nz = normal_vector
        ox, oy, oz = outgoing_vector
        
        dot = ix*nx + iy*ny + iz*nz
        rx = 2 * dot * nx - ix
        ry = 2 * dot * ny - iy
        rz = 2 * dot * nz - iz

        dot = rx*ox + ry*oy + rz*oz
        dot = np.clip(dot, 0, 1)

        brdf = Kd / np.pi + Ks * (n + 2) / (2 * np.pi) * dot**n

        brdf = np.where(ix*nx + iy*ny + iz*nz < 0, 0, brdf)
        brdf = np.where(ox*nx + oy*ny + oz*nz < 0, 0, brdf)
        
        return brdf
    
    return BRDF

def BINOMIAL(B, C, d, l1):
    """
    Binomial scattering model proposed by Greynolds.

    Described here: `Physically Realistic BRDF Models <https://www.researchgate.net/publication/300345613_General_physically-realistic_BRDF_models_for_computing_stray_light_from_arbitrary_isotropic_surfaces>`_

    :param B: b_i model coefficients
    :type B: :class:`np.ndarray`
    :param C: c_i model coefficients
    :type C: :class:`np.ndarray`
    :param d: specularity constant
    :type d: float
    :param l1: Minimum Gaussian index
    :type l1: int
    :return: BRDF function, where f = BRDF(incident vector, normal_vector, outgoing_vector)
    :rtype: function
    """

    def BRDF(incident_vector, normal_vector, outgoing_vector):
        ix, iy, iz = incident_vector
        nx, ny, nz = normal_vector
        ox, oy, oz = outgoing_vector

        dot = ix*nx + iy*ny + iz*nz
        rx = 2 * dot * nx - ix
        ry = 2 * dot * ny - iy
        rz = 2 * dot * nz - iz

        dot = ox*nx + oy*ny + oz*nz
        rho_x = ox - dot * nx
        rho_y = oy - dot * ny
        rho_z = oz - dot * nz

        dot = rx*nx + ry*ny + rz*nz
        rho0_x = rx - dot * nx
        rho0_y = ry - dot * ny
        rho0_z = rz - dot * nz

        D = np.sqrt( (rho_x - rho0_x)**2 + (rho_y - rho0_y)**2 + (rho_z - rho0_z)**2 )
        V = rho_x * rho0_x + rho_y * rho0_y + rho_z * rho0_z

        log_brdf = np.zeros_like(D)
        
        n, m = B.shape
        _, j = C.shape

        for k in range(n):
            term_1 = np.zeros_like(D)

            for i in range(m):
                term_1 = term_1 + B[k, i] * D ** i
            
            term_2 = np.zeros_like(D)
            for i in range(j):
                term_2 = term_2 + C[k, i] * np.log10(1 + d**(i + l1) * D**2)
            
            log_brdf = log_brdf + (term_1 + 0.5 * term_2) * V**k
        
        brdf = 10**log_brdf

        brdf = np.where(ix*nx + iy*ny + iz*nz < 0, 0, brdf)
        brdf = np.where(ox*nx + oy*ny + oz*nz < 0, 0, brdf)
        
        return brdf
    
    return BRDF