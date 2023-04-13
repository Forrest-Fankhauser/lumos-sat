import numpy as np

def LAMBERTIAN(TIS : float = 1.0) -> callable:
    """
    Lambertian scattering model

    Parameters:
        TIS (float) : Total integrated scatter, should be between 0.0 and 1.0
    Returns:
        BRDF (callable) : BRDF function which can be used to create a lumos.geometry.Surface object
    """
    def BRDF(incident_vector, normal_vector, outgoing_vector):
        ix, iy, iz = incident_vector
        nx, ny, nz = normal_vector
        ox, oy, oz = outgoing_vector

        brdf = np.where(ix*nx + iy*ny + iz*nz < 0, 0, TIS / np.pi)
        brdf = np.where(ox*nx + oy*ny + oz*nz < 0, 0, TIS / np.pi)

        return brdf
    
    return BRDF

def ABG(A : float = 0.1, B : float = 0.1, g : float = 1.0) -> callable:
    """
    ABg scattering model

    Parameters:
        A (float) : Amplitude parameter, A/B is the specular peak
        B (float) : Knee parameter, B determines where the BRDF transitions to exponential decay
        g (float) : Slope parameter, g determines the slope of the BRDF in log-log space
    Returns:
        BRDF (callable) : BRDF function which can be used to create a lumos.geometry.Surface object
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

def GAUSSIAN(A = 1.0, sigma = 1.0) -> callable:
    """
    Gaussian scattering model

    Parameters:
        A (float) : Specular peak of BRDF
        sigma (float) : Width of specular peak
    Returns:
        BRDF (callable) : BRDF function which can be used to create a lumos.geometry.Surface object
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

def PHONG(Kd : float = 0.5, Ks : float = 0.5, n : float = 1) -> callable:
    """
    Phong BRDF model

    Parameters:
        Kd (float) : Diffuse coefficient. Ensure Kd + Ks <= 1
        Ks (float) : Specular coefficient. Ensure Kd + Ks <= 1
        n (int) : Specularity. Controls width of specular peak.
    Returns:
        BRDF (callable) : BRDF function which can be used to create a lumos.geometry.Surface object
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

def BINOMIAL(B : np.ndarray, C : np.ndarray, d : float, l1 : int) -> callable:
    """
    Binomial scattering model proposed by Greynolds.

    Parameters:
        B (np.ndarray) : b_i model coefficients
        C (np.ndarray) : c_i model coefficients
        d (np.ndarray) : specularity constant
        l1 (int) : Minimum Gaussian index
    Returns:
        BRDF (callable) : BRDF function which can be used to create a lumos.geometry.Surface object
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

class BinomialHelper:
    def __init__(self, n, m, l1, l2):
        self.n = n
        self.m = m
        self.l1 = l1
        self.l2 = l2
        self.N_params = n * m + (l2 - l1) * n + 1

    def pack_params(self, *params):
        params = np.array(params)
        B = np.reshape( params[:self.n * self.m], (self.n, self.m) )
        C = np.reshape( params[self.n * self.m:-1], (self.n, self.l2 - self.l1))
        d = np.abs(params[-1])
        return B, C, d