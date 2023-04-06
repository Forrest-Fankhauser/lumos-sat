""" Simple calculations that are used throughout Lumos """

import numpy as np

def Rx(theta : float, x : float, y : float, z : float) -> tuple[float, float, float]:
    """
    Rotates points (x, y, z) about x-axis by angle theta

    Parameters:
        theta (float) : Angle of rotation (radians)
        x, y, z (float) : Coordinates of starting point
    
    Returns:
        xp, yp, zp (float) : Coordinates of rotated point
    """
    xp = x
    yp = np.cos(theta) * y - np.sin(theta) * z
    zp = np.sin(theta) * y + np.cos(theta) * z
    return xp, yp, zp

def Ry(theta : float, x : float, y : float, z : float) -> tuple[float, float, float]:
    """
    Rotates points (x, y, z) about y-axis by angle theta

    Parameters:
        theta (float) : Angle of rotation (radians)
        x, y, z (float) : Coordinates of starting point
    
    Returns:
        xp, yp, zp (float) : Coordinates of rotated point
    """
    xp = np.cos(theta) * x + np.sin(theta) * z
    yp = y
    zp = -np.sin(theta) * x + np.cos(theta) * z
    return xp, yp, zp

def Rz(theta : float, x : float, y : float, z : float) -> tuple[float, float, float]:
    """
    Rotates points (x, y, z) about z-axis by angle theta

    Parameters:
        theta (float) : Angle of rotation (radians)
        x, y, z (float) : Coordinates of starting point
    
    Returns:
        xp, yp, zp (float) : Coordinates of rotated point
    """
    xp = np.cos(theta) * x - np.sin(theta) * y
    yp = np.sin(theta) * x + np.cos(theta) * y
    zp = z
    return xp, yp, zp

def det_2(a11, a12,
          a21, a22):
    """
    det(A)

    A = |a11 a12|
        |a21 a22|
    """
    return a11 * a22 - a21 * a12

def det_3(a11, a12, a13,
          a21, a22, a23,
          a31, a32, a33):
    """
    det(A)

    A = |a11 a12 a13|
        |a21 a22 a23|
        |a31 a32 a33|
    """
    return a11 * det_2(a22, a23, a32, a33) - a12 * det_2(a21, a23, a31, a33) + a13 * det_2(a21, a22, a31, a32)

def inv_3(a11, a12, a13,
          a21, a22, a23,
          a31, a32, a33):
    """
    A^-1

    A = |a11 a12 a13|
        |a21 a22 a23|
        |a31 a32 a33|
    """
    c = 1 / det_3(a11, a12, a13,
                  a21, a22, a23,
                  a31, a32, a33)

    I11, I12, I13 = c * det_2(a22, a23, a32, a33), c * det_2(a13, a12, a33, a32), c * det_2(a12, a13, a22, a23)
    I21, I22, I23 = c * det_2(a23, a21, a33, a31), c * det_2(a11, a13, a31, a33), c * det_2(a13, a11, a23, a21)
    I31, I32, I33 = c * det_2(a21, a22, a31, a32), c * det_2(a12, a11, a32, a31), c * det_2(a11, a12, a21, a22)

    return I11, I12, I13, \
           I21, I22, I23, \
           I31, I32, I33