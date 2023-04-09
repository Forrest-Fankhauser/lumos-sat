""" Simple mathematical calculations that are used throughout Lumos """

import numpy as np

def Rx(theta, x, y, z):
    """
    Rotates points :math:`(x, y, z)` about x-axis by angle :math:`\\theta`.
    
    :param theta:
    :type theta: :class:`np.ndarray` or float
    :param x:
    :type x: :class:`np.ndarray` or float
    :param y:
    :type y: :class:`np.ndarray` or float
    :param z:
    :type z: :class:`np.ndarray` or float
    :return: Rotated points :math:`(x^{\\prime}, y^{\\prime}, z^{\\prime})`
    """
    xp = x
    yp = np.cos(theta) * y - np.sin(theta) * z
    zp = np.sin(theta) * y + np.cos(theta) * z
    return xp, yp, zp

def Ry(theta, x, y, z):
    """
    Rotates points :math:`(x, y, z)` about y-axis by angle :math:`\\theta`.
    
    :param theta:
    :type theta: :class:`np.ndarray` or float
    :param x:
    :type x: :class:`np.ndarray` or float
    :param y:
    :type y: :class:`np.ndarray` or float
    :param z:
    :type z: :class:`np.ndarray` or float
    :return: Rotated points :math:`(x^{\\prime}, y^{\\prime}, z^{\\prime})`
    """
    xp = np.cos(theta) * x + np.sin(theta) * z
    yp = y
    zp = -np.sin(theta) * x + np.cos(theta) * z
    return xp, yp, zp

def Rz(theta, x, y, z):
    """
    Rotates points :math:`(x, y, z)` about z-axis by angle :math:`\\theta`.
    
    :param theta:
    :type theta: :class:`np.ndarray` or float
    :param x:
    :type x: :class:`np.ndarray` or float
    :param y:
    :type y: :class:`np.ndarray` or float
    :param z:
    :type z: :class:`np.ndarray` or float
    :return: Rotated points :math:`(x^{\\prime}, y^{\\prime}, z^{\\prime})`
    """
    xp = np.cos(theta) * x - np.sin(theta) * y
    yp = np.sin(theta) * x + np.cos(theta) * y
    zp = z
    return xp, yp, zp

def det_2(a11, a12,
          a21, a22):
    """
    Calculates determinate of 2x2 matrix. Matrix components can be of type
    :class:`np.ndarray` or float.
    
    Returns :math:`\\det \\begin{pmatrix} a_{11} & a_{12} \\\\ a_{21} & a_{22}\\end{pmatrix}`
    """
    return a11 * a22 - a21 * a12

def det_3(a11, a12, a13,
          a21, a22, a23,
          a31, a32, a33):
    """
    Calculates determinate of 3x3 matrix. Matrix components can be of type
    :class:`np.ndarray` or float.

    Returns :math:`\\det \\begin{pmatrix} a_{11} & a_{12} & a_{13} \\\\ a_{21} & a_{22} & a_{23} \\\\ a_{31} & a_{32} & a_{33} \\end{pmatrix}`
    """
    return a11 * det_2(a22, a23, a32, a33) - a12 * det_2(a21, a23, a31, a33) + a13 * det_2(a21, a22, a31, a32)

def inv_3(a11, a12, a13,
          a21, a22, a23,
          a31, a32, a33):
    """
    Calculates inverse of 3x3 matrix. Matrix components can be of type
    :class:`np.ndarray` or float.

    Returns :math:`\\begin{pmatrix} a_{11} & a_{12} & a_{13} \\\\ a_{21} & a_{22} & a_{23} \\\\ a_{31} & a_{32} & a_{33} \\end{pmatrix}^{-1}`
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