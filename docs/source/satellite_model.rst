Creating a Satellite Model
--------------------------

Once we have found the BRDFs for our satellite surfaces, we can create a brightness model
for our satellite. This file simply lists all the primary surfaces on the satellite.
Each surface is defined by its area, normal vector, and BRDF. Note that normal vectors
must be given in the satellite-centered frame, where :math:`\hat{z}` points towards geodetic zenith
and :math:`\hat{y}` is perpendicular to :math:`\hat{z}` and towards the sun.

.. literalinclude:: ../../examples/satellite.py
    :linenos:

Voila! We've created a simple brightness model for SATELLITE!

Actuating Surfaces
==================

It is worth mentioning that normal vectors can also be defined as functions of the
angle past terminator. For example, to make our solar array track directly towards
the sun, we would define the solar array surface as follows:

.. code-block::

    def solar_array_normal(angle_past_terminator):
        return np.cos(angle_past_terminator) * y_hat - np.sin(angle_past_terminator) * z_hat

    solar_array = Surface(
        area = 1.0, # meters^2
        normal = solar_array_normal, # Directly towards the sun
        brdf = lumos.brdf.library.PHONG(1, 1, 1) # BRDF of solar array
        )