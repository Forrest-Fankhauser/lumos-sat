Brightness in the Satellite Frame
---------------------------------

Next, we want to find the brightness in the satellite-centered reference frame.
This frame is useful for satellite operators because it shows how different designs
change brightness very intuitively. It is also helpful to use this frame to make sure your
satellite brightness model is set up correctly.

There are 4 steps to do this:

    1. Set up a mesh of observers on the Earth's surface.
    2. Calculate the intensity seen by each observer
    3. Convert intensity to AB Magnitude.
    4. Plot AB Magnitude.

This is implemented as follows. We calculate brightness for 4 angles
past the terminator, (-10, 0, 10, 20):

.. literalinclude:: ../../examples/satellite_frame.py
    :linenos:

Understanding the Satellite Frame
=================================
In the above images, the sunlight is coming from the right side. The terminator is shown as a 
vertical line. The daytime side of Earth is to the right side and the nighttime side is on the left.
The circle shows the horizon of Earth as seen by the satellite.

Rotating a Surface
==================
To get more intution for the satellite-centered frame, we can see what effect rotating the chassis
has on brightness. First, we angle the chassis at 45 degrees towards the sun:

.. code-block::

    satellite.SURFACES[0].normal = np.array([0, 0.707, -0.707])
    plot_brightness()

Next, we rotate the chassis 45 degrees on the y-axis:

.. code-block::

    satellite.SURFACES[0].normal = np.array([0.707, 0, -0.707])
    plot_brightness()