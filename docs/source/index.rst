.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Getting Started

   Installation <user_install>
   Bidirectional Reflectance Distribution Functions (BRDFs) <brdfs>
   Creating a Satellite Model <satellite_model>
   Brightness in the Satellite Frame <satellite_frame>
   Brightness in the Observer Frame <observer_frame>

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Development

   Installation <dev_install>

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Reference

   API Reference <_autosummary/lumos>

Welcome to the Lumos documentation
==================================

**Lumos** is a Python library for calculating the optical brightness of artificial satellites.

More specifically, Lumos allows the user to:

- Fit analytic BRDFs (bidirectional reflectance distribution functions) to measured data.
- Create brightness models of satellites.
- Calculate satellite optical brightness for observers on the ground.

Possible applications include:

- Allowing optical brightness to be used as a design constraint by satellite operators.
- Lets astronomers plan observations based on when and where in the sky satellites are brightest.