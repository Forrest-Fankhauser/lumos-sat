.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Getting Started

   Basic Tutorial <tutorial>
   Advanced Topics <advanced_topics>

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Reference

   API Reference <_autosummary/lumos>

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: For Developers

   Installation <dev_install>


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

For more information, see: 
`Satellite Optical Brightness (Fankhauser et. al 2023) <https://arxiv.org/abs/2305.11123>`_

Why is Satellite Brightness Important?
======================================

After dusk and before dawn, LEO satellites scatter
sunlight onto the the surface of Earth. This scattered light can interfere 
with both casual stargazing and science from large ground-based observatories.
As of 2023, there are more than 6,000 LEO satellites in operation, a 6-fold 
increase over just two years. Many more planned satellite constellations are on the way!
Predicting satellite brightness can help astronomers remove satellite streaks from telescope
images or quantify the impact of satellites on astronomical discovery and can help
satellite operators engineer dimmer satellites!