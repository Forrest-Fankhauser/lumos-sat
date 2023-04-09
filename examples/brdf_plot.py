import lumos.brdf.library
import lumos.plot.brdf
import matplotlib.pyplot as plt

brdf = lumos.brdf.library.PHONG(0.2, 0.2, 20)

fig = plt.figure()
ax = plt.subplot(polar = True)

lumos.plot.brdf.plot2D(ax, brdf, 50)

plt.show()