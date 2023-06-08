import numpy as np
import lumos.brdf.tools
import matplotlib.pyplot as plt
import lumos.plot
from lumos.brdf.library import ABG

A, B, g = lumos.brdf.tools.fit(
    "examples/aluminum_brdf.csv",
    ABG,
    bounds = (0, 1e4),
    p0 = (1, 1, 1))

print(f"{A = :0.3e}")
print(f"{B = :0.3e}")
print(f"{g = :0.3f}")

fig, ax = plt.subplots()

phi_i, theta_i, phi_o, theta_o, brdf = lumos.brdf.tools.read_brdf("examples/aluminum_brdf.csv")
ax.semilogy(phi_o, brdf, 'k.', alpha = 0.5, zorder = 2)
lumos.plot.BRDF_1D(ax, ABG(A, B, g), incident_angles = np.unique(phi_i))

plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.xlabel("Outgoing Zenith Angle", fontsize = 16)
plt.ylabel(r"BRDF $(sr^{-1})$", fontsize = 16)
plt.title("Aluminum BRDF", fontsize = 20)
plt.legend(fontsize = 14)
plt.show()