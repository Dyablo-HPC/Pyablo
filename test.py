import pydy
import numpy as np
import matplotlib.pyplot as plt
import os

r = pydy.XdmfReader()
filename = '../examples/blast_3d.xmf'
if not os.path.exists(filename):
    print('ARG')
    exit(1)

s = r.readSnapshot(filename)
s.print()

y  = np.linspace(0.01, 0.99, 1000)
xz = np.ones_like(y) * 0.5
v = np.stack((xz.T, y.T, xz.T)).T
rho = s.probeDensities(v)
plt.plot(y, rho, '-+')
plt.show()

