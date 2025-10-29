import os

import matplotlib.pyplot as plt
import numpy as np

import pyablo

r = pyablo.XdmfReader()
filename = 'examples/blast.xmf'
if not os.path.exists(filename):
    print(f"File {filename} does not exist.")
    exit(1)

s = r.readSnapshot(filename)
s.print()

y  = np.linspace(0.01, 0.99, 1000)
xz = np.ones_like(y) * 0.5
v = np.stack((xz.T, y.T, xz.T)).T
rho = s.probeDensity(v)
plt.plot(y, rho, '-+')
plt.show()

