'''
Program 
'''

import math
import scipy.constants
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from scipy import stats
from scipy import optimize
import mesa_reader as mr
import cmath

h = mr.MesaData('LOGS1/history.data')

m1 = 3.5 * 1.9 * 10 ** 33 
m2 = 1.3 * 1.9 * 10 ** 33 
q = m2 / m1

a = 6.68 * 6.96 * 10 ** 10#6.68

G = 6.655 * 10 ** (-8)  # put my units!
omega2 = G * (m1 + m2) / a ** 3
def potential(x, y, z):
	omega2 = G * (m1 + m2) / a ** 3
	r12 = x ** 2 + y ** 2 + z ** 2
	r22 = (x - a) ** 2 + y ** 2 + z ** 2
	xcm = a * q / (1 + q)
	phi = -(G * m1 / np.sqrt(r12)) - (
		G * m2 / np.sqrt(r22)) - omega2 * ((xcm - x) ** 2 + y ** 2
		) / 2
	return phi

x = np.linspace(-8, 12, 1e3) * 6.96 * 10 ** 10#-78
y = np.linspace(-8, 10, 1e3) * 6.96 * 10 ** 10#-77

x_grid, y_grid = np.meshgrid(x, y)

phi = potential(x_grid, y_grid, 1.1e11)


print('phi.max()', phi.max())
print('phi.min()', phi.min())

vmin = -10e25
vmax = -3e15

levels = np.linspace(vmin, vmax, 10e3)

extent = [x.min(), x.max(), y.min(), y.max()]

fig, ax = plt.subplots()

cs = ax.contour(x_grid , y_grid, phi, 10, colors='k')#, level=levels)
ax.imshow(phi, extent=extent, origin='lower', vmin=vmin, vmax=vmax)


lines = []
for line in cs.collections[0].get_paths():
    lines.append(line.vertices)

plt.clabel(cs, inline = False, fontsize=10)
ax.plot(lines[0][:, 0], lines[0][:, 1])


#plt.plot(y, potential_1(y/a, 0, 0) * (-omega2)/2)


plt.show()


