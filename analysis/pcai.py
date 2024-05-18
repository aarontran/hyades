#!/usr/bin/env python
"""
Inspect PCAI test problem
"""

from __future__ import division, print_function

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt

import h5py

fname = "flds.all.hdf5"
dt = 0.01
idumpf = 50  # beware hardcoded from input parameters
Lz = 10.5

with h5py.File(fname, "r") as f:
    bx = f['bx'][()]
    by = f['by'][()]
    bz = f['bz'][()]
    ex = f['ex'][()]
    ey = f['ey'][()]
    ez = f['ez'][()]
    jfx = f['jfx'][()]
    jfy = f['jfy'][()]
    jfz = f['jfz'][()]
    rhof = f['rhof'][()]

print("field array shape", bx.shape)
assert bx.shape[1] == bx.shape[2] < bx.shape[3]

bx   = np.mean(bx  , axis=(1,2))
by   = np.mean(by  , axis=(1,2))
bz   = np.mean(bz  , axis=(1,2))
ex   = np.mean(ex  , axis=(1,2))
ey   = np.mean(ey  , axis=(1,2))
ez   = np.mean(ez  , axis=(1,2))
jfx  = np.mean(jfx , axis=(1,2))
jfy  = np.mean(jfy , axis=(1,2))
jfz  = np.mean(jfz , axis=(1,2))
rhof = np.mean(rhof, axis=(1,2))

zvec = np.arange(0, bx.shape[-1]) * (Lz/(bx.shape[-1]-1))
time = np.arange(0, bx.shape[0]) * dt * idumpf

print("time span", time[0], time[-1], "in units of Omega_ci^-1")
print("z span", zvec[0], zvec[-1], "in units of d_i")

plt.figure(figsize=(4,2))
plt.imshow(
    bx.T,
    origin='lower',
    extent=(
        time[0] - 0.5*np.diff(time)[0], time[-1] + 0.5*np.diff(time)[-1],
        zvec[0] - 0.5*np.diff(zvec)[0], zvec[-1] + 0.5*np.diff(zvec)[-1],
    ),
    interpolation='none',
    vmin=-0.4, vmax=+0.4, cmap='RdBu',
    #norm=mpl.colors.LogNorm(),
)
plt.colorbar()
plt.gca().set_aspect('auto')
plt.title(r'$\delta B_x$')
plt.xlabel(r'$t \Omega_\mathrm{ci}$')
plt.ylabel(r'$z / d_\mathrm{i}$')
plt.subplots_adjust(top=0.8, bottom=0.2, left=0.1, right=0.9)
plt.savefig("media/pcai-dbx2d.pdf", dpi=300, bbox_inches='tight')
plt.show()

plt.figure(figsize=(4,2))
plt.plot(time, np.mean(bx**2, axis=1), '.-', label=r'$|\delta {B_x}^2|$')
plt.plot(time, np.mean(by**2, axis=1), '.-', label=r'$|\delta {B_y}^2|$')
#plt.plot(time, np.mean(ex**2, axis=1), '.-', label=r'$|\delta {E_x}^2|$')
#plt.plot(time, np.mean(ey**2, axis=1), '.-', label=r'$|\delta {E_y}^2|$')
plt.plot(time, np.mean(jfx**2, axis=1), '.-', label=r'$|\delta {J_x}^2|$')
plt.plot(time, np.mean(jfy**2, axis=1), '.-', label=r'$|\delta {J_y}^2|$')
plt.yscale('log')
plt.xlabel(r'$t \Omega_\mathrm{ci}$')
plt.legend()
plt.subplots_adjust(top=0.9, bottom=0.2, left=0.1, right=0.9)
plt.savefig("media/pcai-dbsq1d.pdf", dpi=300, bbox_inches='tight')
plt.show()

plt.figure(figsize=(4,2))
plt.plot(time, np.mean(np.abs(bx), axis=1), '.-', label=r'$|\delta B_x|$')
plt.plot(time, np.mean(np.abs(by), axis=1), '.-', label=r'$|\delta B_y|$')
plt.yscale('log')
plt.xlabel(r'$t \Omega_\mathrm{ci}$')
plt.legend()
plt.subplots_adjust(top=0.9, bottom=0.2, left=0.1, right=0.9)
plt.savefig("media/pcai-db1d.pdf", dpi=300, bbox_inches='tight')
plt.show()
