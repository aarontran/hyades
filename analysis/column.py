#!/usr/bin/env python
"""
Inspect column problem
"""

from __future__ import division, print_function

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt

import h5py

fname = "flds.all.hdf5"
dt = 0.01
idumpf = 1000  # beware hardcoded from input parameters
Ly = 6.667
Lz = 6.667

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

yvec = np.arange(0, bx.shape[-2]) * (Ly/(bx.shape[-2]-1))
zvec = np.arange(0, bx.shape[-1]) * (Lz/(bx.shape[-1]-1))

for it in range(bx.shape[0]):

    plt.figure(figsize=(4,2))
    plt.imshow(
        rhof[it,0,:,:].T,
        origin='lower',
        extent=(
            yvec[0] - 0.5*np.diff(yvec)[0], yvec[-1] + 0.5*np.diff(yvec)[-1],
            zvec[0] - 0.5*np.diff(zvec)[0], zvec[-1] + 0.5*np.diff(zvec)[-1],
        ),
        interpolation='none',
        vmin=0, vmax=1,
        cmap='turbo',
        #norm=mpl.colors.LogNorm(),
    )
    plt.colorbar()
    plt.title(r'$t\Omega_\mathrm{ci} =' + '{:.1f}$'.format(it*dt*idumpf))
    plt.xlabel(r'$y / d_\mathrm{i}$')
    plt.ylabel(r'$z / d_\mathrm{i}$')
    #plt.subplots_adjust(top=0.8, bottom=0.2, left=0.1, right=0.9)
    #plt.show()
    plt.savefig(f"media-frames/column_{it:03d}.png", dpi=300, bbox_inches='tight')
    plt.clf()
    plt.close()

    print("plotted", it)
