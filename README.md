README
======

Hyades is a hybrid plasma simulation code (kinetic ions, fluid electrons).
It's a simpler, stripped-down fork of the
[Hybrid-VPIC](https://github.com/lanl/vpic-kokkos/tree/hybridVPIC)
simulation code developed by
[Le et al. (2023)](https://doi.org/10.1063/5.0146529),
which uses an explicit time-stepping algorithm similar to the H3D code
([Karimabadi et al., 2014;](https://doi.org/10.1063/1.4882875)
 [Le et al., 2016](https://doi.org/10.1063/1.4943893)).

__Hyades is under construction and may have bugs.  It has been tested and shown
to work on one case (1D proton cyclotron anisotropy instability problem).__

Compared to Hybrid-VPIC,
* Hyades supports _only_ periodic boundary conditions.
* Hyades uses >=2 ghost cells on all coordinate axes.
* Hyades is not parallelized at all.

Hyades takes 3e-8 to 5e-8 seconds/particle/step using clang or Intel compilers
with typical optimization flags (see Makefile) on reasonably modern (for 2024)
hardware.  Performance was tested on a 10^3 domain with 1000 particles/cell for
100 steps and 1 subcycle/field advance using commit `cdebbc6b`.

Usage
-----

You will need a C++ compiler with HDF5 libraries.
Edit the Makefile and do:

    make -j
    ./hyades

to perform a simulation, which dumps HDF5 field and particle data to `output/`
in the current working directory.  Analysis code is provided in:

    analysis/
        Makefile
        collate.py
        pcai.py

which uses Python, HDF5, and h5py.

The current input parameters (as of 2024 Feb 25) will setup a 1D simulation of
B-field aligned proton cyclotron anisotropy instability (PCAI) that takes a few
minutes to run on one core.

Developer notes
---------------

Hyades is implemented in C++, and its source code loosely aims to follow the
[Google style guide](https://google.github.io/styleguide/cppguide.html).

Hyades' software license is probably whatever applies to Hybrid-VPIC and VPIC.

Some design principles:
* Avoid macros as much as possible
* Major data structures
  + field.h    -> class FieldArray    -> struct field
  + interp.h   -> class InterpArray   -> struct interp
  + particle.h -> class ParticleArray -> struct particle
* Auxiliary data structures
  + param.h    -> struct param
  + random.h   -> class Random
* typedef structs have `_t` lowercase, classes use CamelCase.
* all filenames singular, except for "hyades.cc"
* interpolator holds coefficients + stencils for how to map fields to
  particles, so that we can easily swap in/out different schemes.
* Pointer type declarations cling to the data type, i.e.,

        float* ux = 1;      // yes
        float *ux = 1;      // avoid
        float * ux = 1;     // avoid
