README
======

Hyades is a hybrid plasma simulation code (kinetic ions, fluid electrons).
It's a simpler, stripped-down fork of the
[Hybrid-VPIC](https://github.com/lanl/vpic-kokkos/tree/hybridVPIC)
simulation code developed by
[Le et al. (2023, Physics of Plasmas)](https://doi.org/10.1063/5.0146529).

__Hyades is under construction and not yet working.__

Compared to Hybrid-VPIC,
* Hyades supports _only_ periodic boundary conditions.
* Hyades can use any number of ghost cells on each coordinate axis, in order to
  support higher-order particle shapes.
* Hyades will be initially designed to run on shared-memory architectures,
  targeting AMD Epyc 7763 nodes.
* Hyades will not be MPI parallelized to start, but thread- and pipeline-level
  parallelism (likely starting with OpenMP) or GPU acceleration wmay be tried
  as time permits.

Hyades' initial goal is to implement a variety of hybrid plasma simulation
algorithms to study and compare their performance.  This may include, in rough
priority order:
* Varying particle deposit and field-to-particle interpolation schemes for
  density and current moments, enabled by arbitrary ghost cell count.
* E/B field time-advance algorithms and mesh staggering.
* Boundary conditions for fields and particles.
* Electron inertia

Hyades is implemented in C++, and its source code loosely aims to follow the
[Google style guide](https://google.github.io/styleguide/cppguide.html).

Hyades' software license is probably whatever applies to Hybrid-VPIC and VPIC.


Developer notes
---------------

Code should follow the
* TODO update class member naming convention to use underscores per
  [https://google.github.io/styleguide/cppguide.html#Variable_Names](here),
  I am already getting bit by naming scheme...
* avoid macros to the largest extent possible

Some design principles:
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
