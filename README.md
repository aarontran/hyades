README
======

Hyades is an under-construction hybrid plasma simulation code (kinetic ions,
fluid electrons) implemented in C++.

It is a stripped-down, highly simplified and lean fork of the Hybrid-VPIC and
VPIC simulation codes, keeping the overall framework and design, but omitting
many features.


Developer notes
---------------

Code should follow the [Google style guide](https://google.github.io/styleguide/cppguide.html).
* TODO update class member naming convention to use underscores per
  [https://google.github.io/styleguide/cppguide.html#Variable_Names](here),
  I am already getting bit by naming scheme...
* avoid macros to the largest extent possible

Take programming cues from VPIC, Hybrid-VPIC, TRISTAN-MP, Pegasus, WarpX,
Torch, and other codes on the market.

TODO: consider only implementing shared-memory parallelism and GPU
acceleration?
Make this code work for 1 node ONLY, periodic boundaries ONLY.
MPI communication is a big job.


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


About Hyades
------------

It has several goals:

1. Implement different algorithms, in a modular way, to study and compare their
   performance.  This includes, in rough priority order:

   + Field time-advance algorithm.
   + Different mesh staggering.
   + Varying particle deposit schemes for density and current moments.
   + Communication: 6 versus 26 neighbors
   + Ghost cell count to allow testing of higher-order algorithms.
   + Boundary conditions for fields.

2. Reasonably modern/clean code (by the standards of scientific programmers).

3. Specific: designed for 3D Cartesian grid, one ion species only.

4. Able to use multiple levels of parallelism.

It's possible not all goals will achieved.  This is a project to teach myself
C++ parallel programming and hopefully also get some science done concerning
hybrid plasma algorithms.

