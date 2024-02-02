README
======

Hyades is an under-construction hybrid plasma simulation code (kinetic ions,
fluid electrons) implemented in C++.


Developer notes
---------------

Code should follow the [Google style guide](https://google.github.io/styleguide/cppguide.html).

Take programming cues from VPIC, Hybrid-VPIC, TRISTAN-MP, Pegasus, WarpX,
Torch, and other codes on the market.

TODO: consider only implementing shared-memory parallelism and GPU
acceleration?
Make this code work for 1 node ONLY, periodic boundaries ONLY.
MPI communication is a big job.

Convention: typedef structs have `_t` lowercase, classes use CamelCase.


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

