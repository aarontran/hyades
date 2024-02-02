#include <stdlib.h>  // for malloc
#include <stdio.h>  // warn about buffer overflows

#include "particles.h"
#include "random.h"

#ifndef M_PI
#define M_PI 3.1415926535897932384626
#endif

// Constructor
ParticleArray::ParticleArray(float q_, float m_, int npmax_) {
  p = (particle_t*) malloc( npmax_*sizeof(particle_t) );
  q = q_;
  m = m_;
  npmax = npmax_;
  np = 0;
}

// Sort the particle list
void ParticleArray::sort() {
  return;
}

// Sort the particle list
void ParticleArray::maxwellian(float vth, int count, int seed) {
  if (count > npmax) {
    printf("WARNING: buffer overflow, aborting maxwellian init!!\n");
    return;
  }
  Random rng = Random(seed);
  particle_t* p0 = p;
  for (int ii=0; ii<count; ++ii) {
    p0->ux = rng.normal(0, vth);
    p0->uy = rng.normal(0, vth);
    p0->uz = rng.normal(0, vth);
    ++p0;
  }
}
