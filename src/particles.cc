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
void ParticleArray::maxwellian(float vth, float vdrx, float vdry, float vdrz, int count, int seed) {
  if (np+count > npmax) {
    printf("WARNING: buffer overflow, aborting maxwellian init!!\n");
    return;
  }
  Random rng = Random(seed);
  particle_t* p0 = &(p[np]);
  for (int ii=np; ii<np+count; ++ii) {
    p0->x   = 0;  // todo distribute on grid
    p0->y   = 0;
    p0->z   = 0;
    p0->ux  = rng.normal(vdrx, vth);
    p0->uy  = rng.normal(vdry, vth);
    p0->uz  = rng.normal(vdrz, vth);
    p0->w   = 1;
    p0->ind = ii;
    ++p0;
  }
  np = np+count;
}

// Compute vx moment on entire domain
float ParticleArray::meanq_vx() {
  particle_t* p0 = p;
  double v = 0;  // High prec to reduce round-off error in summing many terms
  for (int ii=0; ii<np; ++ii) {
    v = v + p0->ux;
    ++p0;
  }
  return v/np;
}

// Compute vy moment on entire domain
float ParticleArray::meanq_vy() {
  particle_t* p0 = p;
  double v = 0;  // High prec to reduce round-off error in summing many terms
  for (int ii=0; ii<np; ++ii) {
    v = v + p0->uy;
    ++p0;
  }
  return v/np;
}

// Compute vz moment on entire domain
float ParticleArray::meanq_vz() {
  particle_t* p0 = p;
  double v = 0;  // High prec to reduce round-off error in summing many terms
  for (int ii=0; ii<np; ++ii) {
    v = v + p0->uz;
    ++p0;
  }
  return v/np;
}

// Compute v^2 moment on entire domain
float ParticleArray::meanq_vsq() {
  particle_t* p0 = p;
  double vsq = 0;  // High prec to reduce round-off error in summing many terms
  for (int ii=0; ii<np; ++ii) {
    vsq = vsq + ((p0->ux)*(p0->ux) + (p0->uy)*(p0->uy) + (p0->uz)*(p0->uz)) ;
    ++p0;
  }
  return vsq/np;
}
