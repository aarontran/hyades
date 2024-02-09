#include <stdlib.h>  // for malloc
#include <stdio.h>  // warn about buffer overflows

#include "particles.h"
#include "random.h"

#ifndef M_PI
#define M_PI 3.1415926535897932384626
#endif

// Constructor
// Have to specify the rng in constructor init list to start!!
// Otherwise C++ tries to call a default constructor which I haven't
// provided...
// https://stackoverflow.com/q/31488756
ParticleArray::ParticleArray(float q_, float m_, int npmax_, Random rng_) :
  q(q_), m(m_), npmax(npmax_), rng(rng_)
{
  //q = q_;
  //m = m_;
  //npmax = npmax_;
  np = 0;
  p = (particle_t*) malloc( npmax_*sizeof(particle_t) );
  //rng = rng_;
}

// Sort the particle list
void ParticleArray::sort() {
  return;
}

// Initialize particles with zero position and velocity
void ParticleArray::initialize(int count) {
  if (np+count > npmax) {
    printf("WARNING: buffer overflow in particle init!\n");
    return;
  }
  particle_t* p0 = &(p[np]);
  for (int ii=np; ii<np+count; ++ii) {
    p0->x   = 0;
    p0->y   = 0;
    p0->z   = 0;
    p0->ux  = 0;
    p0->uy  = 0;
    p0->uz  = 0;
    p0->w   = 1;
    p0->ind = ii;
    ++p0;
  }
  np = np+count;
}

// Initialize particle velocities as Maxwellian
void ParticleArray::maxwellian(int i0, int i1, float vth, float vdrx,
                               float vdry, float vdrz) {

  if (i0 >= np || i1 <= i0) {
    printf("WARNING: out of bounds in particle Maxwellian init!\n");
    return;
  }
  particle_t* p0 = &(p[i0]);
  for (int ii=i0; ii<i1; ++ii) {
    p0->ux = rng.normal(vdrx, vth);
    p0->uy = rng.normal(vdry, vth);
    p0->uz = rng.normal(vdrz, vth);
    ++p0;
  }
}

// Initialize particle positions uniformly on a grid
void ParticleArray::uniform(int i0, int i1, float x0, float x1, float y0,
                            float y1, float z0, float z1) {

  if (i0 >= np || i1 <= i0) {
    printf("WARNING: out of bounds in particle uniform init!\n");
    return;
  }
  particle_t* p0 = &(p[i0]);
  for (int ii=i0; ii<i1; ++ii) {
    p0->x = x0 + (x1-x0)*rng.drand();
    p0->y = y0 + (y1-y0)*rng.drand();
    p0->z = z0 + (z1-z0)*rng.drand();
    ++p0;
  }
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
