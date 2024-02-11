#include <assert.h>
#include <stdlib.h>  // for malloc
#include <stdio.h>  // warn about buffer overflows

#include "field.h"
#include "interp.h"
#include "particle.h"
#include "random.h"

#ifndef M_PI
#define M_PI 3.1415926535897932384626
#endif

// Constructor
// Have to specify the rng in constructor init list to start!!
// Otherwise C++ tries to call a default constructor which I haven't
// provided...
// https://stackoverflow.com/q/31488756
ParticleArray::ParticleArray(float qsp_, float msp_, int npmax_, FieldArray fa_,
                             InterpArray ia_, Random rng_
):
  qsp(qsp_), msp(msp_), npmax(npmax_), fa(fa_), ia(ia_), rng(rng_)
{
  //q = q_;
  //m = m_;
  //npmax = npmax_;
  np = 0;
  p0 = (particle_t*) malloc( npmax_*sizeof(particle_t) );
  //rng = rng_;
}

// Sort the particle list
void ParticleArray::sort() {
  return;
}

// Initialize particles with zero position and velocity
void ParticleArray::initialize(int count) {

  if (np+count > npmax) {
    printf("ERROR: particle overflow requested %d max %d\n", np+count, npmax);
    //return;
  }
  assert(np+count <= npmax);

  particle_t* p = &(p0[np]);
  for (int ii=np; ii<np+count; ++ii) {
    p->x   = 0;
    p->y   = 0;
    p->z   = 0;
    p->ux  = 0;
    p->uy  = 0;
    p->uz  = 0;
    p->w   = 1;
    p->ind = ii;
    ++p;
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
  particle_t* p = &(p0[i0]);
  for (int ii=i0; ii<i1; ++ii) {
    p->ux = rng.normal(vdrx, vth);
    p->uy = rng.normal(vdry, vth);
    p->uz = rng.normal(vdrz, vth);
    ++p;
  }
}

// Initialize particle positions uniformly on a grid
void ParticleArray::uniform(int i0, int i1, float x0, float x1, float y0,
                            float y1, float z0, float z1) {

  if (i0 >= np || i1 <= i0) {
    printf("WARNING: out of bounds in particle uniform init!\n");
    return;
  }
  particle_t* p = &(p0[i0]);
  for (int ii=i0; ii<i1; ++ii) {
    // global real coordinates
    p->x = x0 + (x1-x0)*rng.drand();
    p->y = y0 + (y1-y0)*rng.drand();
    p->z = z0 + (z1-z0)*rng.drand();
    // convert to cell-centered grid coordinates
    p->x = (p->x)/fa.hx + fa.ng - 0.5;
    p->y = (p->y)/fa.hy + fa.ng - 0.5;
    p->z = (p->z)/fa.hz + fa.ng - 0.5;
    ++p;
  }
}

// Compute vx moment on entire domain
float ParticleArray::meanq_vx() {
  particle_t* p = p0;
  double v = 0;  // High prec to reduce round-off error in summing many terms
  for (int ii=0; ii<np; ++ii) {
    v = v + p->ux;
    ++p;
  }
  return v/np;
}

// Compute vy moment on entire domain
float ParticleArray::meanq_vy() {
  particle_t* p = p0;
  double v = 0;  // High prec to reduce round-off error in summing many terms
  for (int ii=0; ii<np; ++ii) {
    v = v + p->uy;
    ++p;
  }
  return v/np;
}

// Compute vz moment on entire domain
float ParticleArray::meanq_vz() {
  particle_t* p = p0;
  double v = 0;  // High prec to reduce round-off error in summing many terms
  for (int ii=0; ii<np; ++ii) {
    v = v + p->uz;
    ++p;
  }
  return v/np;
}

// Compute v^2 moment on entire domain
float ParticleArray::meanq_vsq() {
  particle_t* p = p0;
  double vsq = 0;  // High prec to reduce round-off error in summing many terms
  for (int ii=0; ii<np; ++ii) {
    vsq = vsq + ((p->ux)*(p->ux) + (p->uy)*(p->uy) + (p->uz)*(p->uz)) ;
    ++p;
  }
  return vsq/np;
}
