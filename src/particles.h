#ifndef PARTICLES_H
#define PARTICLES_H

#include "fields.h"
#include "random.h"

typedef struct particle {
  float  x,  y,  z;
  float ux, uy, uz;
  float  w;
  int32_t ind;
  // Don't implement offsetting until you need it.
  // Keep the code simple.
  //float dx, dy, dz;   // position in cell coordinates (on [-1,1])
  //int32_t i;          // voxel containing the particle
} particle_t;


class ParticleArray {

  public:
    ParticleArray(float q_, float m_, int npmax_, FieldArray fa_,
                  InterpArray ia_, Random rng_);

    float q;    // charge
    float m;    // mass
    int npmax;  // max number of particles
    int np;     // current number of particles
    particle_t* p0;  // particle array pointer

    void sort();
    void initialize(int count);
    void maxwellian(int i0, int i1, float vth, float vdrx, float vdry, float vdrz);
    void uniform(int i0, int i1, float x0, float x1, float y0, float y1, float z0, float z1);
    float meanq_vx();
    float meanq_vy();
    float meanq_vz();
    float meanq_vsq();

    void move_boris();

  private:
    FieldArray   fa;
    InterpArray  ia;
    Random      rng;

};


#endif  // PARTICLES_H
