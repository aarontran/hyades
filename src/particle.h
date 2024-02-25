#ifndef PARTICLE_H
#define PARTICLE_H

#include "hdf5.h"

#include "field.h"
#include "random.h"

typedef struct particle {
  // If you change struct layout, also update ParticleArray subroutines
  // dump(...), pseek_fkey(...), pseek_ikey(...), and any other code that
  // accesses struct members by explicit index or name
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
    ParticleArray(float qsp_, float msp_, int npmax_, FieldArray fa_,
                  InterpArray ia_, Random rng_);

    float qsp;  // charge
    float msp;  // mass
    int npmax;  // max number of particles
    int np;     // current number of particles
    particle_t* p0;  // particle array pointer

    // --------------------------------------------------
    // Low-level methods to get/set particle attribute values
    float*   pseek_fkey(const char* name, particle_t* p);
    int32_t* pseek_ikey(const char* name, particle_t* p);

    // --------------------------------------------------
    // High-level methods for top-level hybrid algorithm
    void sort();
    void initialize(int count, float weight);
    void maxwellian(int i0, int i1, float vth, float vdrx, float vdry, float vdrz);
    void uniform(int i0, int i1, float x0, float x1, float y0, float y1, float z0, float z1);
    float meanq_vx();
    float meanq_vy();
    float meanq_vz();
    float meanq_vsq();

    void move();
    void move_uncenter();
    void deposit(int unwind);
    void boundary_teleport();

    // High-level dump methods
    void dump(int step, int stride);
    void hputf(hid_t file_id, const char* attr_name, int stride);
    void hputi(hid_t file_id, const char* attr_name, int stride);

  private:
    FieldArray   fa;
    InterpArray  ia;
    Random      rng;

};


#endif  // PARTICLE_H
