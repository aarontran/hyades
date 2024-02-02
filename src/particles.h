#ifndef PARTICLES_H
#define PARTICLES_H

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

  private:

  public:
    ParticleArray(float q, float m, int npmax);

    float q;    // charge
    float m;    // mass
    int npmax;  // max number of particles
    int np;     // current number of particles
    particle_t* p;

    void sort();
    void maxwellian(float vth, int count, int rngseed);

};


#endif  // PARTICLES_H
