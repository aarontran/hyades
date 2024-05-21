// ----------------------------------------------------------------------------
// Test problem: 1D expansion of plasma into vacuum
// following Amano+ (2014) Section 4.4 test
// ----------------------------------------------------------------------------

#include <cmath>  // for sqrt

#include "../src/field.h"
#include "../src/param.h"
#include "../src/particle.h"

void user_param(param_t* par) {

  // WARNING: >=2 ghost cells are required for
  // * particle move+deposit scheme to not deal with boundary crossings
  // * Ohm's law Laplacian

  par->idumpf = 1000;     // fields dump interval, 0.5 Omci^-1
  par->idumpp =   0;      // particles dump interval, skip for now
  //par->isort = 20;
  par->ilast = 4800;      // run duration 48 Omci^-1
  par->Lx = 256;
  par->Ly = 0.5;
  par->Lz = 0.5;
  par->nx = 512;  // number of x cells in domain (not counting ghosts)
  par->ny = 1;
  par->nz = 1;
  par->ng = 2;    // number of ghost cells
  par->seed = 1;

  par->stridep = 100;  // particle output stride
  //par->stridef = 100;  // fields output stride  // not implemented

  par->dt = 0.01;

  par->nppc  =    100; // low PPC for prototyping
  par->npmax = 100000;

}

void user_initialize(FieldArray* fa, ParticleArray* ions, param_t* par) {

  double vthi = 0.0707107;  // sqrt(kB*Ti/mi) = sqrt(beta,i/2)
  double TiTe = 1;

  // Hybrid algorithm parameters
  fa->hyb_te_ref_   = vthi*vthi/TiTe;  // Electron temperature kB*Te/(mi*vA0^2) scaled to ION mass
  fa->hyb_ne_ref_   = 1.;              // Electron density ref value for EoS
  fa->hyb_ne_floor_ = 0.02;            // Electron density floor
  fa->hyb_eta_      = 0;               // Resistivity
  fa->hyb_hypereta_ = 0;               // Hyper-resistivity

  fa->nsubcycle_    = 1;               // Number of field-advance subcycles

  // Boundary conditions
  fa->particle_bc_x = 1;  // periodic = 0, reflecting = 1
  fa->particle_bc_y = 0;  // periodic = 0, reflecting = 1 (not implemented for y/z yet)
  fa->particle_bc_z = 0;  // periodic = 0, reflecting = 1 (not implemented for y/z yet)
  fa->field_bc_x = 0;  // periodic = 0, conducting zero E = 1 (not implemented yet)
  fa->field_bc_y = 0;  // periodic = 0, conducting zero E = 1 (not implemented yet)
  fa->field_bc_z = 0;  // periodic = 0, conducting zero E = 1 (not implemented yet)

  printf("***user init got BC setting x y z %d %d %d\n", // TODO TESTING
      fa->particle_bc_x,
      fa->particle_bc_y,
      fa->particle_bc_z);

  // Initialize fields
  fa->mesh_set_all(0.);
  fa->mesh_set_b(1., 0., 0.);

  // Charge per macro ion
  int npart_max = par->nppc*(par->nx*par->ny*par->nz);
  double weight = par->Lx*par->Ly*par->Lz/npart_max;

  float fstep = 0.5;
  int npart = (int)(fstep*npart_max);

  ions->initialize(npart, weight);
  ions->maxwellian1d(0, npart, "ux", vthi, 0);
  ions->maxwellian1d(0, npart, "uy", vthi, 0);
  ions->maxwellian1d(0, npart, "uz", vthi, 0);
  // Only initialize on part of domain to make discontinuity
  // NOTE initialize with +0.5 because of reflecting wall boundary...
  ions->uniform(0, npart, 0.5, fstep*par->Lx, 0., par->Ly, 0., par->Lz);

}
