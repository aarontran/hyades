// ----------------------------------------------------------------------------
// Test problem: 1D domain with anisotropic protons to drive cyclotron waves
// matched to Hybrid-VPIC deck "examples/pcai"
// ----------------------------------------------------------------------------

#include <cmath>  // for sqrt

#include "../src/field.h"
//#include "interp.h"
#include "../src/param.h"
#include "../src/particle.h"
//#include "random.h"
//#include "timer.h"

void user_param(param_t* par) {

  // WARNING: >=2 ghost cells are required for
  // * particle move+deposit scheme to not deal with boundary crossings
  // * Ohm's law Laplacian

  par->idumpf =  50;     // fields dump interval, 0.5 Omci^-1
  par->idumpp =2000;     // particles dump interval, 20 Omci^-1
  //par->isort = 20;
  par->ilast = 6000;     // run duration 60 Omci^-1
  par->Lx = (10.5/64.);  // use cubical cells
  par->Ly = (10.5/64.);
  par->Lz = 10.5;
  par->nx = 1;     // number of x cells in domain (not counting ghosts)
  par->ny = 1;
  par->nz = 64;
  par->ng =  2;    // number of ghost cells
  par->seed = 1;   // zero-th particle gets teleported so useful test case

  par->stridep = 100;  // particle output stride
  //par->stridef = 100;  // fields output stride  // not implemented

  par->dt = 0.01;

  par->nppc  = 10000;
  par->npmax = 1000000;

}

void user_initialize(FieldArray* fa, ParticleArray* ions, param_t* par) {

  double vthi = 0.707;  // sqrt(kB*Ti/mi) = sqrt(beta,i/2)
  double TiTe = 1;
  double Tianiso = 3;  // perp/parallel temperature ratio
  int npart = par->nppc*(par->nx*par->ny*par->nz);

  // Hybrid algorithm parameters
  fa->hyb_te_ref_   = vthi*vthi/TiTe;  // Electron temperature kB*Te/(mi*vA0^2) scaled to ION mass
  fa->hyb_ne_ref_   = 1.;              // Electron density ref value for EoS
  fa->hyb_ne_floor_ = 0.05;            // Electron density floor
  fa->hyb_eta_      = 0;//1e-4;        // Resistivity
  fa->hyb_hypereta_ = 0;//1e-4;        // Hyper-resistivity

  fa->nsubcycle_    = 1;               // Number of field-advance subcycles

  // Boundary conditions
  fa->particle_bc_x = 0;  // periodic = 0, reflecting = 1
  fa->particle_bc_y = 0;  // periodic = 0, reflecting = 1 (not implemented for y/z yet)
  fa->particle_bc_z = 0;  // periodic = 0, reflecting = 1 (not implemented for y/z yet)
  fa->field_bc_x = 0;  // periodic = 0, conducting zero E = 1 (not implemented yet)
  fa->field_bc_y = 0;  // periodic = 0, conducting zero E = 1 (not implemented yet)
  fa->field_bc_z = 0;  // periodic = 0, conducting zero E = 1 (not implemented yet)

  // Initialize fields
  fa->mesh_set_all(0.);
  fa->mesh_set_b(0, 0, 1.);

  // Charge per macro ion
  double weight = par->Lx*par->Ly*par->Lz/(par->nppc*par->nx*par->ny*par->nz);

  ions->initialize(npart, weight);
  ions->maxwellian1d(0, npart, "ux", vthi*sqrt(Tianiso), 0);
  ions->maxwellian1d(0, npart, "uy", vthi*sqrt(Tianiso), 0);
  ions->maxwellian1d(0, npart, "uz", vthi,               0);
  ions->uniform(0, npart, 0., par->Lx, 0., par->Ly, 0., par->Lz);

}
