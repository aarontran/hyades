// ----------------------------------------------------------------------------
// Test problem: 2D plasma column cross section mimicking WHAM device
// ----------------------------------------------------------------------------

#include <cmath>  // for sqrt

#include "../src/field.h"
//#include "interp.h"
#include "../src/param.h"
#include "../src/particle.h"
#include "../src/random.h"

void user_param(param_t* par) {

  // WARNING: >=2 ghost cells are required for
  // * particle move+deposit scheme to not deal with boundary crossings
  // * Ohm's law Laplacian

  par->idumpf = 1000;     // fields dump interval, 0.5 Omci^-1
  par->idumpp =   0;      // particles dump interval, skip for now
  //par->isort = 20;
  par->ilast = 33000;     // run duration 330 Omci^-1
  par->Lx = (6.667/96);
  par->Ly = 6.667;
  par->Lz = 6.667;
  par->nx = 1;     // number of x cells in domain (not counting ghosts)
  par->ny = 96;
  par->nz = 96;
  par->ng =  2;    // number of ghost cells
  par->seed = 1;   // zero-th particle gets teleported so useful test case

  par->stridep = 100;  // particle output stride
  //par->stridef = 100;  // fields output stride  // not implemented

  par->dt = 0.01;

  par->nppc  = 200;
  par->npmax = 2000000;

}

void user_initialize(FieldArray* fa, ParticleArray* ions, param_t* par) {

  double vthi = 0.404145;  // sqrt(kB*Ti/mi) = sqrt(beta,i/2)
  double TiTe = 8;

  double rLim = 1.7;  // units of d_i

  // Hybrid algorithm parameters
  fa->hyb_te_ref_   = vthi*vthi/TiTe;  // Electron temperature kB*Te/(mi*vA0^2) scaled to ION mass
  fa->hyb_ne_ref_   = 1.;              // Electron density ref value for EoS
  fa->hyb_ne_floor_ = 0.05;            // Electron density floor
  fa->hyb_eta_      = 1e-4;            // Resistivity
  fa->hyb_hypereta_ = 1e-4;            // Hyper-resistivity

  fa->nsubcycle_    = 50;              // Number of field-advance subcycles

  // Boundary conditions
  fa->particle_bc_x = 0;  // periodic = 0, reflecting = 1
  fa->particle_bc_y = 0;  // periodic = 0, reflecting = 1 (not implemented for y/z yet)
  fa->particle_bc_z = 0;  // periodic = 0, reflecting = 1 (not implemented for y/z yet)
  fa->field_bc_x = 0;  // periodic = 0, conducting zero E = 1
  fa->field_bc_y = 0;  // periodic = 0, conducting zero E = 1 (not implemented yet)
  fa->field_bc_z = 0;  // periodic = 0, conducting zero E = 1 (not implemented yet)

  // Initialize fields
  fa->mesh_set_all(0.);
  fa->mesh_set_b(1., 0., 0.);

  // Charge per macro ion
  int npart_max = par->nppc*(par->nx*par->ny*par->nz);
  double weight = par->Lx*par->Ly*par->Lz/npart_max;

  // just grab another RNG kinda silly but whatever
  Random rng = Random(par->seed);

  int np = 0;
  int np_attempt = 0;
  while (np_attempt < npart_max) {

    // global real coordinates
    double x = -0.5*par->Lx + par->Lx * rng.drand();
    double y = -0.5*par->Ly + par->Ly * rng.drand();
    double z = -0.5*par->Lz + par->Lz * rng.drand();

    double rho = sqrt(y*y + z*z);

    double f = 0;
    if (rho < rLim) {
      f = cos( (rho*rho*rho*rho)/(rLim*rLim*rLim*rLim) * M_PI/2. );
    }
    double coin = rng.drand();
    if (coin < f) {
      // make a new particle
      ions->initialize(1, weight);
      ions->maxwellian1d(np, np+1, "ux", vthi, 0);
      ions->maxwellian1d(np, np+1, "uy", vthi, 0);
      ions->maxwellian1d(np, np+1, "uz", vthi, 0);

      particle_t* pnew = &(ions->p0[np]);
      // adjust to current grid layout 0 to Lz, 0 to Ly, etc..
      // convert to cell-centered grid coordinates
      pnew->x = (x + 0.5*par->Lx)/fa->hx + fa->ng - 0.5;
      pnew->y = (y + 0.5*par->Ly)/fa->hy + fa->ng - 0.5;
      pnew->z = (z + 0.5*par->Lz)/fa->hz + fa->ng - 0.5;

      ++np;
    }

    ++np_attempt;
  }
  printf("Created %d particles of %d attempts\n", np, np_attempt);

  //ions->initialize(npart, weight);
  //ions->maxwellian1d(0, npart, "ux", vthi, 0);
  //ions->maxwellian1d(0, npart, "uy", vthi, 0);
  //ions->maxwellian1d(0, npart, "uz", vthi, 0);
  //ions->uniform(0, npart, 0., par->Lx, 0., par->Ly, 0., par->Lz);

}
