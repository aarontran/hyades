#include <assert.h>
#include <stdlib.h>  // for malloc, free
#include <stdio.h>  // for printf
#include <omp.h>

#include "field.h"
#include "interp.h"
#include "param.h"
#include "particle.h"
#include "random.h"

int main(int argc, char* argv[]) {

  printf("Max number of threads %d\n", omp_get_max_threads());

  // --------------------------------------------------------------------------
  // Initialization

  // WARNING: >=2 ghost cells are required for
  // * particle move+deposit scheme to not deal with boundary crossings
  // * Ohm's law Laplacian

  param_t        par;
  //par.idump = 100;
  //par.isort = 20;
  par.ilast = 10;
  par.Lx = 10;
  par.Ly = 10;
  par.Lz = 10;
  par.nx = 10;
  par.ny = 10;
  par.nz = 10;
  par.ng =  2;    // number of ghost cells
  par.seed = 1;   // zero-th particle gets teleported so useful test case

  par.nppc  = 100;
  par.npmax = 200000;

  Random         rng = Random(par.seed);
  FieldArray      fa = FieldArray(par.nx, par.ny, par.nz, par.ng,  // nx,ny,nz,ng
                                  par.Lx/par.nx,  // hx
                                  par.Ly/par.ny,  // hy
                                  par.Lz/par.nz,  // hz
                                  0.001);  // dt
  InterpArray     ia = InterpArray(fa);
  ParticleArray ions = ParticleArray(1, 1, par.npmax, fa, ia, rng);

  double vthi = 0.404;
  double TiTe = 8;
  int npart = par.nppc*(par.nx*par.ny*par.nz);

  // Hybrid algorithm parameters
  fa.hyb_te_ref_   = vthi*vthi/TiTe;  // Electron temperature kB*Te/(mi*vA0^2) scaled to ION mass
  fa.hyb_ne_ref_   = 1.;              // Electron density ref value for EoS
  fa.hyb_ne_floor_ = 0.05;            // Electron density floor
  fa.hyb_eta_      = 1e-4;            // Resistivity
  fa.hyb_hypereta_ = 1e-4;            // Hyper-resistivity

  int field_nsub   = 10;              // Number of field-advance subcycles

  fa.mesh_set_all(0.);
  fa.mesh_set_b(0, 0, 1);
  fa.mesh_set_e(0, 0, 0);

  ions.initialize(npart);
  ions.maxwellian(0, npart, vthi, 10., 0, 0);
  ions.uniform(0, npart, 0., par.Lx, 0., par.Ly, 0., par.Lz);

  // Set initial E/B values on grid
  // This requires deposition of ion moments
  //field_advance(fa, ions);
  fa.mesh_set_jrho(0., 0., 0., 0.);     // TODO not correct, need to deposit
  fa.mesh_set_jrho0();                  // ion moments --ATr,2024feb24

  //Simulation* sim = Simulation(...);  // may want for checkpoints eventually
  //sim->par  = par
  //sim->fa   = fa
  //sim->ions = ions

  // --------------------------------------------------------------------------
  // Evolution

  int step = 0;

  while (step < par.ilast) {

    //if (step % par.isort == 0) ions.sort();  // not implemented

    // prepare for particle advance
    fa.mesh_set_jrho0();              // copy old j/rho into j0/rho0
    fa.mesh_set_jrho(0., 0., 0., 0.);
    fa.ghost_copy_eb();             // update E/B in ghosts
    ia.update();                    // compute E/B field interpolation coeffs

    ions.move_deposit();            // r,v advanced
    fa.ghost_deposit_jrho();
    fa.ghost_copy_jrho();           // j advanced on live+ghost
    //fa.smooth_jrho();             // not implemented

    for (int isub=0; isub<field_nsub; ++isub) {
      fa.advance_eb_rk4_ctrmesh(isub, field_nsub);  // E,B advanced on live+ghost
    }
    //fa.smooth_eb();               // not implemented

    step++;

    //if (step % par->idump == 0) diagnostics();  // not implemented

  }

  // At the end of each loop iteration
  // ... particle r advanced from t=n to t=n+1
  // ... particle v advanced from t=n-1/2 to t=n+1/2
  // ... fields j,rho advanced from t=n-1/2 to t=n+1/2
  // ... fields E,B advanced from t=n to t=n+1

  // --------------------------------------------------------------------------
  // Finalize

  printf("interp array dbxdx %f\n", ia.voxel(3,3,3)->dbxdx);
  printf("field array bx %f\n", fa.voxel(5,5,5)->bx);
  printf("field array jfx %f\n", fa.voxel(2,2,2)->jfx);
  printf("field array jfy %f\n", fa.voxel(2,2,2)->jfy);
  printf("field array jfz %f\n", fa.voxel(2,2,2)->jfz);
  printf("field array rhof %f\n",fa.voxel(2,2,2)->rhof);

  free(ions.p0);
  free(ia.ic0);
  free(fa.f0);
  free(fa.ivoxels_ghost);
  free(fa.ivoxels_ghsrc);

  return 0;
}
