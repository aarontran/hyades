#include <assert.h>
#include <stdlib.h>  // for malloc, free
#include <stdio.h>  // for printf
#include <omp.h>

#include <filesystem>

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
  par.idump = 20;
  //par.isort = 20;
  par.ilast = 100;
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

  // At simulation start, the user provides:
  //     particle r at t=0
  //     particle v at t=0
  //     field B at t=0
  //
  // To begin the simulation, we need:
  //     particle r at t=0
  //     particle v at t=-1/2
  //     field B at t=0
  //     field E at t=0
  //     fields j,rho at t=-1/2
  //
  // Initialize E at t=0 and "unwind" the particles to t=-1/2 as follows.
  // Step 1: deposit ion moments at t=0 into both current j/rho and old j/rho
  // Step 2: advance E with frac=0.  This nominally averages old/new ion
  //         moments stored at t=+/-1/2; "cheat" using t=0 for both old/new.
  // Step 3: unwind particle velocities to t=-1/2, keep t=0 positions
  // Step 4: deposit particle moments at t=-1/2

  // Step 1.
  fa.mesh_set_jrho(0., 0., 0., 0.);
  ions.deposit(0);                  // arg=0 to deposit using r at t=0
  fa.ghost_deposit_jrho();
  fa.ghost_copy_jrho();             // j,rho at t=0
  fa.mesh_set_jrho0();              // old j,rho at t=0

  // Step 2.
  fa.ghost_copy_b();
  fa.advance_e_ctrmesh(0.);         // E at t=0 using old=new ion moments

  // Step 3.
  ions.move_uncenter();             // v at t=-1/2, keeping r at t=0

  // Step 4.
  fa.mesh_set_jrho(0., 0., 0., 0.);
  ions.deposit(1);                  // arg=1 to deposit using r at t=-1/2
  fa.ghost_deposit_jrho();
  fa.ghost_copy_jrho();             // j,rho at t=-1/2

  //Simulation* sim = Simulation(...);  // may want for checkpoints eventually
  //sim->par  = par
  //sim->fa   = fa
  //sim->ions = ions

  if (par.idump > 0) {
    std::filesystem::create_directory("output");
    fa.dump(0);
    //ions.dump(0, "output");  // not implemented
  }

  {
    //field_t* fv = fa.voxel(2,2,2); // Useful to check deposit at ghosts
    field_t* fv = fa.voxel(5,5,5);
    printf("init field bx %f by %f bz %f\n", fv->bx, fv->by, fv->bz);
    printf("init field ex %f ey %f ez %f\n", fv->ex, fv->ey, fv->ez);
    printf("init field jfx %f jfy %f jfz %f rhof %f\n", fv->jfx, fv->jfy, fv->jfz, fv->rhof);
  }

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

    ions.move();                    // r,v advanced on ghosts
    ions.deposit(1);                // j   advanced on ghosts
    ions.boundary_teleport();       // r,v advanced
    fa.ghost_deposit_jrho();
    fa.ghost_copy_jrho();           // j advanced
    //fa.smooth_jrho();             // not implemented

    for (int isub=0; isub<field_nsub; ++isub) {
      fa.advance_eb_rk4_ctrmesh(isub, field_nsub);  // E,B advanced on live+ghost
    }
    //fa.smooth_eb();               // not implemented

    step++;

    if (par.idump > 0 && step % par.idump == 0) {
      fa.dump(step);
      //ions.dump(step, "output");  // not implemented
    }

  }

  // At the end of each loop iteration
  // ... particle r advanced from t=n to t=n+1
  // ... particle v advanced from t=n-1/2 to t=n+1/2
  // ... fields j,rho advanced from t=n-1/2 to t=n+1/2
  // ... fields E,B advanced from t=n to t=n+1

  // --------------------------------------------------------------------------
  // Finalize

  {
    //printf("final interp array dbxdx %f\n", ia.voxel(3,3,3)->dbxdx);
    //field_t* fv = fa.voxel(2,2,2); // Useful to check deposit at ghosts
    field_t* fv = fa.voxel(5,5,5);
    printf("\n");
    printf("final field bx %f by %f bz %f\n", fv->bx, fv->by, fv->bz);
    printf("final field ex %f ey %f ez %f\n", fv->ex, fv->ey, fv->ez);
    printf("final field jfx %f jfy %f jfz %f rhof %f\n", fv->jfx, fv->jfy, fv->jfz, fv->rhof);
  }

  free(ions.p0);
  free(ia.ic0);
  free(fa.f0);
  free(fa.ivoxels_ghost);
  free(fa.ivoxels_ghsrc);

  return 0;
}
