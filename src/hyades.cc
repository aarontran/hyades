#include <assert.h>
#include <stdlib.h>  // for malloc, free
#include <stdio.h>  // for printf
#include <omp.h>

#include <cmath>  // for sqrt

#include <filesystem>

#include "field.h"
#include "interp.h"
#include "param.h"
#include "particle.h"
#include "random.h"
#include "timer.h"

int main(int argc, char* argv[]) {

  printf("Max number of threads %d\n", omp_get_max_threads());

  // --------------------------------------------------------------------------
  // User initialize
  // --------------------------------------------------------------------------

  // Test problem: 1D domain with anisotropic protons to drive cyclotron waves
  // matched to Hybrid-VPIC deck "examples/pcai"

  // WARNING: >=2 ghost cells are required for
  // * particle move+deposit scheme to not deal with boundary crossings
  // * Ohm's law Laplacian

  param_t        par;
  par.idumpf =  50;     // fields dump interval, 0.5 Omci^-1
  par.idumpp =2000;     // particles dump interval, 20 Omci^-1
  //par.isort = 20;
  par.ilast = 6000;     // run duration 60 Omci^-1
  par.Lx = (10.5/64.);  // use cubical cells
  par.Ly = (10.5/64.);
  par.Lz = 10.5;
  par.nx = 1;     // TODO illegal to have nx,ny=1 with nghost=2,
  par.ny = 1;     // I guess the ghost just keep stacking and get the "other" ghosts so it works
                  //but not a good idea, should guard against this...
                  //ALSO current deposit may be broken! -ATr,2024feb25
  par.nz = 64;
  par.ng =  2;    // number of ghost cells
  par.seed = 1;   // zero-th particle gets teleported so useful test case

  par.stridep = 100;  // particle output stride
  //par.stridef = 100;  // fields output stride  // not implemented

  par.nppc  = 10000;
  par.npmax = 1000000;

  Random         rng = Random(par.seed);
  FieldArray      fa = FieldArray(par.nx, par.ny, par.nz, par.ng,  // nx,ny,nz,ng
                                  par.Lx/par.nx,  // hx
                                  par.Ly/par.ny,  // hy
                                  par.Lz/par.nz,  // hz
                                  0.01);  // dt
  InterpArray     ia = InterpArray(fa);
  ParticleArray ions = ParticleArray(1, 1, par.npmax, fa, ia, rng);

  double vthi = 0.707;  // sqrt(kB*Ti/mi) = sqrt(beta,i/2)
  double TiTe = 1;
  double Tianiso = 3;  // perp/parallel temperature ratio
  int npart = par.nppc*(par.nx*par.ny*par.nz);

  // Hybrid algorithm parameters
  fa.hyb_te_ref_   = vthi*vthi/TiTe;  // Electron temperature kB*Te/(mi*vA0^2) scaled to ION mass
  fa.hyb_ne_ref_   = 1.;              // Electron density ref value for EoS
  fa.hyb_ne_floor_ = 0.05;            // Electron density floor
  fa.hyb_eta_      = 0;//1e-4;        // Resistivity
  fa.hyb_hypereta_ = 0;//1e-4;        // Hyper-resistivity

  int field_nsub   = 1;               // Number of field-advance subcycles

  fa.mesh_set_all(0.);
  fa.mesh_set_b(0, 0, 1.);

  // Charge per macro ion
  double weight = par.Lx*par.Ly*par.Lz/(par.nppc*par.nx*par.ny*par.nz);

  ions.initialize(npart, weight);
  ions.maxwellian1d(0, npart, "ux", vthi*sqrt(Tianiso), 0);
  ions.maxwellian1d(0, npart, "uy", vthi*sqrt(Tianiso), 0);
  ions.maxwellian1d(0, npart, "uz", vthi,               0);
  ions.uniform(0, npart, 0., par.Lx, 0., par.Ly, 0., par.Lz);

  // --------------------------------------------------------------------------
  // Initialize
  // --------------------------------------------------------------------------

  TimerArray clock = TimerArray(100);
  clock.add("total");  // total and init are large blocks
  clock.add("init");
  clock.add("evolve");  // close timing of evolution loop
  clock.add("clearj");
  clock.add("interp");
  clock.add("move");
  clock.add("dep");
  clock.add("bdry");
  clock.add("depgh");
  clock.add("field");

  clock.tic("total");
  clock.tic("init");
  //clock.add("finalize");  // not implemented

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

  std::filesystem::create_directory("output");
  if (par.idumpf > 0) {
    fa.dump(0, "output/flds.%d.hdf5", false);
  }
  if (par.idumpp > 0) {
    ions.dump(0, "output/prtl.%d.hdf5", par.stridep);
  }

  clock.toc("init");
  printf("Initialized in %g seconds.\n", clock.flush("init"));
  printf("Evolving simulation.\n");

  // --------------------------------------------------------------------------
  // Evolve
  // --------------------------------------------------------------------------

  int step = 0;

  while (step < par.ilast) {

    clock.tic("evolve");

    //if (step % par.isort == 0) ions.sort();  // not implemented

    // prepare for particle advance
    clock.tic("clearj");
    fa.mesh_set_jrho0();            // copy old j/rho into j0/rho0
    fa.mesh_set_jrho(0,0,0,0);
    clock.toc("clearj");

    clock.tic("interp");
    //fa.ghost_copy_eb();           // update E/B in ghosts -- already done in field advance
    ia.update();                    // compute E/B field interpolation coeffs
    clock.toc("interp");

    clock.tic("move");
    ions.move();                    // r,v advanced on ghosts
    clock.toc("move");

    clock.tic("dep");
    ions.deposit(1);                // j   advanced on ghosts
    clock.toc("dep");

    clock.tic("bdry");
    ions.boundary_teleport();       // r,v advanced
    clock.toc("bdry");

    clock.tic("depgh");
    fa.ghost_deposit_jrho();
    fa.ghost_copy_jrho();           // j   advanced
    //fa.smooth_jrho();             // not implemented
    clock.toc("depgh");

    clock.tic("field");
    for (int isub=0; isub<field_nsub; ++isub) {
      fa.advance_eb_rk4_ctrmesh(isub, field_nsub);  // E,B advanced on live+ghost
    }
    //fa.smooth_eb();               // not implemented
    clock.toc("field");

    clock.toc("evolve");

    step++;

    if (par.idumpf > 0 && step % par.idumpf == 0) {
      fa.dump(step, "output/flds.%d.hdf5", false);
    }
    if (par.idumpp > 0 && step % par.idumpp == 0) {
      ions.dump(step, "output/prtls.%d.hdf5", par.stridep);
    }

    if (step % 10 == 0) {
      printf("Step %d of %d elapsed %g seconds.\n", step, par.ilast, clock.read("total"));
      printf("    clearj %.2e\n", clock.flush("clearj"));
      printf("    interp %.2e\n", clock.flush("interp"));
      printf("      move %.2e\n", clock.flush("move"));
      printf("       dep %.2e\n", clock.flush("dep"));
      printf("      bdry %.2e\n", clock.flush("bdry"));
      printf("     depgh %.2e\n", clock.flush("depgh"));
      printf("     field %.2e\n", clock.flush("field"));
      printf("     total %.2e\n", clock.flush("evolve"));
    }

  }

  // At the end of each loop iteration
  // ... particle r advanced from t=n to t=n+1
  // ... particle v advanced from t=n-1/2 to t=n+1/2
  // ... fields j,rho advanced from t=n-1/2 to t=n+1/2
  // ... fields E,B advanced from t=n to t=n+1

  // --------------------------------------------------------------------------
  // Finalize
  // --------------------------------------------------------------------------

  free(ions.p0);
  free(ia.ic0);
  free(fa.f0);
  free(fa.ivoxels_ghost);
  free(fa.ivoxels_ghsrc);

  clock.toc("total");

  // in case we didn't end on a step with timer printout
  clock.flush_all();

  printf("All done in %g seconds.\n", clock.read_total("total"));
  printf("Initializing took %g seconds.\n", clock.read_total("init"));
  printf("Evolving took %g seconds.\n", clock.read_total("evolve"));
  printf("    clearj %.2e\n", clock.read_total("clearj"));
  printf("    interp %.2e\n", clock.read_total("interp"));
  printf("      move %.2e\n", clock.read_total("move"));
  printf("       dep %.2e\n", clock.read_total("dep"));
  printf("      bdry %.2e\n", clock.read_total("bdry"));
  printf("     depgh %.2e\n", clock.read_total("depgh"));
  printf("     field %.2e\n", clock.read_total("field"));

  free(clock.t0);

  return 0;
}
