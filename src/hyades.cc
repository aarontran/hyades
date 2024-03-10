#include <assert.h>
#include <stdlib.h>  // for malloc, free
#include <stdio.h>  // for printf
#include <sys/resource.h>  // for getrusage
#include <omp.h>

#include <filesystem>

#include "field.h"
#include "interp.h"
#include "param.h"
#include "particle.h"
#include "random.h"
#include "timer.h"
#include "user.h"

int main(int argc, char* argv[]) {

  printf("Max number of threads %d\n", omp_get_max_threads());

  // --------------------------------------------------------------------------
  // User initialize
  // --------------------------------------------------------------------------

  param_t par;
  user_param(&par);

  FieldArray      fa = FieldArray(par.nx, par.ny, par.nz, par.ng,
                                  par.Lx/par.nx, par.Ly/par.ny, par.Lz/par.nz,
                                  par.dt);
  InterpArray     ia = InterpArray(fa);
  Random         rng = Random(par.seed);
  ParticleArray ions = ParticleArray(1, 1, par.npmax, fa, ia, rng);

  user_initialize(&fa, &ions, &par);

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
  {
    struct rusage result;
    getrusage(RUSAGE_SELF, &result);
    // Mac OS X maxrss in bytes, Linux in kilobytes; see "man getrusage".
#if defined(__APPLE__) && defined(__MACH__)
    printf("Memory (maxRSS) %ld MB\n", result.ru_maxrss/1000/1000);
#else
    printf("Memory (maxRSS) %ld MB\n", result.ru_maxrss/1000);
#endif
  }
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
    fa.advance_eb_rk4_ctrmesh();    // E,B advanced on live+ghost
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

    if (step % 100 == 0) {
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
  printf("Evolution summary.\n");
  printf("    clearj %.2e\n", clock.read_total("clearj"));
  printf("    interp %.2e\n", clock.read_total("interp"));
  printf("      move %.2e\n", clock.read_total("move"));
  printf("       dep %.2e\n", clock.read_total("dep"));
  printf("      bdry %.2e\n", clock.read_total("bdry"));
  printf("     depgh %.2e\n", clock.read_total("depgh"));
  printf("     field %.2e\n", clock.read_total("field"));
  printf("     total %.2e\n", clock.read_total("evolve"));

  free(clock.t0);

  return 0;
}
