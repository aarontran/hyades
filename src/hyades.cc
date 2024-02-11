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
  // Initialization - physics

  param_t        par;  // allocated on stack (not heap) for now...
  //par.idump = 100;
  //par.isort = 20;
  par.ilast = 40;
  par.Lx = 10;
  par.Ly = 10;
  par.Lz = 10;
  //par.nx = 10;
  //par.ny = 10;
  //par.nz = 10;
  par.nx =  5;  // for ghost cell testing
  par.ny = 10;
  par.nz = 20;
  par.ng =  2;  // number of ghost cells
  par.seed = 1;  // zero-th particle gets teleported so useful test case

  par.nppc  = 10;  // total 10,000 particles...
  par.npmax = 20000;


  // WARNING: my particle deposit scheme requires 2 ghost cells
  // so that move+deposit can all be done without any particle aliasing
  // across boundaries

  Random         rng = Random(par.seed);
  FieldArray      fa = FieldArray(par.nx, par.ny, par.nz, par.ng,  // nx,ny,nz,ng
                                  par.Lx/par.nx,  // hx
                                  par.Ly/par.ny,  // hy
                                  par.Lz/par.nz,  // hz
                                  0.001);  // dt
  InterpArray     ia = InterpArray(fa);
  ParticleArray ions = ParticleArray(1, 1, par.npmax, fa, ia, rng);

  double vth = 0.404;
  int npart = par.nppc*(par.nx*par.ny*par.nz);

  fa.mesh_set_all(0.);
  fa.mesh_set_b(1, 1, 1);
  fa.mesh_set_e(0, 0, 0);

  ions.initialize(npart);  // index 0 to 9
  ions.maxwellian(0, npart, vth, 0, 0, 0);
  ions.uniform(0, npart, 0., par.Lx, 0., par.Ly, 0., par.Lz);

  // Set initial E/B values on grid
  // This requires deposition of ion moments
  //field_advance(fa, ions);

  //Simulation* sim = Simulation(...);  // may want for checkpoints eventually
  //sim->par  = par
  //sim->fa   = fa
  //sim->ions = ions

  // --------------------------------------------------------------------------
  // Print basic information about data structures
  // and show how the user can interact with them.

  // Example: field voxel pointers in 
  //printf("voxel(0,0,0) is %p\n", fa.voxel(0,0,0));
  //printf("voxel(1,0,0) is %p\n", fa.voxel(1,0,0));
  //printf("voxel(2,0,0) is %p\n", fa.voxel(2,0,0));

  // Example: field access via voxel indexing
  //printf("bx(0,2,2) is %f\n", fa.voxel(0,2,2)->bx);
  //printf("bx(1,2,2) is %f\n", fa.voxel(1,2,2)->bx);
  //printf("bx(2,2,2) is %f\n", fa.voxel(2,2,2)->bx);
  //printf("bx(11,11,11) is %f\n", fa.voxel(11,11,11)->bx);

  // Example: field access via pointer arithmetic
  //field_t* f = fa.voxel(0,0,0);
  //printf("bx(0,0,0) is %f\n", f->bx);
  //printf("bx(1,0,0) is %f\n", (f+1)->bx);
  //printf("bx(2,0,0) is %f\n", (f+2)->bx);

  // Example: interpolation coefficient access
  //interp_t* ic = ia.voxel(2,2,2);
  //printf("interp at (2,2,2) got bx %f\n",    ic->bx);
  //printf("interp at (2,2,2) got dbxdx %f\n", ic->dbxdx);
  //printf("interp at (2,2,2) got dbxdy %f\n", ic->dbxdy);
  //printf("interp at (2,2,2) got dbxdz %f\n", ic->dbxdz);
  //printf("interp at (2,2,2) got d2bxdx %f\n", ic->d2bxdx);
  //printf("interp at (2,2,2) got d2bxdy %f\n", ic->d2bxdy);
  //printf("interp at (2,2,2) got d2bxdz %f\n", ic->d2bxdz);

  // Example: particle access
  //printf("total particles %d of max %d\n", ions.np, ions.npmax);
  //printf("p(0)->x is %f\n", ions.p->x);
  //printf("p(0)->y is %f\n", ions.p->y);
  //printf("p(0)->z is %f\n", ions.p->z);
  //printf("p(0)->ux is %f\n", ions.p->ux);
  //printf("p(0)->uy is %f\n", ions.p->uy);
  //printf("p(0)->uz is %f\n", ions.p->uz);
  //printf("p(0)->ind is %d\n", ions.p->ind);

  //printf("np is %d\n", ions.np);
  //printf("p(np-1)->x is %f\n", (&ions.p[ions.np-1])->x);
  //printf("p(np-1)->y is %f\n", (&ions.p[ions.np-1])->y);
  //printf("p(np-1)->z is %f\n", (&ions.p[ions.np-1])->z);
  //printf("p(np-1)->ux is %f\n", (&ions.p[ions.np-1])->ux);
  //printf("p(np-1)->uy is %f\n", (&ions.p[ions.np-1])->uy);
  //printf("p(np-1)->uz is %f\n", (&ions.p[ions.np-1])->uz);
  //printf("p(np-1)->ind is %d\n", (&ions.p[ions.np-1])->ind);

  //printf("v^2 mean is %f, expected %f\n", ions.meanq_vsq(), 3*0.404*0.404);
  //printf("vx mean is %f\n", ions.meanq_vx());
  //printf("vy mean is %f\n", ions.meanq_vy());
  //printf("vz mean is %f\n", ions.meanq_vz());

  // --------------------------------------------------------------------------
  // Evolution

  int step = 0;

  while (step < par.ilast) {

    // fields at B(t=n), E(t=n)
    // particle at r(t=n), v(t=n-1/2)

    //printf(
    //  "step %d p[0] ind %d; w = %f; x,y,z = %f,%f,%f; ux,uy,uz = %f,%f,%f\n",
    //  step, (ions.p0[0]).ind, (ions.p0[0]).w,
    //  (ions.p0[0]).x,(ions.p0[0]).y,(ions.p0[0]).z,
    //  (ions.p0[0]).ux,(ions.p0[0]).uy,(ions.p0[0]).uz
    //);

    //if (step % par.isort == 0) ions.sort();

    // TODO fill ghost cells for E/B required before interp
    // I think better to display the ghost cell filling/resetting operations
    // all at the top level b/c important for evolution loop understanding.

    // update field interpolation coefficients
    ia.update();

    // particles advance to r(t=n+1), v(t=n+1/2);
    // deposit charge, current at n(t=n+1/2), j(t=n+1/2);
    // teleport across periodic grid boundaries
    ions.move_deposit();
    // smooth deposited particle moments if needed
    //ions.smooth();

    // update ghost cells for particle currents
    //fa.unload_ghost(...);

    // fields advance to B(t=n+1), E(t=n+1)
    //fa.advance_b(ions);
    // smooth deposited fields if needed
    //fa.smooth();

    step++;

    //if (step % par->idump == 0) diagnostics();

  }

  printf("interp array dbxdx %f\n", ia.voxel(3,3,3)->dbxdx);
  printf("field array bx %f\n", fa.voxel(5,5,5)->bx);
  printf("field array jfx %f\n", fa.voxel(5,5,5)->jfx);
  printf("field array jfy %f\n", fa.voxel(5,5,5)->jfy);
  printf("field array jfz %f\n", fa.voxel(5,5,5)->jfz);
  printf("field array rhof %f\n", fa.voxel(5,5,5)->rhof);

  // Pointers malloc'ed in class constructors cannot be freed in main?
  // clang compiler says "pointer being freed was not allocated"
  // don't understand why this doesn't work --ATr,2024feb10
  //free(ions.p0);
  //free(ia.ic0);
  //free(fa.f);

  return 0;
}
