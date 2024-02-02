#include <stdio.h>  // for printf
#include <omp.h>

#include "fields.h"
#include "params.h"
#include "particles.h"

int main(int argc, char* argv[]) {

  printf("Max number of threads %d\n", omp_get_max_threads());

  // --------------------------------------------------------------------------
  // Initialization - physics
  // Eventually, let user customize using arbitrary code blocks via macros

  //Simulation* sim = Simulation(...);  // may need for checkpoints
  param_t        par;  // allocated on stack (not heap) for now...
  //Grid*         gr = Grid(...);
  FieldArray      fa = FieldArray(10, 10, 10, 1);
  ParticleArray ions = ParticleArray(1, 1, 1000);

  //sim->par  = par
  //sim->grid = gr
  //sim->fa   = fa
  //sim->ions = ions

  par.isort = 20;
  par.idump = 100;
  par.last  = 1000;

  fa.uniform_b(1, 1, 1);
  fa.uniform_e(0, 0, 0);

  ions.maxwellian(0.404, 1000, 1);

  // Set initial E/B values on grid
  // This requires deposition of ion moments
  //field_advance(gr, fa, ions);

  // --------------------------------------------------------------------------
  // Print basic information about data structures
  // and show how the user can interact with them.

  // Example: field voxel pointers in 
  //printf("voxel(0,0,0) is %p\n", fa.voxel(0,0,0));
  //printf("voxel(1,0,0) is %p\n", fa.voxel(1,0,0));
  //printf("voxel(2,0,0) is %p\n", fa.voxel(2,0,0));

  // Example: field access via voxel indexing
  //printf("bx(0,0,0) is %f\n", fa.voxel(0,0,0)->bx);
  //printf("bx(1,0,0) is %f\n", fa.voxel(1,0,0)->bx);
  //printf("bx(2,0,0) is %f\n", fa.voxel(2,0,0)->bx);
  //printf("bx(11,11,11) is %f\n", fa.voxel(11,11,11)->bx);

  // Example: field access via pointer arithmetic
  field_t* f0 = fa.voxel(0,0,0);
  printf("bx(0,0,0) is %f\n", f0->bx);
  printf("bx(1,0,0) is %f\n", (f0+1)->bx);
  printf("bx(2,0,0) is %f\n", (f0+2)->bx);

  // Example: particle access
  printf("p(0)->ux is %f\n", ions.p->ux); // todo test that this is
  printf("p(0)->uy is %f\n", ions.p->uy); // actually maxwellian
  printf("p(0)->uz is %f\n", ions.p->uz);

  // --------------------------------------------------------------------------
  // Evolution
  int step = 0;

  while (step < par.last) {

    if (step % par.isort == 0) {
      ions.sort();
    }

    //move(gr, fa, ions);

    //field_advance(gr, fa, ions);

    step++;

    //if (step % par->idump == 0) {
    //  diagnostics(sim);
    //}

  }

  return 0;
}
