#include <stdlib.h>  // for malloc, free
#include <stdio.h>  // for printf
#include <omp.h>

#include "fields.h"
#include "interp.h"
#include "params.h"
#include "particles.h"
#include "random.h"

int main(int argc, char* argv[]) {

  printf("Max number of threads %d\n", omp_get_max_threads());

  // --------------------------------------------------------------------------
  // Initialization - physics
  // Eventually, let user customize using arbitrary code blocks via macros

  Random rng = Random(1);  // choose random seed

  //Simulation* sim = Simulation(...);  // may need for checkpoints
  param_t        par;  // allocated on stack (not heap) for now...
  FieldArray      fa = FieldArray(10, 10, 10, 1, 1.0, 1.0, 1.0, 0.01);
  InterpArray     ia = InterpArray(fa);
  ParticleArray ions = ParticleArray(1, 1, 1000000, fa, ia, rng);

  double vth = 0.404;

  //sim->par  = par
  //sim->fa   = fa
  //sim->ions = ions

  //par.isort = 20;
  //par.idump = 100;
  par.last  = 1000;

  fa.uniform_b(1, 1, 1);
  fa.uniform_e(0, 0, 0);
  //fa.update_ghost();  // TODO populate ghost cells

  ia.update();  // Re-compute field interpolation coefficients

  ions.initialize(1000);  // index 0 to 999
  ions.maxwellian(  0,  500, vth,  600, 0, 0);
  ions.maxwellian(500, 1000, vth, -300, 0, 0);
  ions.uniform(0, 1000, -5, 5, -5, 5, -5, 5);

  // Set initial E/B values on grid
  // This requires deposition of ion moments
  //field_advance(fa, ions);

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
  //field_t* f0 = fa.voxel(0,0,0);
  //printf("bx(0,0,0) is %f\n", f0->bx);
  //printf("bx(1,0,0) is %f\n", (f0+1)->bx);
  //printf("bx(2,0,0) is %f\n", (f0+2)->bx);

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

  //printf("p(1)->x is %f\n", (&ions.p[1])->x);
  //printf("p(1)->y is %f\n", (&ions.p[1])->y);
  //printf("p(1)->z is %f\n", (&ions.p[1])->z);
  //printf("p(1)->ux is %f\n", (&ions.p[1])->ux);
  //printf("p(1)->uy is %f\n", (&ions.p[1])->uy);
  //printf("p(1)->uz is %f\n", (&ions.p[1])->uz);
  //printf("p(1)->ind is %d\n", (&ions.p[1])->ind);

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

  while (step < par.last) {

    //if (step % par.isort == 0) {
    //  ions.sort();
    //}

    ions.move_boris();

    //field_advance(fa, ions);

    step++;

    //if (step % par->idump == 0) {
    //  diagnostics(sim);
    //}

  }

  printf("ions array %p\n", ions.p0);
  printf("interp array %p\n", ia.ic0);
  printf("field array %p\n", fa.f);

  printf("ions array %f\n", ions.p0[10].x);
  printf("interp array %f\n", ia.ic0[10].dbxdx);
  printf("field array %f\n", fa.f[10].bx);

  // Pointers malloc'ed in class constructors cannot be freed in main?
  // clang compiler says "pointer being freed was not allocated"
  // don't understand why this doesn't work --ATr,2024feb10
  //free(ions.p0);
  //free(ia.ic0);
  //free(fa.f);

  return 0;
}
