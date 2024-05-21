#include <stdlib.h> // need stdlib or cstdint to get int32_t at compile time???,
                    // don't understand why --ATr,2024feb09
#include <stdio.h>

#include <omp.h>

#include "field.h"
#include "interp.h"
#include "particle.h"

// =============================================
// Particle Boris pusher copied from Hybrid-VPIC
// src/species_advance/standard/hyb_advance_p.cc
// =============================================

void ParticleArray::move() {

  const float one            = 1.;
  const float two            = 2.;
  const float one_third      = 1./3.;
  const float one_half       = 1./2.;
  const float two_fifteenths = 2./15.;
  const float qdt_2mc = (qsp*fa.dt)/(2*msp);  // q*dt/(2*m*c) with c=1
  const float cdt_dx  = fa.dt/fa.hx;          // c*dt/dx with c=1
  const float cdt_dy  = fa.dt/fa.hy;
  const float cdt_dz  = fa.dt/fa.hz;

  // TODO HACKY HACK HACK!!!!!
#ifdef SHAPE_CIC
  const float fnx = fa.nx;
  const float fny = fa.ny;
  const float fnz = fa.nz;
  const float fng = fa.ng;
  // Linear array increments
  const int ox = 1;
  const int oy = (fnx+2*fng);
  const int oz = (fnx+2*fng) * (fny+2*fng);
#endif
  // TODO HACKY HACK HACK!!!!!
#ifdef SHAPE_TSC
  const float one_eighth = 1./8.;
  const float one_fourth = 1./4.;
  const float three = 3.;
#endif

  // Note: some arithmetic ops are NOT optimized; beware compiler re-ordering
  // of math b/c you may lose precision.

#pragma omp parallel for
  for (int ip=0; ip<np; ++ip) {

    particle_t* p = &(p0[ip]);

#if defined SHAPE_NGP || SHAPE_QS
    // Voxel indices
    int ix = (int)(p->x + one_half);  // particles use cell-centered coordinates
    int iy = (int)(p->y + one_half);
    int iz = (int)(p->z + one_half);
    // Voxel offsets on interval [-1,1]
    float dx = two*((p->x) - ix);
    float dy = two*((p->y) - iy);
    float dz = two*((p->z) - iz);

    interp_t* ic = ia.voxel(ix,iy,iz);

    float hax  = qdt_2mc * ia.exloc(ic, dx, dy, dz);
    float hay  = qdt_2mc * ia.eyloc(ic, dx, dy, dz);
    float haz  = qdt_2mc * ia.ezloc(ic, dx, dy, dz);
    float cbx  = ia.bxloc(ic, dx, dy, dz);
    float cby  = ia.byloc(ic, dx, dy, dz);
    float cbz  = ia.bzloc(ic, dx, dy, dz);
#else
#ifdef SHAPE_CIC
    // This is a very hacky implementation, we can do use interpolator too but
    // it would require big data structure change.

    // Voxel indices nearest (below/left) of streak midpoint
    int ix = (int)(p->x);  // particles use cell-centered coordinates
    int iy = (int)(p->y);
    int iz = (int)(p->z);
    // "Lower/left" voxel offsets on interval [0,1]
    // unlike other parts of VPIC code
    float dx = p->x - ix;
    float dy = p->y - iy;
    float dz = p->z - iz;

    ///////////////////
    // First attempt

    // CIC = area weighting = trilinear interpolation
//    field_t* f0   = fa.voxel( ix  , iy  , iz   );
//    field_t* fz   = fa.voxel( ix  , iy  , iz+1 );
//    field_t* fy   = fa.voxel( ix  , iy+1, iz   );
//    field_t* fyz  = fa.voxel( ix  , iy+1, iz+1 );
//    field_t* fx   = fa.voxel( ix+1, iy  , iz   );
//    field_t* fxz  = fa.voxel( ix+1, iy  , iz+1 );
//    field_t* fxy  = fa.voxel( ix+1, iy+1, iz   );
//    field_t* fxyz = fa.voxel( ix+1, iy+1, iz+1 );
//
//    float hax = qdt_2mc * (
//        (1.-dx)*(1.-dy)*(1.-dz) *   f0->ex
//      + (1.-dx)*(1.-dy)*    dz  *   fz->ex
//      + (1.-dx)*    dy *(1.-dz) *   fy->ex
//      + (1.-dx)*    dy *    dz  *  fyz->ex
//      +     dx *(1.-dy)*(1.-dz) *   fx->ex
//      +     dx *(1.-dy)*    dz  *  fxz->ex
//      +     dx *    dy *(1.-dz) *  fxy->ex
//      +     dx *    dy *    dz  * fxyz->ex
//    );
//    float hay = qdt_2mc * (
//        (1.-dx)*(1.-dy)*(1.-dz) *   f0->ey
//      + (1.-dx)*(1.-dy)*    dz  *   fz->ey
//      + (1.-dx)*    dy *(1.-dz) *   fy->ey
//      + (1.-dx)*    dy *    dz  *  fyz->ey
//      +     dx *(1.-dy)*(1.-dz) *   fx->ey
//      +     dx *(1.-dy)*    dz  *  fxz->ey
//      +     dx *    dy *(1.-dz) *  fxy->ey
//      +     dx *    dy *    dz  * fxyz->ey
//    );
//    float haz = qdt_2mc * (
//        (1.-dx)*(1.-dy)*(1.-dz) *   f0->ez
//      + (1.-dx)*(1.-dy)*    dz  *   fz->ez
//      + (1.-dx)*    dy *(1.-dz) *   fy->ez
//      + (1.-dx)*    dy *    dz  *  fyz->ez
//      +     dx *(1.-dy)*(1.-dz) *   fx->ez
//      +     dx *(1.-dy)*    dz  *  fxz->ez
//      +     dx *    dy *(1.-dz) *  fxy->ez
//      +     dx *    dy *    dz  * fxyz->ez
//    );
//
//    float cbx = (
//        (1.-dx)*(1.-dy)*(1.-dz) *   f0->bx
//      + (1.-dx)*(1.-dy)*    dz  *   fz->bx
//      + (1.-dx)*    dy *(1.-dz) *   fy->bx
//      + (1.-dx)*    dy *    dz  *  fyz->bx
//      +     dx *(1.-dy)*(1.-dz) *   fx->bx
//      +     dx *(1.-dy)*    dz  *  fxz->bx
//      +     dx *    dy *(1.-dz) *  fxy->bx
//      +     dx *    dy *    dz  * fxyz->bx
//    );
//    float cby = (
//        (1.-dx)*(1.-dy)*(1.-dz) *   f0->by
//      + (1.-dx)*(1.-dy)*    dz  *   fz->by
//      + (1.-dx)*    dy *(1.-dz) *   fy->by
//      + (1.-dx)*    dy *    dz  *  fyz->by
//      +     dx *(1.-dy)*(1.-dz) *   fx->by
//      +     dx *(1.-dy)*    dz  *  fxz->by
//      +     dx *    dy *(1.-dz) *  fxy->by
//      +     dx *    dy *    dz  * fxyz->by
//    );
//    float cbz = (
//        (1.-dx)*(1.-dy)*(1.-dz) *   f0->bz
//      + (1.-dx)*(1.-dy)*    dz  *   fz->bz
//      + (1.-dx)*    dy *(1.-dz) *   fy->bz
//      + (1.-dx)*    dy *    dz  *  fyz->bz
//      +     dx *(1.-dy)*(1.-dz) *   fx->bz
//      +     dx *(1.-dy)*    dz  *  fxz->bz
//      +     dx *    dy *(1.-dz) *  fxy->bz
//      +     dx *    dy *    dz  * fxyz->bz
//    );

    ///////////////////
    // Second attempt with more direct array access
    // and refactored index calculation
    // TENTATIVE RESULT: doesn't seem to matter more than few percent level...
    int il = ix + (fnx+2*fng)*(iy + (fny+2*fng)*iz);

//    field_t* f0   = fa.voxel( ix  , iy  , iz   );
//    field_t* fz   = fa.voxel( ix  , iy  , iz+1 );
//    field_t* fy   = fa.voxel( ix  , iy+1, iz   );
//    field_t* fyz  = fa.voxel( ix  , iy+1, iz+1 );
//    field_t* fx   = fa.voxel( ix+1, iy  , iz   );
//    field_t* fxz  = fa.voxel( ix+1, iy  , iz+1 );
//    field_t* fxy  = fa.voxel( ix+1, iy+1, iz   );
//    field_t* fxyz = fa.voxel( ix+1, iy+1, iz+1 );

    float hax = qdt_2mc * (
        (1.-dx)*(1.-dy)*(1.-dz) * (fa.f0[ il          ]).ex
      + (1.-dx)*(1.-dy)*    dz  * (fa.f0[ il      +oz ]).ex
      + (1.-dx)*    dy *(1.-dz) * (fa.f0[ il   +oy    ]).ex
      + (1.-dx)*    dy *    dz  * (fa.f0[ il   +oy+oz ]).ex
      +     dx *(1.-dy)*(1.-dz) * (fa.f0[ il+ox       ]).ex
      +     dx *(1.-dy)*    dz  * (fa.f0[ il+ox   +oz ]).ex
      +     dx *    dy *(1.-dz) * (fa.f0[ il+ox+oy    ]).ex
      +     dx *    dy *    dz  * (fa.f0[ il+ox+oy+oz ]).ex
    );
    float hay = qdt_2mc * (
        (1.-dx)*(1.-dy)*(1.-dz) * (fa.f0[ il          ]).ey
      + (1.-dx)*(1.-dy)*    dz  * (fa.f0[ il      +oz ]).ey
      + (1.-dx)*    dy *(1.-dz) * (fa.f0[ il   +oy    ]).ey
      + (1.-dx)*    dy *    dz  * (fa.f0[ il   +oy+oz ]).ey
      +     dx *(1.-dy)*(1.-dz) * (fa.f0[ il+ox       ]).ey
      +     dx *(1.-dy)*    dz  * (fa.f0[ il+ox   +oz ]).ey
      +     dx *    dy *(1.-dz) * (fa.f0[ il+ox+oy    ]).ey
      +     dx *    dy *    dz  * (fa.f0[ il+ox+oy+oz ]).ey
    );
    float haz = qdt_2mc * (
        (1.-dx)*(1.-dy)*(1.-dz) * (fa.f0[ il          ]).ez
      + (1.-dx)*(1.-dy)*    dz  * (fa.f0[ il      +oz ]).ez
      + (1.-dx)*    dy *(1.-dz) * (fa.f0[ il   +oy    ]).ez
      + (1.-dx)*    dy *    dz  * (fa.f0[ il   +oy+oz ]).ez
      +     dx *(1.-dy)*(1.-dz) * (fa.f0[ il+ox       ]).ez
      +     dx *(1.-dy)*    dz  * (fa.f0[ il+ox   +oz ]).ez
      +     dx *    dy *(1.-dz) * (fa.f0[ il+ox+oy    ]).ez
      +     dx *    dy *    dz  * (fa.f0[ il+ox+oy+oz ]).ez
    );

    float cbx = (
        (1.-dx)*(1.-dy)*(1.-dz) * (fa.f0[ il          ]).bx
      + (1.-dx)*(1.-dy)*    dz  * (fa.f0[ il      +oz ]).bx
      + (1.-dx)*    dy *(1.-dz) * (fa.f0[ il   +oy    ]).bx
      + (1.-dx)*    dy *    dz  * (fa.f0[ il   +oy+oz ]).bx
      +     dx *(1.-dy)*(1.-dz) * (fa.f0[ il+ox       ]).bx
      +     dx *(1.-dy)*    dz  * (fa.f0[ il+ox   +oz ]).bx
      +     dx *    dy *(1.-dz) * (fa.f0[ il+ox+oy    ]).bx
      +     dx *    dy *    dz  * (fa.f0[ il+ox+oy+oz ]).bx
    );
    float cby = (
        (1.-dx)*(1.-dy)*(1.-dz) * (fa.f0[ il          ]).by
      + (1.-dx)*(1.-dy)*    dz  * (fa.f0[ il      +oz ]).by
      + (1.-dx)*    dy *(1.-dz) * (fa.f0[ il   +oy    ]).by
      + (1.-dx)*    dy *    dz  * (fa.f0[ il   +oy+oz ]).by
      +     dx *(1.-dy)*(1.-dz) * (fa.f0[ il+ox       ]).by
      +     dx *(1.-dy)*    dz  * (fa.f0[ il+ox   +oz ]).by
      +     dx *    dy *(1.-dz) * (fa.f0[ il+ox+oy    ]).by
      +     dx *    dy *    dz  * (fa.f0[ il+ox+oy+oz ]).by
    );
    float cbz = (
        (1.-dx)*(1.-dy)*(1.-dz) * (fa.f0[ il          ]).bz
      + (1.-dx)*(1.-dy)*    dz  * (fa.f0[ il      +oz ]).bz
      + (1.-dx)*    dy *(1.-dz) * (fa.f0[ il   +oy    ]).bz
      + (1.-dx)*    dy *    dz  * (fa.f0[ il   +oy+oz ]).bz
      +     dx *(1.-dy)*(1.-dz) * (fa.f0[ il+ox       ]).bz
      +     dx *(1.-dy)*    dz  * (fa.f0[ il+ox   +oz ]).bz
      +     dx *    dy *(1.-dz) * (fa.f0[ il+ox+oy    ]).bz
      +     dx *    dy *    dz  * (fa.f0[ il+ox+oy+oz ]).bz
    );

    ///////////////////
    // The mover logic below uses nearest voxel index and
    // offsets so adjust (ix,iy,iz) and (dx,dy,dz) back to "usual"

    // Voxel indices
    ix = (int)(p->x + one_half);  // particles use cell-centered coordinates
    iy = (int)(p->y + one_half);
    iz = (int)(p->z + one_half);
    // Voxel offsets on interval [-1,1]
    dx = two*((p->x) - ix);
    dy = two*((p->y) - iy);
    dz = two*((p->z) - iz);
#else
#ifdef SHAPE_TSC

    // Voxel indices
    int ix = (int)(p->x + one_half);  // particles use cell-centered coordinates
    int iy = (int)(p->y + one_half);
    int iz = (int)(p->z + one_half);
    // Voxel offsets on interval [-1,1]
    float dx = two*((p->x) - ix);
    float dy = two*((p->y) - iy);
    float dz = two*((p->z) - iz);

    //interp_t* ic = ia.voxel(ix,iy,iz);

    // TSC = 27 points required
    // Hockney/Eastwood Eqn (5-88)
    //float qw  = p->w * qsp_rV;
    //float wmx = 1/2 * (3/2 - (2+dx)/2) * (3/2 - (2+dx)/2)
    //float w0  = 3/4 - (dx/2)*(dx/2);
    //float wx  = 1/2 * (3/2 - (2-dx)/2) * (3/2 - (2-dx)/2)
    // Refactored form
    float wxs[3] = { one_eighth * (one-dx)*(one-dx),
                     one_fourth * (three - dx*dx),
                     one_eighth * (one+dx)*(one+dx) };
    float wys[3] = { one_eighth * (one-dy)*(one-dy),
                     one_fourth * (three - dy*dy),
                     one_eighth * (one+dy)*(one+dy) };
    float wzs[3] = { one_eighth * (one-dz)*(one-dz),
                     one_fourth * (three - dz*dz),
                     one_eighth * (one+dz)*(one+dz) };
    float hax = 0.;
    float hay = 0.;
    float haz = 0.;
    float cbx = 0.;
    float cby = 0.;
    float cbz = 0.;
    for (int ii=-1; ii<=1; ++ii) {
    for (int jj=-1; jj<=1; ++jj) {
    for (int kk=-1; kk<=1; ++kk) {
      field_t* fv = fa.voxel( ix+ii, iy+jj, iz+kk );
      float w3 = wxs[ii+1] * wys[jj+1] * wzs[kk+1];
      hax += w3 * fv->ex;
      hay += w3 * fv->ey;
      haz += w3 * fv->ez;
      cbx += w3 * fv->bx;
      cby += w3 * fv->by;
      cbz += w3 * fv->bz;
    }}}

    hax *= qdt_2mc;
    hay *= qdt_2mc;
    haz *= qdt_2mc;

// endifdef shape_TSC
#endif
// endifdef shape_CIC
#endif
// endifdef shape_NGP || SHAPE_QS
#endif

    //if (ip == 0) {
    //  printf("mover ind %d xyz %.3f %.3f %.3f uxyz % .3f % .3f % .3f . . . ixyz %d %d %d dxyz % .3f % .3f % .3f haxyz % f % f % f cbxyz % f % f % f\n",
    //      p->ind, p->x,p->y,p->z, p->ux,p->uy,p->uz,
    //      ix,iy,iz, dx,dy,dz, hax,hay,haz, cbx,cby,cbz
    //  );
    //}

    // Update momentum t-1/2 to t+1/2
    float ux   = p->ux;                       // Load momentum
    float uy   = p->uy;
    float uz   = p->uz;
    ux  += hax;                               // Half advance E
    uy  += hay;
    uz  += haz;
    float v0   = qdt_2mc;
    float v1   = cbx*cbx + (cby*cby + cbz*cbz);
    float v2   = (v0*v0)*v1;
    float v3   = v0*(one+v2*(one_third+v2*two_fifteenths));
    float v4   = v3/(one+v1*(v3*v3));
    v4  += v4;
    v0   = ux + v3*( uy*cbz - uz*cby );       // Boris - uprime
    v1   = uy + v3*( uz*cbx - ux*cbz );
    v2   = uz + v3*( ux*cby - uy*cbx );
    ux  += v4*( v1*cbz - v2*cby );            // Boris - rotation
    uy  += v4*( v2*cbx - v0*cbz );
    uz  += v4*( v0*cby - v1*cbx );
    ux  += hax;                               // Half advance E
    uy  += hay;
    uz  += haz;
    p->ux = ux;                               // Store momentum
    p->uy = uy;
    p->uz = uz;

    // Update position t to t+1
    ux  *= cdt_dx;                            // Get norm displacement
    uy  *= cdt_dy;
    uz  *= cdt_dz;
    //v0   = dx + ux;                         // Streak midpoint (inbnds)
    //v1   = dy + uy;
    //v2   = dz + uz;
    //dx   = v0 + ux;                         // New position
    //dy   = v1 + uy;
    //dz   = v2 + uz;
    dx   = dx + two*ux;                       // New position
    dy   = dy + two*uy;
    dz   = dz + two*uz;
    p->x = ix + one_half*dx;                  // Store new position
    p->y = iy + one_half*dy;
    p->z = iz + one_half*dz;

  }  // end particle mover loop

} // end ParticleArray::move()


// Similar to ParticleArray::move(), but use fields at t=0 to push particle
// momenta BACKWARDS 1/2 step in time and do not update particle positions.
// In the Boris push,
// * skip the first "half E accel"
// * apply B rotation using a half timestep
// Used only to initialize the simulation.
// This corresponds to Ari's "hyb_uncenter_p(...)"
void ParticleArray::move_uncenter() {

  const float one            = 1.;
  const float two            = 2.;
  const float one_third      = 1./3.;
  const float one_half       = 1./2.;
  const float two_fifteenths = 2./15.;
  const float qdt_2mc = -1 * (qsp*fa.dt)/(2*msp); // -1*q*dt/(2*m*c) with c=1
  const float qdt_4mc = -1 * (qsp*fa.dt)/(4*msp); // -1*q*dt/(4*m*c) with c=1

  // Note: some arithmetic ops are NOT optimized; beware compiler re-ordering
  // of math b/c you may lose precision.

  for (int ip=0; ip<np; ++ip) {

    particle_t* p = &(p0[ip]);

    // Voxel indices
    int ix = (int)(p->x + one_half);  // particles use cell-centered coordinates
    int iy = (int)(p->y + one_half);
    int iz = (int)(p->z + one_half);
    // Voxel offsets on interval [-1,1]
    float dx = two*((p->x) - ix);
    float dy = two*((p->y) - iy);
    float dz = two*((p->z) - iz);

    interp_t* ic = ia.voxel(ix,iy,iz);

    float hax  = qdt_2mc * ia.exloc(ic, dx, dy, dz);
    float hay  = qdt_2mc * ia.eyloc(ic, dx, dy, dz);
    float haz  = qdt_2mc * ia.ezloc(ic, dx, dy, dz);
    float cbx  = ia.bxloc(ic, dx, dy, dz);
    float cby  = ia.byloc(ic, dx, dy, dz);
    float cbz  = ia.bzloc(ic, dx, dy, dz);

    // Update momentum t-1/2 to t+1/2
    float ux   = p->ux;                       // Load momentum
    float uy   = p->uy;
    float uz   = p->uz;
    //ux  += hax;                             // Half advance E SKIPPED
    //uy  += hay;
    //uz  += haz;
    //float v0   = qdt_2mc;
    float v0   = qdt_4mc;                     // HALVED from usual Boris push
    float v1   = cbx*cbx + (cby*cby + cbz*cbz);
    float v2   = (v0*v0)*v1;
    float v3   = v0*(one+v2*(one_third+v2*two_fifteenths));
    float v4   = v3/(one+v1*(v3*v3));
    v4  += v4;
    v0   = ux + v3*( uy*cbz - uz*cby );       // Boris - uprime
    v1   = uy + v3*( uz*cbx - ux*cbz );
    v2   = uz + v3*( ux*cby - uy*cbx );
    ux  += v4*( v1*cbz - v2*cby );            // Boris - rotation
    uy  += v4*( v2*cbx - v0*cbz );
    uz  += v4*( v0*cby - v1*cbx );
    ux  += hax;                               // Half advance E
    uy  += hay;
    uz  += haz;
    p->ux = ux;                               // Store momentum
    p->uy = uy;
    p->uz = uz;

  }  // end particle mover loop

} // end ParticleArray::move()


// ==========================================
// Deposit charge and current density on grid
// ==========================================
//
// ------------------------------
// Notes on the deposit algorithm
// ------------------------------
//
// Hybrid-VPIC uses accumulator arrays with 28 entries
// to allow pipelining/threading of the current accumulation.
//
//   hyb_reduce_accumulators.cc does some magic for threading...
//   hyb_unload_accumulator.cc  dumps into the field arrays
//
// I am NOT GONNA DO THAT.  DO THE DUMB THING until we can prove the
// performance gain is needed.
// Deposit on field arrays instead of using separate accumulator
// HOWEVER, that means we MUST the deposit step outside the mover
// loop because the deposit cannot be easily parallelized.
//
// ------------------------------------
// Notes on cell/grid boundary handling
// ------------------------------------
//
// Numerical approaches for particles crossing cell boundaries and crossing
// grid boundaries?
//
// 1. TRISTAN: split the move/deposit loops and handle particle
//             communication during deposit step, where you need
//             to check the mid-point and final-points both.
//
// 2. Hybrid-VPIC: move/deposit together for the "common" case of
//                 particles staying entirely inside cell.  The floating
//                 point calculation is highly optimized when combined,
//                 e.g. you don't need to "unwind" the trajectory.
//                 So this works GREAT for cold plasma or small timestep.
//
// It can be useful to make exception lists for the particles that need to
// be communicated b/c that clobbers cache locality... VPIC takes this so far
// as to effectively exception lists for particles crossing CELL boundaries,
// not just GRID boundaries.
//
// Deposit before teleport, how to make it work?
// Use TWO ghost cells instead of ONE.
// That way, any particle crossing off grid edge can still deposit, we don't
// need to clobber...
//
// Teleport then deposit?  Would that work better?
//   For NGP + QS shapes, then only need one ghost cell and don't need
//   corner ghosts so a little easier to do communication.
// With other higher-order shapes, you need to deposit in corners regardless
// of whether you teleport before or after the deposit.
//
// Not obvious to me whether deposit->teleport or teleport->deposit is
// better.  Shouldn't matter if you are careful.
//
// For my use case, I think it's best to split the loops so the mover
// can be easily vectorized; the deposit is hard to vectorize without
// the accumulator arrays or careful OpenMP threading...

void ParticleArray::deposit(int unwind) {

  const float cdt_dx = fa.dt/fa.hx;  // implicit c=1 in c*dt/dx
  const float cdt_dy = fa.dt/fa.hy;
  const float cdt_dz = fa.dt/fa.hz;
  // Deposit at particle streak midpoint,
  // or deposit at current particle position?
  const float frac = (unwind == 1) ? 0.5 : 0;  // For QS scheme

  for (int ip=0; ip<np; ++ip) {

    particle_t* p = &(p0[ip]);

    // Unwind to get streak midpoint, so deposit at t=n+1/2 (not t=n+1 !!)
    float xmh = (p->x - frac*(p->ux)*cdt_dx);  // xmh = x minus half step
    float ymh = (p->y - frac*(p->uy)*cdt_dy);
    float zmh = (p->z - frac*(p->uz)*cdt_dz);

#ifdef SHAPE_QS
    deposit_one_qs( xmh, ymh, zmh, p->ux, p->uy, p->uz, p->w );
#endif
#ifdef SHAPE_NGP
    deposit_one_ngp( xmh, ymh, zmh, p->ux, p->uy, p->uz, p->w );
#endif
#ifdef SHAPE_CIC
    deposit_one_cic( xmh, ymh, zmh, p->ux, p->uy, p->uz, p->w );
#endif
#ifdef SHAPE_TSC
    deposit_one_tsc( xmh, ymh, zmh, p->ux, p->uy, p->uz, p->w );
#endif

  }  // end particle deposit loop

} // end ParticleArray::deposit()


// =====================================
// Particle deposit algorithms
// =====================================

// shape needs to be matched to interpolator
// borrowing function interface style from TRISTAN's zigzag(...)

inline void ParticleArray::deposit_one_qs(
  float xx, float yy, float zz,  // position
  float ux, float uy, float uz,  // velocity
  float wt  // charge*weight
) {

  const static float one       = 1.;
  const static float two       = 2.;
  const static float three     = 3.;
  const static float one_half  = 1./2.;
  const static float qsp_rV_12 = qsp/(12*fa.hx*fa.hy*fa.hz);

  // Voxel indices
  int ix = (int)(xx + one_half);  // particles use cell-centered coordinates
  int iy = (int)(yy + one_half);
  int iz = (int)(zz + one_half);
  // Voxel offsets on interval [-1,1]
  float v0 = two*(xx - ix);
  float v1 = two*(yy - iy);
  float v2 = two*(zz - iz);

  // Combined (charge x particle weight x inverse cell size) factor
  float qw = wt * qsp_rV_12;

  // Ari's quadratic spline deposit scheme
  // weights sum to 1, without the factor of qw
  float w0 =  qw*two*( three - v0*v0 - v1*v1 - v2*v2 );
  float wx =  qw*( v0 + one )*( v0 + one );
  float wy =  qw*( v1 + one )*( v1 + one );
  float wz =  qw*( v2 + one )*( v2 + one );
  float wmx = qw*( v0 - one )*( v0 - one );
  float wmy = qw*( v1 - one )*( v1 - one );
  float wmz = qw*( v2 - one )*( v2 - one );

  field_t* f0  = fa.voxel(ix,  iy,  iz  );
  field_t* fx  = fa.voxel(ix+1,iy,  iz  );
  field_t* fy  = fa.voxel(ix,  iy+1,iz  );
  field_t* fz  = fa.voxel(ix,  iy,  iz+1);
  field_t* fmx = fa.voxel(ix-1,iy,  iz  );
  field_t* fmy = fa.voxel(ix,  iy-1,iz  );
  field_t* fmz = fa.voxel(ix,  iy,  iz-1);

   f0->jfx +=  w0*ux;
   fx->jfx +=  wx*ux;
   fy->jfx +=  wy*ux;
   fz->jfx +=  wz*ux;
  fmx->jfx += wmx*ux;
  fmy->jfx += wmy*ux;
  fmz->jfx += wmz*ux;

   f0->jfy +=  w0*uy;
   fx->jfy +=  wx*uy;
   fy->jfy +=  wy*uy;
   fz->jfy +=  wz*uy;
  fmx->jfy += wmx*uy;
  fmy->jfy += wmy*uy;
  fmz->jfy += wmz*uy;

   f0->jfz +=  w0*uz;
   fx->jfz +=  wx*uz;
   fy->jfz +=  wy*uz;
   fz->jfz +=  wz*uz;
  fmx->jfz += wmx*uz;
  fmy->jfz += wmy*uz;
  fmz->jfz += wmz*uz;

   f0->rhof +=  w0;
   fx->rhof +=  wx;
   fy->rhof +=  wy;
   fz->rhof +=  wz;
  fmx->rhof += wmx;
  fmy->rhof += wmy;
  fmz->rhof += wmz;

} // end ParticleArray::deposit_one_qs()


inline void ParticleArray::deposit_one_ngp(
  float xx, float yy, float zz,  // position
  float ux, float uy, float uz,  // velocity
  float wt  // charge*weight
) {

  const static float one_half = 1./2.;
  const static float qsp_rV   = qsp/(fa.hx*fa.hy*fa.hz);

  // Voxel indices
  int ix = (int)(xx + one_half);  // particles use cell-centered coordinates
  int iy = (int)(yy + one_half);
  int iz = (int)(zz + one_half);

  // Combined (charge x particle weight x inverse cell size) factor
  float w0 = wt * qsp_rV;

  field_t* f0 = fa.voxel(ix,iy,iz);
  f0->jfx  += w0 * ux;
  f0->jfy  += w0 * uy;
  f0->jfz  += w0 * uz;
  f0->rhof += w0;

} // end ParticleArray::deposit_one_ngp()


inline void ParticleArray::deposit_one_cic(
  float xx, float yy, float zz,  // position
  float ux, float uy, float uz,  // velocity
  float wt  // charge*weight
) {

  const static float one      = 1.;
  const static float one_half = 1./2.;
  const static float qsp_rV   = qsp/(fa.hx*fa.hy*fa.hz);

  // Voxel indices nearest (below/left) of streak midpoint
  int ix = (int)(xx);  // particles use cell-centered coordinates
  int iy = (int)(yy);
  int iz = (int)(zz);
  // "Lower/left" voxel offsets on interval [0,1]
  // unlike other parts of VPIC code
  float dx = xx - ix;
  float dy = yy - iy;
  float dz = zz - iz;

  // CIC = area weighting = trilinear interpolation
  field_t* f0   = fa.voxel( ix  , iy  , iz   );
  field_t* fz   = fa.voxel( ix  , iy  , iz+1 );
  field_t* fy   = fa.voxel( ix  , iy+1, iz   );
  field_t* fyz  = fa.voxel( ix  , iy+1, iz+1 );
  field_t* fx   = fa.voxel( ix+1, iy  , iz   );
  field_t* fxz  = fa.voxel( ix+1, iy  , iz+1 );
  field_t* fxy  = fa.voxel( ix+1, iy+1, iz   );
  field_t* fxyz = fa.voxel( ix+1, iy+1, iz+1 );
  // No interpolation is actually required...
  // what I could do with interp array is have each voxel store its
  // 7 neighbor values?

  // Combined charge x particle weight x inverse cell size
  float qw   = wt * qsp_rV;
  float w0   = qw * (1.-dx)*(1.-dy)*(1.-dz);
  float wz   = qw * (1.-dx)*(1.-dy)*    dz ;
  float wy   = qw * (1.-dx)*    dy *(1.-dz);
  float wyz  = qw * (1.-dx)*    dy *    dz ;
  float wx   = qw *     dx *(1.-dy)*(1.-dz);
  float wxz  = qw *     dx *(1.-dy)*    dz ;
  float wxy  = qw *     dx *    dy *(1.-dz);
  float wxyz = qw *     dx *    dy *    dz ;

    f0->jfx +=   w0 * ux;
    fz->jfx +=   wz * ux;
    fy->jfx +=   wy * ux;
   fyz->jfx +=  wyz * ux;
    fx->jfx +=   wx * ux;
   fxz->jfx +=  wxz * ux;
   fxy->jfx +=  wxy * ux;
  fxyz->jfx += wxyz * ux;

    f0->jfy +=   w0 * uy;
    fz->jfy +=   wz * uy;
    fy->jfy +=   wy * uy;
   fyz->jfy +=  wyz * uy;
    fx->jfy +=   wx * uy;
   fxz->jfy +=  wxz * uy;
   fxy->jfy +=  wxy * uy;
  fxyz->jfy += wxyz * uy;

    f0->jfz +=   w0 * uz;
    fz->jfz +=   wz * uz;
    fy->jfz +=   wy * uz;
   fyz->jfz +=  wyz * uz;
    fx->jfz +=   wx * uz;
   fxz->jfz +=  wxz * uz;
   fxy->jfz +=  wxy * uz;
  fxyz->jfz += wxyz * uz;

    f0->rhof +=   w0;
    fz->rhof +=   wz;
    fy->rhof +=   wy;
   fyz->rhof +=  wyz;
    fx->rhof +=   wx;
   fxz->rhof +=  wxz;
   fxy->rhof +=  wxy;
  fxyz->rhof += wxyz;

} // end ParticleArray::deposit_one_cic()


inline void ParticleArray::deposit_one_tsc(
  float xx, float yy, float zz,  // position
  float ux, float uy, float uz,  // velocity
  float wt  // charge*weight
) {

  const static float one      = 1.;
  const static float one_half = 1./2.;
  const static float two      = 2.;
  const static float qsp_rV   = qsp/(fa.hx*fa.hy*fa.hz);

  // Midpoint voxel offsets on interval [-1,1]

  // Voxel indices
  int ix = (int)(xx + one_half);  // particles use cell-centered coordinates
  int iy = (int)(yy + one_half);
  int iz = (int)(zz + one_half);
  // Voxel offsets on interval [-1,1]
  float dx = two*(xx - ix);
  float dy = two*(yy - iy);
  float dz = two*(zz - iz);

  // TSC = 27 points required
  // Hockney/Eastwood Eqn (5-88)
  //float qw  = p->w * qsp_rV;
  //float wmx = 1/2 * (3/2 - (2+dx)/2) * (3/2 - (2+dx)/2)
  //float w0  = 3/4 - (dx/2)*(dx/2);
  //float wx  = 1/2 * (3/2 - (2-dx)/2) * (3/2 - (2-dx)/2)
  // Refactored form
  float qw     = wt * qsp_rV / 512.;
  float wxs[3] = {     (1-dx)*(1-dx),
                   2 * (3 - dx*dx),
                       (1+dx)*(1+dx) };
  float wys[3] = {     (1-dy)*(1-dy),
                   2 * (3 - dy*dy),
                       (1+dy)*(1+dy) };
  float wzs[3] = {     (1-dz)*(1-dz),
                   2 * (3 - dz*dz),
                       (1+dz)*(1+dz) };

  for (int ii=-1; ii<=1; ++ii) {
  for (int jj=-1; jj<=1; ++jj) {
  for (int kk=-1; kk<=1; ++kk) {
    field_t* fv = fa.voxel( ix+ii, iy+jj, iz+kk );
    float w3 = qw * wxs[ii+1] * wys[jj+1] * wzs[kk+1];
    fv-> jfx += w3 * ux;
    fv-> jfy += w3 * uy;
    fv-> jfz += w3 * uz;
    fv->rhof += w3;
  }}}


} // end ParticleArray::deposit_one_tsc()


// =====================================
// Teleport particle positions if needed
// =====================================
void ParticleArray::boundary_teleport() {

  // Try to be consistent about usage of >= versus > etc...
  // but maybe moot given vagaries of floating point...
  const float ng = fa.ng;
  const float nx = fa.nx;
  const float ny = fa.ny;
  const float nz = fa.nz;

  // grid boundaries at cell left/right edges
  // for nx=10, ng=2, the valid domain is
  // ghost [-0.5, 1.5]
  // live  [ 1.5,11.5]
  // ghost [11.5,13.5]
  const float x0 =    ng-0.5;
  const float x1 = nx+ng-0.5;
  const float y0 =    ng-0.5;
  const float y1 = ny+ng-0.5;
  const float z0 =    ng-0.5;
  const float z1 = nz+ng-0.5;

#pragma omp parallel for
  for (int ip=0; ip<np; ++ip) {
    particle_t* p = &(p0[ip]);
    if (p->x <  x0) { p->x += nx; };
    if (p->x >= x1) { p->x -= nx; };
    if (p->y <  y0) { p->y += ny; };
    if (p->y >= y1) { p->y -= ny; };
    if (p->z <  z0) { p->z += nz; };
    if (p->z >= z1) { p->z -= nz; };
  }  // end particle teleport loop

} // end ParticleArray::boundary_teleport()
