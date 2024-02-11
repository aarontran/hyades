#include <stdlib.h> // need stdlib to get int32_t at compile time???,
                    // don't understand why --ATr,2024feb09
#include <stdio.h>
#include "fields.h"
#include "interp.h"
#include "particles.h"

void ParticleArray::move_deposit() {

  const float one            = 1.;
  const float two            = 2.;
  const float three          = 3.;
  const float one_third      = 1./3.;
  const float one_half       = 1./2.;
  const float two_fifteenths = 2./15.;
  const float cdt_dx         = fa.dt/fa.hx;  // implicit c=1 in c*dt/dx
  const float cdt_dy         = fa.dt/fa.hy;
  const float cdt_dz         = fa.dt/fa.hz;

  float qdt_2mc = (qsp*fa.dt)/(2*msp);  // implicit c=1 in q*dt/(2*m*c)

  // =============================================
  // Particle Boris pusher copied from Hybrid-VPIC
  // src/species_advance/standard/hyb_advance_p.cc
  // =============================================

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

    if (ip == 0) {
      printf("mover ind %d xyz %f %f %f uxyz %f %f %f . . . ixyz %d %d %d dxyz %f %f %f ha %f %f %f cb %f %f %f\n",
          p->ind, p->x,p->y,p->z, p->ux,p->uy,p->uz,
          ix,iy,iz, dx,dy,dz, hax,hay,haz, cbx,cby,cbz
      );
    }

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
  // Use TWO ghost cells instead of ONE.
  // That way, any particle crossing off grid edge can still deposit, we don't
  // need to clobber...
  //
  // For my use case, I think it's best to split the loops so the mover
  // can be easily vectorized; the deposit is hard to vectorize without
  // the accumulator arrays or careful OpenMP threading...

  for (int ip=0; ip<np; ++ip) {

    particle_t* p = &(p0[ip]);

    // Unwind to get streak midpoint
    float xmh = (p->x - one_half*(p->ux)*cdt_dx);  // xmh = x minus half step
    float ymh = (p->y - one_half*(p->uy)*cdt_dx);
    float zmh = (p->z - one_half*(p->uz)*cdt_dx);

    // Midpoint voxel indices
    int ix = (int)(xmh + one_half);  // particles use cell-centered coordinates
    int iy = (int)(ymh + one_half);
    int iz = (int)(zmh + one_half);
    // Midpoint voxel offsets on interval [-1,1]
    float v0 = two*(xmh - ix);
    float v1 = two*(ymh - iy);
    float v2 = two*(zmh - iz);

    float q = p->w * qsp;

    // Ari's quadratic spline deposit scheme
    float w0 =  q*two*( three - v0*v0 - v1*v1 - v2*v2 );
    float wx =  q*( v0 + one )*( v0 + one );
    float wy =  q*( v1 + one )*( v1 + one );
    float wz =  q*( v2 + one )*( v2 + one );
    float wmx = q*( v0 - one )*( v0 - one );
    float wmy = q*( v1 - one )*( v1 - one );
    float wmz = q*( v2 - one )*( v2 - one );

    if (ip == 0) {
      printf("depst ind %d xmh %f %f %f uxyz %f %f %f . . . ixmh %d %d %d v012 %f %f %f weights %f %f %f %f %f %f %f\n",
          p->ind, xmh,ymh,zmh, p->ux,p->uy,p->uz,
          ix,iy,iz, v0,v1,v2,
          w0, wx, wy, wz, wmx, wmy, wmz
      );
    }

    field_t* f0  = fa.voxel(ix,  iy,  iz  );
    field_t* fx  = fa.voxel(ix+1,iy,  iz  );
    field_t* fy  = fa.voxel(ix,  iy+1,iz  );
    field_t* fz  = fa.voxel(ix,  iy,  iz+1);
    field_t* fmx = fa.voxel(ix-1,iy,  iz  );
    field_t* fmy = fa.voxel(ix,  iy-1,iz  );
    field_t* fmz = fa.voxel(ix,  iy,  iz-1);

    // Ari's quadratic spline deposit scheme current
     f0->jfx +=  w0*p->ux;
     fx->jfx +=  wx*p->ux;
     fy->jfx +=  wy*p->ux;
     fz->jfx +=  wz*p->ux;
    fmx->jfx += wmx*p->ux;
    fmy->jfx += wmy*p->ux;
    fmz->jfx += wmz*p->ux;

     f0->jfy +=  w0*p->uy;
     fx->jfy +=  wx*p->uy;
     fy->jfy +=  wy*p->uy;
     fz->jfy +=  wz*p->uy;
    fmx->jfy += wmx*p->uy;
    fmy->jfy += wmy*p->uy;
    fmz->jfy += wmz*p->uy;

     f0->jfz +=  w0*p->uz;
     fx->jfz +=  wx*p->uz;
     fy->jfz +=  wy*p->uz;
     fz->jfz +=  wz*p->uz;
    fmx->jfz += wmx*p->uz;
    fmy->jfz += wmy*p->uz;
    fmz->jfz += wmz*p->uz;

    // Ari's quadratic spline deposit scheme density
     f0->rhof +=  w0;
     fx->rhof +=  wx;
     fy->rhof +=  wy;
     fz->rhof +=  wz;
    fmx->rhof += wmx;
    fmy->rhof += wmy;
    fmz->rhof += wmz;

    if (ip == 0) {
      printf("depst ind %d xmh %f %f %f uxyz %f %f %f . . . jrho %f %f %f %f\n",
          p->ind, xmh,ymh,zmh, p->ux,p->uy,p->uz,
          f0->jfx, f0->jfy, f0->jfz, f0->rhof
      );
    }

  }  // end particle deposit loop

  // =====================================
  // Teleport particle positions if needed
  // =====================================

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

  for (int ip=0; ip<np; ++ip) {
    particle_t* p = &(p0[ip]);
    if (p->x <  x0) { p->x += nx; };
    if (p->x >= x1) { p->x -= nx; };
    if (p->y <  y0) { p->y += ny; };
    if (p->y >= y1) { p->y -= ny; };
    if (p->z <  z0) { p->z += nz; };
    if (p->z >= z1) { p->z -= nz; };
  }  // end particle teleport loop

  return;
}
