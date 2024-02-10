#include <stdlib.h> // need stdlib to get int32_t at compile time???,
                    // don't understand why --ATr,2024feb09
#include <stdio.h>
#include "fields.h"
#include "interp.h"
#include "particles.h"

// Implement Boris pusher from HybridVPIC
void ParticleArray::move_boris() {

  const float one            = 1.;
  const float one_third      = 1./3.;
  const float two_fifteenths = 2./15.;
  const float cdt_dx         = fa.dt/fa.hx;  // implicit c=1 in c*dt/dx
  const float cdt_dy         = fa.dt/fa.hy;
  const float cdt_dz         = fa.dt/fa.hz;

  float qdt_2mc = (q*fa.dt)/(2*m);  // implicit c=1 in q*dt/(2*m*c)

  for (int ip=0; ip<np; ++ip) {

    particle_t* p = &(p0[ip]);

    int ix = (int)p->x;
    int iy = (int)p->y;
    int iz = (int)p->z;
    float dx = (p->x) - ix;
    float dy = (p->y) - iy;
    float dz = (p->z) - iz;

    interp_t* ic = ia.voxel(ix,iy,iz);

    float hax  = qdt_2mc * ia.exloc(ic, dx, dy, dz);
    float hay  = qdt_2mc * ia.eyloc(ic, dx, dy, dz);
    float haz  = qdt_2mc * ia.ezloc(ic, dx, dy, dz);
    float cbx  = ia.bxloc(ic, dx, dy, dz);
    float cby  = ia.byloc(ic, dx, dy, dz);
    float cbz  = ia.bzloc(ic, dx, dy, dz);

    //if (ip == 0) {
    //  printf("got ip %d ixyz %d %d %d dxyz %f %f %f ha %f %f %f cb %f %f %f\n",
    //      ip, ix,iy,iz, dx,dy,dz, hax,hay,haz, cbx,cby,cbz
    //  );
    //}

    // ---------------------------------------------
    // Boris pusher copied from Hybrid-VPIC
    // src/species_advance/standard/hyb_advance_p.cc
    // ---------------------------------------------

    float ux   = p->ux;                       // Load momentum
    float uy   = p->uy;
    float uz   = p->uz;
    //float q    = p->w;
    ux  += hax;                               // Half advance E
    uy  += hay;
    uz  += haz;
    float v0   = qdt_2mc;///sqrtf(one + (ux*ux + (uy*uy + uz*uz)));
    /**/                                      // Boris - scalars
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
    //v0   = one;///sqrtf(one + (ux*ux+ (uy*uy + uz*uz)));
    /**/                                      // Get norm displacement
    ux  *= cdt_dx;
    uy  *= cdt_dy;
    uz  *= cdt_dz;
    //ux  *= v0;
    //uy  *= v0;
    //uz  *= v0;
    v0   = dx + ux;                           // Streak midpoint (inbnds)
    v1   = dy + uy;
    v2   = dz + uz;
    dx   = v0 + ux;                         // New position
    dy   = v1 + uy;
    dz   = v2 + uz;

    // ---------------------------------------------------------
    // ACCUMULATOR CODE TO BE PORTED OVER FOR DENSITY + CURRENTS
    // ---------------------------------------------------------

    // ---------------------------------------------------------
    // BOUNDARY CODE TO BE PORTED OVER AS WELL.
    // ---------------------------------------------------------

  }

  return;
}
