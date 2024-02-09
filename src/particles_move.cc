#include <stdlib.h> // need stdlib to get int32_t at compile time???,
                    // don't understand why --ATr,2024feb09
#include "fields.h"
#include "interp.h"
#include "particles.h"

// Implement Boris pusher from HybridVPIC
void ParticleArray::move_boris() {

  float hax, hay, haz;
  float cbx, cby, cbz;

  float qdt_2mc = (q*fa.dt)/(2*m);  // implicit c=1

  particle_t* p = p0;
  for (int ip=0; ip<np; ++ip) {

    int ix = (int)p->x;
    int iy = (int)p->y;
    int iz = (int)p->z;
    float dx = (p->x) - ix;
    float dy = (p->y) - iy;
    float dz = (p->z) - iz;

    interp_t* ic = ia.voxel(ix,iy,iz);

    /////////////////////////////////////////
    // BEGIN OLD CODE PORTED FROM Hybrid-VPIC
    // Hybrid-VPIC field interpolation to prtl uses 2nd derivatives!
    // f = interpolator_t instance for EACH particle
   // Interpolate E 
    hax  = qdt_2mc*( ic->ex + fa.hx*( ic->dexdx + dx*ic->d2exdx )
                            + fa.hy*( ic->dexdy + dy*ic->d2exdy )
                            + fa.hz*( ic->dexdz + dz*ic->d2exdz ) );
    hay  = qdt_2mc*( ic->ey + fa.hx*( ic->deydx + dx*ic->d2eydx )
                            + fa.hy*( ic->deydy + dy*ic->d2eydy )
                            + fa.hz*( ic->deydz + dz*ic->d2eydz ) );
    haz  = qdt_2mc*( ic->ez + fa.hx*( ic->dezdx + dx*ic->d2ezdx )
                            + fa.hy*( ic->dezdy + dy*ic->d2ezdy )
                            + fa.hz*( ic->dezdz + dz*ic->d2ezdz ) );
    cbx  = ic->bx + fa.hx*( ic->dbxdx + dx*ic->d2bxdx )
                  + fa.hy*( ic->dbxdy + dy*ic->d2bxdy )
                  + fa.hz*( ic->dbxdz + dz*ic->d2bxdz );
    cby  = ic->by + fa.hx*( ic->dbydx + dx*ic->d2bydx )
                  + fa.hy*( ic->dbydy + dy*ic->d2bydy )
                  + fa.hz*( ic->dbydz + dz*ic->d2bydz );
    cbz  = ic->bz + fa.hx*( ic->dbzdx + dx*ic->d2bzdx )
                  + fa.hy*( ic->dbzdy + dy*ic->d2bzdy )
                  + fa.hz*( ic->dbzdz + dz*ic->d2bzdz );
//
//    ux   = p->ux;                             // Load momentum
//    uy   = p->uy;
//    uz   = p->uz;
//    q    = p->w;
//    ux  += hax;                               // Half advance E
//    uy  += hay;
//    uz  += haz;
//    v0   = qdt_2mc;///sqrtf(one + (ux*ux + (uy*uy + uz*uz)));
//    /**/                                      // Boris - scalars
//    v1   = cbx*cbx + (cby*cby + cbz*cbz);
//    v2   = (v0*v0)*v1;
//    v3   = v0*(one+v2*(one_third+v2*two_fifteenths));
//    v4   = v3/(one+v1*(v3*v3));
//    v4  += v4;
//    v0   = ux + v3*( uy*cbz - uz*cby );       // Boris - uprime
//    v1   = uy + v3*( uz*cbx - ux*cbz );
//    v2   = uz + v3*( ux*cby - uy*cbx );
//    ux  += v4*( v1*cbz - v2*cby );            // Boris - rotation
//    uy  += v4*( v2*cbx - v0*cbz );
//    uz  += v4*( v0*cby - v1*cbx );
//    ux  += hax;                               // Half advance E
//    uy  += hay;
//    uz  += haz;
//    p->ux = ux;                               // Store momentum
//    p->uy = uy;
//    p->uz = uz;
//    //v0   = one;///sqrtf(one + (ux*ux+ (uy*uy + uz*uz)));
//    /**/                                      // Get norm displacement
//    ux  *= cdt_dx;
//    uy  *= cdt_dy;
//    uz  *= cdt_dz;
//    //ux  *= v0;
//    //uy  *= v0;
//    //uz  *= v0;
//    v0   = dx + ux;                           // Streak midpoint (inbnds)
//    v1   = dy + uy;
//    v2   = dz + uz;
//    dx   = v0 + ux;                         // New position
//    dy   = v1 + uy;
//    dz   = v2 + uz;

    // ---------------------------------------------------------
    // ACCUMULATOR CODE TO BE PORTED OVER FOR DENSITY + CURRENTS
    // ---------------------------------------------------------

    // ---------------------------------------------------------
    // BOUNDARY CODE TO BE PORTED OVER AS WELL.
    // ---------------------------------------------------------

    // END OLD CODE PORTED FROM Hybrid-VPIC
    /////////////////////////////////////////

    ++p0;
  }

  return;
}
