#include <stdlib.h>  // for malloc
#include <stdio.h>  // for printf

#include "field.h"
#include "interp.h"

// Constructor
InterpArray::InterpArray(FieldArray fa_) :
  fa(fa_)
{
  ic0 = (interp_t*) malloc( fa.nvall*sizeof(interp_t) );
  return;
}

// Access interpolator using 3D integer tuple of grid mesh indices, 0-indexed
interp_t* InterpArray::voxel(int ii, int jj, int kk) {
  return &ic0[ fa.ivoxel(ii,jj,kk) ];
}

#ifdef SHAPE_QS
// Update the interpolation coefficients from field quantities
void InterpArray::update() {

  const float twelfth = 1./12.;
  //const float eighth = 0.125;
  const float sixth = 1./6.;
  //const float fourth = 0.25;
  //const float three_fourths = 0.75;
  const float two = 2.0;
  const float six = 6.0;

  // Loop over non-ghost cells only
  for (int kk=fa.ng; kk < (fa.nz+fa.ng); ++kk) {
    for (int jj=fa.ng; jj < (fa.ny+fa.ng); ++jj) {
      for (int ii=fa.ng; ii < (fa.nx+fa.ng); ++ii) {

        interp_t* ic = voxel(ii, jj, kk);

        // Field interpolation code is ported from
        // sf_interface/hyb_interpolator_array.cc
        field_t* pf0  = fa.voxel( ii   , jj   , kk   );
        field_t* pfx  = fa.voxel( ii+1 , jj   , kk   );
        field_t* pfy  = fa.voxel( ii   , jj+1 , kk   );
        field_t* pfz  = fa.voxel( ii   , jj   , kk+1 );
        field_t* pfmx = fa.voxel( ii-1 , jj   , kk   );
        field_t* pfmy = fa.voxel( ii   , jj-1 , kk   );
        field_t* pfmz = fa.voxel( ii   , jj   , kk-1 );

        // Copy-pasting is not elegant, but we might experiment with different
        // field interpolation schemes for E and B fields, so better to not
        // refactor anything yet.

        // Electric field interpolation, 2nd-order accurate centered stencils
        float w0   = pf0->ex;
        float wx   = pfx->ex;
        float wy   = pfy->ex;
        float wz   = pfz->ex;
        float wmx  = pfmx->ex;
        float wmy  = pfmy->ex;
        float wmz  = pfmz->ex;
        ic->ex     = twelfth*(six*w0 + wx + wy + wz + wmx + wmy + wmz);
        ic->dexdx  = sixth*(wx - wmx);
        ic->dexdy  = sixth*(wy - wmy);
        ic->dexdz  = sixth*(wz - wmz);
        ic->d2exdx = twelfth*(wx + wmx - two*w0);
        ic->d2exdy = twelfth*(wy + wmy - two*w0);
        ic->d2exdz = twelfth*(wz + wmz - two*w0);

        w0   = pf0->ey;
        wx   = pfx->ey;
        wy   = pfy->ey;
        wz   = pfz->ey;
        wmx  = pfmx->ey;
        wmy  = pfmy->ey;
        wmz  = pfmz->ey;
        ic->ey     = twelfth*(six*w0 + wx + wy + wz + wmx + wmy + wmz);
        ic->deydx  = sixth*(wx - wmx);
        ic->deydy  = sixth*(wy - wmy);
        ic->deydz  = sixth*(wz - wmz);
        ic->d2eydx = twelfth*(wx + wmx - two*w0);
        ic->d2eydy = twelfth*(wy + wmy - two*w0);
        ic->d2eydz = twelfth*(wz + wmz - two*w0);

        w0   = pf0->ez;
        wx   = pfx->ez;
        wy   = pfy->ez;
        wz   = pfz->ez;
        wmx  = pfmx->ez;
        wmy  = pfmy->ez;
        wmz  = pfmz->ez;
        ic->ez     = twelfth*(six*w0 + wx + wy + wz + wmx + wmy + wmz);
        ic->dezdx  = sixth*(wx - wmx);
        ic->dezdy  = sixth*(wy - wmy);
        ic->dezdz  = sixth*(wz - wmz);
        ic->d2ezdx = twelfth*(wx + wmx - two*w0);
        ic->d2ezdy = twelfth*(wy + wmy - two*w0);
        ic->d2ezdz = twelfth*(wz + wmz - two*w0);

        // Magnetic field interpolation, 2nd-order accurate centered stencils
        w0   = pf0->bx;
        wx   = pfx->bx;
        wy   = pfy->bx;
        wz   = pfz->bx;
        wmx  = pfmx->bx;
        wmy  = pfmy->bx;
        wmz  = pfmz->bx;
        ic->bx     = twelfth*(six*w0 + wx + wy + wz + wmx + wmy + wmz);
        ic->dbxdx  = sixth*(wx - wmx);
        ic->dbxdy  = sixth*(wy - wmy);
        ic->dbxdz  = sixth*(wz - wmz);
        ic->d2bxdx = twelfth*(wx + wmx - two*w0);
        ic->d2bxdy = twelfth*(wy + wmy - two*w0);
        ic->d2bxdz = twelfth*(wz + wmz - two*w0);

        w0   = pf0->by;
        wx   = pfx->by;
        wy   = pfy->by;
        wz   = pfz->by;
        wmx  = pfmx->by;
        wmy  = pfmy->by;
        wmz  = pfmz->by;
        ic->by     = twelfth*(six*w0 + wx + wy + wz + wmx + wmy + wmz);
        ic->dbydx  = sixth*(wx - wmx);
        ic->dbydy  = sixth*(wy - wmy);
        ic->dbydz  = sixth*(wz - wmz);
        ic->d2bydx = twelfth*(wx + wmx - two*w0);
        ic->d2bydy = twelfth*(wy + wmy - two*w0);
        ic->d2bydz = twelfth*(wz + wmz - two*w0);

        w0   = pf0->bz;
        wx   = pfx->bz;
        wy   = pfy->bz;
        wz   = pfz->bz;
        wmx  = pfmx->bz;
        wmy  = pfmy->bz;
        wmz  = pfmz->bz;
        ic->bz     = twelfth*(six*w0 + wx + wy + wz + wmx + wmy + wmz);
        ic->dbzdx  = sixth*(wx - wmx);
        ic->dbzdy  = sixth*(wy - wmy);
        ic->dbzdz  = sixth*(wz - wmz);
        ic->d2bzdx = twelfth*(wx + wmx - two*w0);
        ic->d2bzdy = twelfth*(wy + wmy - two*w0);
        ic->d2bzdz = twelfth*(wz + wmz - two*w0);

      } // end for ii
    } // end for jj
  } // end for kk
} // end InterpArray::update()
#else
#ifdef SHAPE_NGP
// Update the interpolation coefficients from field quantities
void InterpArray::update() {
  // Loop over non-ghost cells only
  for (int kk=fa.ng; kk < (fa.nz+fa.ng); ++kk) {
    for (int jj=fa.ng; jj < (fa.ny+fa.ng); ++jj) {
      for (int ii=fa.ng; ii < (fa.nx+fa.ng); ++ii) {

        interp_t* ic = voxel(ii, jj, kk);

        field_t* pf0  = fa.voxel(ii, jj, kk);

        // Electric field interpolation, NGP
        ic->ex     = pf0->ex;
        ic->ey     = pf0->ey;
        ic->ez     = pf0->ez;
        // Magnetic field interpolation, NGP
        ic->bx     = pf0->bx;
        ic->by     = pf0->by;
        ic->bz     = pf0->bz;

      } // end for ii
    } // end for jj
  } // end for kk
} // end InterpArray::update()
#else
#if defined SHAPE_CIC || defined SHAPE_TSC

// No-op, fields are accessed directly in the particle mover
void InterpArray::update() {
  return;
}

// endifdef SHAPE_CIC || defined SHAPE_TSC
#endif
// endifdef SHAPE_NGP
#endif
// endifdef SHAPE_QS
#endif

// Interpolate field(s) to (x,y,z) location given by voxel + offsets.
// Voxel offsets dx,dy,dz are within interval [-1,1].
#ifdef SHAPE_QS
float InterpArray::exloc(interp_t* ic, float dx, float dy, float dz) {
    return ( ic->ex + dx*( ic->dexdx + dx*ic->d2exdx )
                    + dy*( ic->dexdy + dy*ic->d2exdy )
                    + dz*( ic->dexdz + dz*ic->d2exdz ) );
}
float InterpArray::eyloc(interp_t* ic, float dx, float dy, float dz) {
    return ( ic->ey + dx*( ic->deydx + dx*ic->d2eydx )
                    + dy*( ic->deydy + dy*ic->d2eydy )
                    + dz*( ic->deydz + dz*ic->d2eydz ) );
}
float InterpArray::ezloc(interp_t* ic, float dx, float dy, float dz) {
    return ( ic->ez + dx*( ic->dezdx + dx*ic->d2ezdx )
                    + dy*( ic->dezdy + dy*ic->d2ezdy )
                    + dz*( ic->dezdz + dz*ic->d2ezdz ) );
}
float InterpArray::bxloc(interp_t* ic, float dx, float dy, float dz) {
    return ( ic->bx + dx*( ic->dbxdx + dx*ic->d2bxdx )
                    + dy*( ic->dbxdy + dy*ic->d2bxdy )
                    + dz*( ic->dbxdz + dz*ic->d2bxdz ) );
}
float InterpArray::byloc(interp_t* ic, float dx, float dy, float dz) {
    return ( ic->by + dx*( ic->dbydx + dx*ic->d2bydx )
                    + dy*( ic->dbydy + dy*ic->d2bydy )
                    + dz*( ic->dbydz + dz*ic->d2bydz ) );
}
float InterpArray::bzloc(interp_t* ic, float dx, float dy, float dz) {
    return ( ic->bz + dx*( ic->dbzdx + dx*ic->d2bzdx )
                    + dy*( ic->dbzdy + dy*ic->d2bzdy )
                    + dz*( ic->dbzdz + dz*ic->d2bzdz ) );
}
#else
#ifdef SHAPE_NGP
float InterpArray::exloc(interp_t* ic, float dx, float dy, float dz) {
    return ic->ex;
}
float InterpArray::eyloc(interp_t* ic, float dx, float dy, float dz) {
    return ic->ey;
}
float InterpArray::ezloc(interp_t* ic, float dx, float dy, float dz) {
    return ic->ez;
}
float InterpArray::bxloc(interp_t* ic, float dx, float dy, float dz) {
    return ic->bx;
}
float InterpArray::byloc(interp_t* ic, float dx, float dy, float dz) {
    return ic->by;
}
float InterpArray::bzloc(interp_t* ic, float dx, float dy, float dz) {
    return ic->bz;
}
#else
#if defined SHAPE_CIC || defined SHAPE_TSC

// No-op, fields are accessed directly in the particle mover
float InterpArray::exloc(interp_t* ic, float dx, float dy, float dz) { return 0; }
float InterpArray::eyloc(interp_t* ic, float dx, float dy, float dz) { return 0; }
float InterpArray::ezloc(interp_t* ic, float dx, float dy, float dz) { return 0; }
float InterpArray::bxloc(interp_t* ic, float dx, float dy, float dz) { return 0; }
float InterpArray::byloc(interp_t* ic, float dx, float dy, float dz) { return 0; }
float InterpArray::bzloc(interp_t* ic, float dx, float dy, float dz) { return 0; }

// endifdef SHAPE_CIC || defined SHAPE_TSC
#endif
// endifdef SHAPE_NGP
#endif
// endifdef SHAPE_QS
#endif
