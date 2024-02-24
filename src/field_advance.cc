//#include <assert.h>
//#include <stdlib.h>  // for malloc, qsort
//#include <stdio.h>  // for printf

#include "field.h"

// RK4 scheme for magnetic field advance
// assuming all fields and ion moments are defined at cell centers.
// TODO subcycling not implemented yet
void FieldArray::advance_eb_rk4_ctrmesh() {

  // Current state:
  // E at time n
  // B at time n
  // particles r at n+1
  // particles v at n+1/2

  // RK4 scheme for B-field advance.
  // 1. Compute k1 using B_n, E_n
  // 2. Compute k2 using B_{n+1/2)',  E_{n+1/2}'  <- depends on k1
  // 3. Compute k3 using B_{n+1/2)'', E_{n+1/2}'' <- depends on k2
  // 4. Compute k4 using B_{n+1}',    E_{n+1}'    <- depends on k3
  // 5. Compute final    B_{n+1},     E_{n+1}
  //
  // How to do the memory?
  // Need 3 tmp vars (bx0,by0,bz0)    to hold B_n
  // Need 3 tmp vars (tmpx,tmpy,tmpz) to accumulate k1,k2,k3,k4 for final B_{n+1}
  // Need 3 tmp vars (smex,smey,smez) for hyper-resistivity

  // Setup
  field_t fv;
  for (int ii = 0; ii < nvall; ++ii) {
    fv = f0[ii];
    fv.bx0  = fv.bx;
    fv.by0  = fv.by;
    fv.bz0  = fv.bz;
  }

  // Step 1.
  advance_e_ctrmesh(0.);       // E_n using B_n
//  k1 = -1*curl(e);
//  tmp = 1/6 * k1;
//  b = b0 + 1/2 * k1;  // B_{n+1/2}'
//
//  // Step 2.
//  advance_e_ctrmesh(0.5);     // E_{n+1/2}'
//  k2 = -1*curl(e);
//  tmp += 1/3 * k2;
//  b = b0 + 1/2 * k2;  // B_{n+1/2}''
//
//  // Step 3.
//  advance_e_ctrmesh(0.5);     // E_{n+1/2}''
//  k3 = -1*curl(e);
//  tmp += 1/3 * k3;
//  b = b0 + k3;        // B_{n+1}'
//
//  // Step 4.
//  advance_e_ctrmesh(1);       // E_{n+1}'
//  k4 = -1*curl(e);
//  tmp += 1/6 * k4;
//  b = b0 + tmp;       // B_{n+1} final
//
//  // Step 5.
//  advance_e_ctrmesh(1);       // E_{n+1} final

}

// ----------------------------------------------------------------------------
// advance_e_ctrmesh() implements Ohm's law assuming that (E, B, rhof, jf) are
// all stored on a cell-centered mesh.
// Inputs:
//    frac = 0 to 1 chooses a time between t_n to t_{n+1} for ion moments.
//    The stored B field must updated to this time already.
// Result:
//    E field updated.
// ----------------------------------------------------------------------------
void FieldArray::advance_e_ctrmesh(float frac) {

  // for finite difference stencils, "r" = reciprocal
  const float rhx = (nx > 1) ? 1/hx : 0;
  const float rhy = (ny > 1) ? 1/hy : 0;
  const float rhz = (nz > 1) ? 1/hz : 0;
  const float r2hx = (nx > 1) ? 1/(2*hx) : 0;
  const float r2hy = (ny > 1) ? 1/(2*hy) : 0;
  const float r2hz = (nz > 1) ? 1/(2*hz) : 0;
  const float rhxhx = (nx > 1) ? 1/(hx*hx) : 0;
  const float rhyhy = (ny > 1) ? 1/(hy*hy) : 0;
  const float rhzhz = (nz > 1) ? 1/(hz*hz) : 0;

  // for electron pressure
  const float hyb_te_ne_ref = hyb_te_ref_ / hyb_ne_ref_;

  // Compute E field on live cells only
  for (int kk = ng; kk < (nz+ng); ++kk) {
  for (int jj = ng; jj < (ny+ng); ++jj) {
  for (int ii = ng; ii < (nx+ng); ++ii) {

    field_t* fv  = voxel(ii  ,jj  ,kk  );
    field_t* fmx = voxel(ii-1,jj  ,kk  );
    field_t* fmy = voxel(ii  ,jj-1,kk  );
    field_t* fmz = voxel(ii  ,jj  ,kk-1);
    field_t* fx  = voxel(ii+1,jj  ,kk  );
    field_t* fy  = voxel(ii  ,jj+1,kk  );
    field_t* fz  = voxel(ii  ,jj  ,kk+1);

    // In Ari's current scheme, both rhof and jf are saved at t_{n+1/2},
    // and the "old" quantities are at t_{n-1/2}.
    // When frac=0, we get t_n     = 1/2*t_{n+1/2}) + 1/2*t_{n-1/2}
    // When frac=1, we get t_{n+1} = 3/2*t_{n+1/2}  - 1/2*t_{n-1/2}
    // Interpolate
    float rho = 0.5*((1.-frac)*(fv->rhof + fv->rhof0) + frac*(3*fv->rhof - fv->rhof0));
    float jfx = 0.5*((1.-frac)*(fv->jfx  + fv->jfx0 ) + frac*(3*fv->jfx  - fv->jfx0 ));
    float jfy = 0.5*((1.-frac)*(fv->jfy  + fv->jfy0 ) + frac*(3*fv->jfy  - fv->jfy0 ));
    float jfz = 0.5*((1.-frac)*(fv->jfz  + fv->jfz0 ) + frac*(3*fv->jfz  - fv->jfz0 ));

    // Electron pressure needs density on all face-adjacent neighbors.
    // TODO (low priority) inefficient b/c re-computed several times per voxel.
    // Try separate loop to cache time-interpolated rho for each voxel?
    // --ATr,2024feb24
    float rhox  = 0.5*((1.-frac)*( fx->rhof +  fx->rhof0) + frac*(3* fx->rhof -  fx->rhof0));
    float rhoy  = 0.5*((1.-frac)*( fy->rhof +  fy->rhof0) + frac*(3* fy->rhof -  fy->rhof0));
    float rhoz  = 0.5*((1.-frac)*( fz->rhof +  fz->rhof0) + frac*(3* fz->rhof -  fz->rhof0));
    float rhomx = 0.5*((1.-frac)*(fmx->rhof + fmx->rhof0) + frac*(3*fmx->rhof - fmx->rhof0));
    float rhomy = 0.5*((1.-frac)*(fmy->rhof + fmy->rhof0) + frac*(3*fmy->rhof - fmy->rhof0));
    float rhomz = 0.5*((1.-frac)*(fmz->rhof + fmz->rhof0) + frac*(3*fmz->rhof - fmz->rhof0));
    // Floor to avoid division by zero in vacuum regions
    rho   = ( rho  > hyb_ne_floor_) ?  rho  : hyb_ne_floor_;
    rhox  = ( rhox > hyb_ne_floor_) ?  rhox : hyb_ne_floor_;
    rhoy  = ( rhoy > hyb_ne_floor_) ?  rhoy : hyb_ne_floor_;
    rhoz  = ( rhoz > hyb_ne_floor_) ?  rhoz : hyb_ne_floor_;
    rhomx = (rhomx > hyb_ne_floor_) ? rhomx : hyb_ne_floor_;
    rhomy = (rhomy > hyb_ne_floor_) ? rhomy : hyb_ne_floor_;
    rhomz = (rhomz > hyb_ne_floor_) ? rhomz : hyb_ne_floor_;

    // Compute all finite-difference stencils together,
    // before going into Ohm's law

    // Electron pressure is hard-coded as isothermal for now,
    // for simplicity and efficiency.
    float divPx = hyb_te_ne_ref * r2hx*(rhox - rhomx);
    float divPy = hyb_te_ne_ref * r2hy*(rhoy - rhomy);
    float divPz = hyb_te_ne_ref * r2hz*(rhoz - rhomz);

    float curlbx = r2hy*(fy->bz - fmy->bz) - r2hz*(fz->by - fmz->by);  // ∂y Bz - ∂z By
    float curlby = r2hz*(fz->bx - fmz->bx) - r2hx*(fx->bz - fmx->bz);  // ∂z Bx - ∂x Bz
    float curlbz = r2hx*(fx->by - fmx->by) - r2hy*(fy->bx - fmy->bx);  // ∂x By - ∂y Bx

    // Compute E via Ohm's law
    // Try to minimize the number of operations by
    // (1) factoring arithmetic expressions,
    // (2) pulling out all pieces that require multi-point stencils
    //     (i.e., possibly non-local memory access)
    // and not duplicating any of said non-local accesses.

    float invrho = 1./rho;

    fv->ex = (
      invrho * ( -    jfy * fv->bz +    jfz * fv->by  // -(u x B)
                 + curlby * fv->bz - curlbz * fv->by  // + curl(B) x B / n_e
                 - divPx )                            // - div(P_e) / n_e
      + hyb_eta_ * curlbx                             // eta*J
    );

    fv->ey = (
      invrho * ( -    jfz * fv->bx +    jfx * fv->bz  // -(u x B)
                 + curlbz * fv->bx - curlbx * fv->bz  // + curl(B) x B / n_e
                 - divPy )                            // - div(P_e) / n_e
      + hyb_eta_ * curlby                             // eta*J
    );

    fv->ez = (
      invrho * ( -    jfx * fv->by +    jfy * fv->bx  // -(u x B)
                 + curlbx * fv->by - curlby * fv->bx  // + curl(B) x B / n_e
                 - divPz )                            // - div(P_e) / n_e
      + hyb_eta_ * curlbz                             // eta*J
    );

  }}}

  if (hyb_hypereta_ <= 0) {
    return;
  }

  // Cache del^2(B) for hyper-resistivity,
  // including one layer of ghost cells(!).
  for (int kk = ng-1; kk < (nz+ng+1); ++kk) {
  for (int jj = ng-1; jj < (ny+ng+1); ++jj) {
  for (int ii = ng-1; ii < (nx+ng+1); ++ii) {

    field_t* fv  = voxel(ii  ,jj  ,kk  );
    field_t* fmx = voxel(ii-1,jj  ,kk  );
    field_t* fmy = voxel(ii  ,jj-1,kk  );
    field_t* fmz = voxel(ii  ,jj  ,kk-1);
    field_t* fx  = voxel(ii+1,jj  ,kk  );
    field_t* fy  = voxel(ii  ,jj+1,kk  );
    field_t* fz  = voxel(ii  ,jj  ,kk+1);

    fv->smex = (  rhxhx*( fx->bx + fmx->bx - 2*fv->bx )
                + rhyhy*( fy->bx + fmy->bx - 2*fv->bx )
                + rhzhz*( fz->bx + fmz->bx - 2*fv->bx ) );
    fv->smey = (  rhxhx*( fx->by + fmx->by - 2*fv->by )
                + rhyhy*( fy->by + fmy->by - 2*fv->by )
                + rhzhz*( fz->by + fmz->by - 2*fv->by ) );
    fv->smez = (  rhxhx*( fx->bz + fmx->bz - 2*fv->bz )
                + rhyhy*( fy->bz + fmy->bz - 2*fv->bz )
                + rhzhz*( fz->bz + fmz->bz - 2*fv->bz ) );
  }}}

  // Compute hyper-resistive term in Ohm's law
  // TODO Ari implements curl(del^2(B)) but the equation is del^2(curl(B)),
  // confirm that operators commute? --ATr,2024feb24
  // TODO do curl(B) stencils need another factor of two? --ATr,2024feb24
  for (int kk = ng; kk < (nz+ng); ++kk) {
  for (int jj = ng; jj < (ny+ng); ++jj) {
  for (int ii = ng; ii < (nx+ng); ++ii) {

    field_t* fv  = voxel(ii  ,jj  ,kk  );
    field_t* fmx = voxel(ii-1,jj  ,kk  );
    field_t* fmy = voxel(ii  ,jj-1,kk  );
    field_t* fmz = voxel(ii  ,jj  ,kk-1);
    field_t* fx  = voxel(ii+1,jj  ,kk  );
    field_t* fy  = voxel(ii  ,jj+1,kk  );
    field_t* fz  = voxel(ii  ,jj  ,kk+1);

    fv->ex -= hyb_hypereta_*( rhy*(fy->smez - fmy->smez) - rhz*(fz->smey - fmz->smey) );
    fv->ey -= hyb_hypereta_*( rhz*(fz->smex - fmz->smex) - rhx*(fx->smez - fmx->smez) );
    fv->ez -= hyb_hypereta_*( rhx*(fx->smey - fmx->smey) - rhy*(fy->smex - fmy->smex) );
  }}}

}
