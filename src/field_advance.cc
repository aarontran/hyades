//#include <assert.h>
//#include <stdlib.h>  // for malloc, qsort
//#include <stdio.h>  // for printf

#include "field.h"

// ----------------------------------------------------------------------------
// advance_eb_rk4_ctrmesh() advances both E/B fields in time using a
// 4th-order Runge Kutta (RK4) algorithm, assuming all fields and ion moments
// are stored on a cell-centered mesh.
// It uses the following fields:
//    E at t_n
//    B at t_n
//    ion density at t_{n+1/2}
//    ion velocity at t_{n+1/2}
//    old ion density at t_{n-1/2}
//    old ion velocity at t_{n-1/2}
// The ion moments are extrapolated forward in time to t_{n+1}.
//
// Inputs:
//    isub = current subcycle number (from 0 to nsub-1)
//    nsub = total subcycle number
// Result:
//    E and B fields advanced 1 substep in time on live AND ghost cells.
// ----------------------------------------------------------------------------
void FieldArray::advance_eb_rk4_ctrmesh(int isub, int nsub) {

  const float r2hx = (nx > 1) ? 1/(2*hx) : 0;
  const float r2hy = (ny > 1) ? 1/(2*hy) : 0;
  const float r2hz = (nz > 1) ? 1/(2*hz) : 0;

  const float isub_ = (float) isub;
  const float nsub_ = (float) nsub;

  const float subdt       = dt/nsub_;
  const float subdt_half  = dt/nsub_ / 2.;
  const float subdt_sixth = dt/nsub_ / 6.;

  // RK4 scheme.
  // 1. Compute k1 using B_n, E_n
  // 2. Compute k2 using B_{n+1/2)',  E_{n+1/2}'
  // 3. Compute k3 using B_{n+1/2)'', E_{n+1/2}''
  // 4. Compute k4 using B_{n+1}',    E_{n+1}'
  // 5. Compute final    B_{n+1},     E_{n+1}
  //
  // Need tmp vars (bx0,by0,bz0)    to hold B_n
  //               (tmpx,tmpy,tmpz) to accumulate k1,k2,k3,k4 for final B_{n+1}
  //               (smex,smey,smez) for hyper-resistivity

  // -------
  // Setup.
  // -------
  field_t* fv = f0;
  for (int ii = 0; ii < nvall; ++ii) {
    fv->bx0 = fv->bx;
    fv->by0 = fv->by;
    fv->bz0 = fv->bz;
    ++fv;
  }

  // -------
  // Step 1.
  // -------

  // E_n
  advance_e_ctrmesh( isub_/nsub_ );

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
    // -curl(E_n)
    float k1x = - r2hy*(fy->ez - fmy->ez) + r2hz*(fz->ey - fmz->ey);  // -(∂y Ez - ∂z Ey)
    float k1y = - r2hz*(fz->ex - fmz->ex) + r2hx*(fx->ez - fmx->ez);  // -(∂z Ex - ∂x Ez)
    float k1z = - r2hx*(fx->ey - fmx->ey) + r2hy*(fy->ex - fmy->ex);  // -(∂x Ey - ∂y Ex)
    // Accumulate final B_{n+1}
    // notice "=" assignment unlike "+=" in following RK4 steps
    fv->tmpx = k1x;
    fv->tmpy = k1y;
    fv->tmpz = k1z;
    // B_{n+1/2}'
    fv->bx = fv->bx0 + subdt_half * k1x;
    fv->by = fv->by0 + subdt_half * k1y;
    fv->bz = fv->bz0 + subdt_half * k1z;
  }}}

  ghost_copy_b();

  // -------
  // Step 2.
  // -------

  // E_{n+1/2}'
  advance_e_ctrmesh( (isub_+0.5)/nsub_ );

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
    // -curl(E_{n+1/2}')
    float k2x = - r2hy*(fy->ez - fmy->ez) + r2hz*(fz->ey - fmz->ey);  // -(∂y Ez - ∂z Ey)
    float k2y = - r2hz*(fz->ex - fmz->ex) + r2hx*(fx->ez - fmx->ez);  // -(∂z Ex - ∂x Ez)
    float k2z = - r2hx*(fx->ey - fmx->ey) + r2hy*(fy->ex - fmy->ex);  // -(∂x Ey - ∂y Ex)
    // Accumulate final B_{n+1}
    fv->tmpx += 2 * k2x;
    fv->tmpy += 2 * k2y;
    fv->tmpz += 2 * k2z;
    // B_{n+1/2}''
    fv->bx = fv->bx0 + subdt_half * k2x;
    fv->by = fv->by0 + subdt_half * k2y;
    fv->bz = fv->bz0 + subdt_half * k2z;
  }}}

  ghost_copy_b();

  // -------
  // Step 3.
  // -------

  // E_{n+1/2}''
  advance_e_ctrmesh( (isub_+0.5)/nsub_ );

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
    // -curl(E_{n+1/2}'')
    float k3x = - r2hy*(fy->ez - fmy->ez) + r2hz*(fz->ey - fmz->ey);  // -(∂y Ez - ∂z Ey)
    float k3y = - r2hz*(fz->ex - fmz->ex) + r2hx*(fx->ez - fmx->ez);  // -(∂z Ex - ∂x Ez)
    float k3z = - r2hx*(fx->ey - fmx->ey) + r2hy*(fy->ex - fmy->ex);  // -(∂x Ey - ∂y Ex)
    // Accumulate final B_{n+1}
    fv->tmpx += 2 * k3x;
    fv->tmpy += 2 * k3y;
    fv->tmpz += 2 * k3z;
    // B_{n+1}'
    fv->bx = fv->bx0 + subdt * k3x;
    fv->by = fv->by0 + subdt * k3y;
    fv->bz = fv->bz0 + subdt * k3z;
  }}}

  ghost_copy_b();

  // -------
  // Step 4.
  // -------

  // E_{n+1}'
  advance_e_ctrmesh( (isub_+1.0)/nsub_ );

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
    // -curl(E_{n+1}')
    float k4x = - r2hy*(fy->ez - fmy->ez) + r2hz*(fz->ey - fmz->ey);  // -(∂y Ez - ∂z Ey)
    float k4y = - r2hz*(fz->ex - fmz->ex) + r2hx*(fx->ez - fmx->ez);  // -(∂z Ex - ∂x Ez)
    float k4z = - r2hx*(fx->ey - fmx->ey) + r2hy*(fy->ex - fmy->ex);  // -(∂x Ey - ∂y Ex)
    // Final B_{n+1}
    fv->bx = fv->bx0 + subdt_sixth * (fv->tmpx + k4x);
    fv->by = fv->by0 + subdt_sixth * (fv->tmpy + k4y);
    fv->bz = fv->bz0 + subdt_sixth * (fv->tmpz + k4z);
  }}}

  ghost_copy_b();

  // -------
  // Step 5.
  // -------

  // Final E_{n+1}
  advance_e_ctrmesh( (isub_+1.0)/nsub_ );

  return;
}

// ----------------------------------------------------------------------------
// advance_e_ctrmesh() implements Ohm's law assuming that (E, B, rhof, jf) are
// all stored on a cell-centered mesh.
// Inputs:
//    frac = 0 to 1 chooses a time between t_n to t_{n+1} by extrapolating
//           ion moments at t_{n-1/2} and t_{n+1/2}.
//    The stored B field must be updated to the "target" time already.
// Result:
//    E field updated on live AND ghost cells.
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
    // Field advance methods are responsible for updating their "own" ghosts.
    // TODO branching returns here is a bad idea (basically a GOTO),
    // split hyper-resistivity into its own method -ATr,2024feb25
    ghost_copy_e();
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

  // Field advance methods are responsible for updating their "own" ghosts.
  // TODO branching returns here is a bad idea (basically a GOTO),
  // split hyper-resistivity into its own method -ATr,2024feb25
  ghost_copy_e();
  return;

}
