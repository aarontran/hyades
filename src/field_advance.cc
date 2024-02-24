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
  // Need 3 temporary variables to hold B_n
  // Need 3 temporary variables to accumulate k1,k2,k3,k4 for final B_{n+1}

  // Setup
  field_t fv;
  for (int ii = 0; ii < nvall; ++ii) {
    fv = f0[ii];
    fv.bx0  = fv.bx;
    fv.by0  = fv.by;
    fv.bz0  = fv.bz;
  }

//  // Step 1.
//  advance_e_ctrmesh(0);       // E_n using B_n
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

void FieldArray::advance_e_ctrmesh(float frac) {
}
