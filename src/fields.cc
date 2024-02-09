#include <stdlib.h>  // for malloc

#include "fields.h"

// Constructor
FieldArray::FieldArray(int nx_, int ny_, int nz_, int ng_,
                       float hx_, float hy_, float hz_, float dt_) {
  nx = nx_;
  ny = ny_;
  nz = nz_;
  ng = ng_;
  nv = (nx+2*ng) * (ny+2*ng) * (nz+2*ng);

  hx = hx_;
  hy = hy_;
  hz = hz_;
  dt = dt_;

  f = (field_t*) malloc( nv*sizeof(field_t) );
}

// Access fields using 3D integer tuple of grid mesh indices, 0-indexed
int FieldArray::ivoxel(int ii, int jj, int kk) {
  return ii + (nx+2*ng)*(jj + (ny+2*ng)*kk);
}
field_t* FieldArray::voxel(int ii, int jj, int kk) {
  return &f[ ivoxel(ii,jj,kk) ];
}

// Set B field to a uniform value
void FieldArray::uniform_b(float bx, float by, float bz) {
  field_t* f0 = f;
  // Pointer arithmetic version
  for (int ii=0; ii < nv; ++ii) {
    f0->bx = bx;
    f0->by = by;
    f0->bz = bz;
    ++f0;
  }
  // Array indexing version
  //for (int kk=0; kk < (nz+2*ng); ++kk) {
  //  for (int jj=0; jj < (ny+2*ng); ++jj) {
  //    for (int ii=0; ii < (nx+2*ng); ++ii) {
  //      f0 = voxel(ii,jj,kk);
  //      f0->bx = bx;
  //      f0->by = by;
  //      f0->bz = bz;
  //    }
  //  }
  //}
}

// Set E field to a uniform value
void FieldArray::uniform_e(float ex, float ey, float ez) {
  field_t* f0 = f;
  // Pointer arithmetic version
  for (int ii=0; ii < nv; ++ii) {
    f0->ex = ex;
    f0->ey = ey;
    f0->ez = ez;
    ++f0;
  }
  // Array indexing version
  //for (int kk=0; kk < (nz+2*ng); ++kk) {
  //  for (int jj=0; jj < (ny+2*ng); ++jj) {
  //    for (int ii=0; ii < (nx+2*ng); ++ii) {
  //      f0 = voxel(ii,jj,kk);
  //      f0->ex = ex;
  //      f0->ey = ey;
  //      f0->ez = ez;
  //    }
  //  }
  //}
}
