#include <stdlib.h>  // for malloc

#include "fields.h"

// Constructor
FieldArray::FieldArray(int nx, int ny, int nz, int ng) {
  int nv = (nx+2*ng) * (ny+2*ng) * (nz+2*ng);
  f = (field_t*) malloc( nv*sizeof(field_t) );
  _nx = nx;
  _ny = ny;
  _nz = nz;
  _ng = ng;
  _nv = nv;
}

// Access fields using 3D integer tuple of grid mesh indices, 0-indexed
int FieldArray::ivoxel(int ii, int jj, int kk) {
  return ii + (_nx+2*_ng)*(jj + (_ny+2*_ng)*kk);
}
field_t* FieldArray::voxel(int ii, int jj, int kk) {
  return &f[ ivoxel(ii,jj,kk) ];
}

// Set B field to a uniform value
void FieldArray::uniform_b(float bx, float by, float bz) {
  field_t* f0;
  for (int kk=0; kk < (_nz+2*_ng); ++kk) {
    for (int jj=0; jj < (_ny+2*_ng); ++jj) {
      for (int ii=0; ii < (_nx+2*_ng); ++ii) {
        f0 = voxel(ii,jj,kk);
        f0->bx = bx;
        f0->by = by;
        f0->bz = bz;
      }
    }
  }
}

// Set E field to a uniform value
void FieldArray::uniform_e(float ex, float ey, float ez) {
  field_t* f0;
  for (int kk=0; kk < (_nz+2*_ng); ++kk) {
    for (int jj=0; jj < (_ny+2*_ng); ++jj) {
      for (int ii=0; ii < (_nx+2*_ng); ++ii) {
        f0 = voxel(ii,jj,kk);
        f0->ex = ex;
        f0->ey = ey;
        f0->ez = ez;
      }
    }
  }
}
