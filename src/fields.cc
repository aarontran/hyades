#include <assert.h>
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

  // assumes float = 4 bytes, not performance portable but it's such a common
  // standard that I think this check is worthwhile
  assert(sizeof(field_t) == 4*nfstruct);

  f0 = (field_t*) malloc( nv*sizeof(field_t) );
}

// ----------------------------------------------------------------------------
// Low-level methods

// Access fields using 3D integer tuple of grid mesh indices, 0-indexed
int FieldArray::ivoxel(int ii, int jj, int kk) {
  return ii + (nx+2*ng)*(jj + (ny+2*ng)*kk);
}

// Access fields using 3D integer tuple of grid mesh indices, 0-indexed
field_t* FieldArray::voxel(int ii, int jj, int kk) {
  return &f0[ ivoxel(ii,jj,kk) ];
}

// Low-level method to access over field_t struct members
float* FieldArray::fget(field_t* ff, int ii) {
  switch(ii) {
    case  0: return &(ff->ex);
    case  1: return &(ff->ey);
    case  2: return &(ff->ez);
    case  3: return &(ff->epsx);
    case  4: return &(ff->bx);
    case  5: return &(ff->by);
    case  6: return &(ff->bz);
    case  7: return &(ff->epsy);
    case  8: return &(ff->bx0);
    case  9: return &(ff->by0);
    case 10: return &(ff->bz0);
    case 11: return &(ff->epsz);
    case 12: return &(ff->jfx);
    case 13: return &(ff->jfy);
    case 14: return &(ff->jfz);
    case 15: return &(ff->rhof);
    default: return NULL;
  }
}

// Low-level method to copy field_t struct members from one voxel to another
// ff = destination, gg = data to copy
void FieldArray::fcopy(field_t* ff, field_t* gg) {
  for (int ii = 0; ii < nfstruct; ++ii) {
    *(fget(ff,ii)) = *(fget(gg,ii));
  }
}

// Low-level method to set ghost and live cell values for all fields
// (including ghost corners), meant for testing and debugging.
void FieldArray::freset(float v) {
  field_t* ff = f0;
  for (int ii=0; ii < nv; ++ii) {
    for (int nn = 0; nn < nfstruct; ++nn) {
      *(fget(ff,nn)) = v;
    }
    ++ff;
  }
}

// ----------------------------------------------------------------------------
// High-level methods

// Set B field to a uniform value, EXCLUDING ghost cells
void FieldArray::uniform_b(float bx, float by, float bz) {
  for (int kk=ng; kk < (nz+ng); ++kk) {
    for (int jj=ng; jj < (ny+ng); ++jj) {
      for (int ii=ng; ii < (nx+ng); ++ii) {
        field_t* ff = voxel(ii,jj,kk);
        ff->bx = bx;
        ff->by = by;
        ff->bz = bz;
      }
    }
  }
}

// Set E field to a uniform value, EXCLUDING ghost cells
void FieldArray::uniform_e(float ex, float ey, float ez) {
  for (int kk=ng; kk < (nz+ng); ++kk) {
    for (int jj=ng; jj < (ny+ng); ++jj) {
      for (int ii=ng; ii < (nx+ng); ++ii) {
        field_t* ff = voxel(ii,jj,kk);
        ff->ex = ex;
        ff->ey = ey;
        ff->ez = ez;
      }
    }
  }
}

// Copy ghost cell values
void FieldArray::update_ghost() {
  field_t* local;
  field_t* remote;

  // Start by updating slabs

  // Left x slab <- right, no corners
//  for (int kk=ng; kk < (nz+ng); ++kk) {
//  for (int jj=ng; jj < (ny+ng); ++jj) {
//  for (int ii= 0; ii <     ng ; ++ii) {
//    local  = voxel(ii,             jj, kk);
//    remote = voxel((nx+ng-1)-ii,   jj, kk);
//    fcpy(local, remote);
//  }}}

}
