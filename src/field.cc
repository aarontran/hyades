#include <assert.h>
#include <stdlib.h>  // for NULL
#include <stdio.h>  // for printf
#include <string.h>  // for strcmp

#include "field.h"

// ----------------------------------------------------------------------------
// Low-level methods

// Access fields using 3D integer tuple of grid mesh indices, 0-indexed
int FieldArray::ivoxel(int ii, int jj, int kk) {
  return ii + (nx+2*ng)*(jj + (ny+2*ng)*kk);
}

// Access fields using 3D integer tuple of grid mesh indices, 0-indexed
field_t* FieldArray::voxel(int ii, int jj, int kk) {
  return &( f0[ivoxel(ii,jj,kk)] );
}

// Low-level method to access field_t struct members by name,
// returning a pointer (called "seek" rather than "get" because caller
// can both read and write).
float* FieldArray::fseek_key(const char* name, field_t* ff) {
  if (strcmp(   "ex",name)==0) return &(ff->ex);
  if (strcmp(   "ey",name)==0) return &(ff->ey);
  if (strcmp(   "ez",name)==0) return &(ff->ez);
  if (strcmp( "tmpx",name)==0) return &(ff->tmpx);
  if (strcmp(   "bx",name)==0) return &(ff->bx);
  if (strcmp(   "by",name)==0) return &(ff->by);
  if (strcmp(   "bz",name)==0) return &(ff->bz);
  if (strcmp( "tmpy",name)==0) return &(ff->tmpy);
  if (strcmp(  "bx0",name)==0) return &(ff->bx0);
  if (strcmp(  "by0",name)==0) return &(ff->by0);
  if (strcmp(  "bz0",name)==0) return &(ff->bz0);
  if (strcmp( "tmpz",name)==0) return &(ff->tmpz);
  if (strcmp(  "jfx",name)==0) return &(ff->jfx);
  if (strcmp(  "jfy",name)==0) return &(ff->jfy);
  if (strcmp(  "jfz",name)==0) return &(ff->jfz);
  if (strcmp( "rhof",name)==0) return &(ff->rhof);
  if (strcmp( "jfx0",name)==0) return &(ff->jfx0);
  if (strcmp( "jfy0",name)==0) return &(ff->jfy0);
  if (strcmp( "jfz0",name)==0) return &(ff->jfz0);
  if (strcmp("rhof0",name)==0) return &(ff->rhof0);
  if (strcmp( "smex",name)==0) return &(ff->smex);
  if (strcmp( "smey",name)==0) return &(ff->smey);
  if (strcmp( "smez",name)==0) return &(ff->smez);
  if (strcmp(  "pad",name)==0) return &(ff->pad);
  assert(0);
}

// Low-level method to access field_t struct members by linear index,
// returning a pointer (called "seek" rather than "get" because caller
// can both read and write).
float* FieldArray::fseek_one(int mm, field_t* ff) {
  switch(mm) {
    case  0: return &(ff->ex);
    case  1: return &(ff->ey);
    case  2: return &(ff->ez);
    case  3: return &(ff->tmpx);
    case  4: return &(ff->bx);
    case  5: return &(ff->by);
    case  6: return &(ff->bz);
    case  7: return &(ff->tmpy);
    case  8: return &(ff->bx0);
    case  9: return &(ff->by0);
    case 10: return &(ff->bz0);
    case 11: return &(ff->tmpz);
    case 12: return &(ff->jfx);
    case 13: return &(ff->jfy);
    case 14: return &(ff->jfz);
    case 15: return &(ff->rhof);
    case 16: return &(ff->jfx0);
    case 17: return &(ff->jfy0);
    case 18: return &(ff->jfz0);
    case 19: return &(ff->rhof0);
    case 20: return &(ff->smex);
    case 21: return &(ff->smey);
    case 22: return &(ff->smez);
    case 23: return &(ff->pad);
    default: return NULL;
  }
}

// Low-level method to set ghost and live cell values for one field
void FieldArray::fset_one(int mm, float v, field_t* ff) {
  *(fseek_one(mm,ff)) = v;
}

// Low-level method to set ghost and live cell values for all fields
void FieldArray::fset_all(float v, field_t* ff) {
  for (int mm = 0; mm < nfstruct; ++mm) {
    fset_one(mm, v, ff);
  }
}

// Low-level method to copy one field value from one voxel to another
// ii = field struct member index, ff = destination voxel, gg = source voxel
void FieldArray::fcopy_one(int mm, field_t* ff, field_t* gg) {
  *(fseek_one(mm,ff)) = *(fseek_one(mm,gg));
}

// Low-level method to copy all field values from one voxel to another
// ff = destination voxel, gg = source voxel
void FieldArray::fcopy_all(field_t* ff, field_t* gg) {
  for (int mm = 0; mm < nfstruct; ++mm) {
    fcopy_one(mm, ff, gg);
  }
}

// ----------------------------------------------------------------------------
// Intermediate-level full/alive-mesh methods

// Set one field struct member mm to value v on full mesh (live+ghost)
//void FieldArray::mesh_set_one(int mm, float v) {
//  for (int kk = 0; kk < (nz+2*ng); ++kk) {
//  for (int jj = 0; jj < (ny+2*ng); ++jj) {
//  for (int ii = 0; ii < (nx+2*ng); ++ii) {
//    field_t* ff = voxel(ii,jj,kk);
//    fset_one(mm, v, ff);
//  }}}
//}

// Set all field struct members to value v on full mesh (live+ghost)
void FieldArray::mesh_set_all(float v) {
  for (int kk = 0; kk < (nz+2*ng); ++kk) {
  for (int jj = 0; jj < (ny+2*ng); ++jj) {
  for (int ii = 0; ii < (nx+2*ng); ++ii) {
    field_t* ff = voxel(ii,jj,kk);
    fset_all(v, ff);
  }}}
}

// Set one field struct member mm to value v on live cells
//void FieldArray::alive_set_one(int mm, float v) {
//  for (int kk = ng; kk < (nz+ng); ++kk) {
//  for (int jj = ng; jj < (ny+ng); ++jj) {
//  for (int ii = ng; ii < (nx+ng); ++ii) {
//    field_t* ff = voxel(ii,jj,kk);
//    fset_one(mm, v, ff);
//  }}}
//}

// Set all field struct members to value v on live cells
//void FieldArray::alive_set_all(float v) {
//  for (int kk = ng; kk < (nz+ng); ++kk) {
//  for (int jj = ng; jj < (ny+ng); ++jj) {
//  for (int ii = ng; ii < (nx+ng); ++ii) {
//    field_t* ff = voxel(ii,jj,kk);
//    fset_all(v, ff);
//  }}}
//}

// ----------------------------------------------------------------------------
// Intermediate-level ghost-mesh methods

// Copy all field struct members from live cells to ghost cells
void FieldArray::ghost_copy_all() {
  int* ivghost = ivoxels_ghost;
  int* ivghsrc = ivoxels_ghsrc;
  for (int ii = 0; ii < nvg; ++ii) {
    field_t* gh  = &( f0[*ivghost] );
    field_t* src = &( f0[*ivghsrc] );
    fcopy_all(gh, src);
    ++ivghost;
    ++ivghsrc;
  }
}

// ----------------------------------------------------------------------------
// High-level methods

// Set B field to a uniform value, including ghost cells
void FieldArray::mesh_set_b(float bx, float by, float bz) {
  for (int kk = 0; kk < (nz+2*ng); ++kk) {
  for (int jj = 0; jj < (ny+2*ng); ++jj) {
  for (int ii = 0; ii < (nx+2*ng); ++ii) {
    field_t* ff = voxel(ii,jj,kk);
    ff->bx = bx;
    ff->by = by;
    ff->bz = bz;
  }}}
}

// Set E field to a uniform value, including ghost cells
void FieldArray::mesh_set_e(float ex, float ey, float ez) {
  for (int kk = 0; kk < (nz+2*ng); ++kk) {
  for (int jj = 0; jj < (ny+2*ng); ++jj) {
  for (int ii = 0; ii < (nx+2*ng); ++ii) {
    field_t* ff = voxel(ii,jj,kk);
    ff->ex = ex;
    ff->ey = ey;
    ff->ez = ez;
  }}}
}

// Set j/rho fields to a uniform value, including ghost cells
void FieldArray::mesh_set_jrho(float jx, float jy, float jz, float rho) {
  for (int kk = 0; kk < (nz+2*ng); ++kk) {
  for (int jj = 0; jj < (ny+2*ng); ++jj) {
  for (int ii = 0; ii < (nx+2*ng); ++ii) {
    field_t* ff = voxel(ii,jj,kk);
    ff->jfx = jx;
    ff->jfy = jy;
    ff->jfz = jz;
    ff->rhof = rho;
  }}}
}

// Set old j/rho fields to current j/rho fields, including ghost cells
void FieldArray::mesh_set_jrho0() {
  for (int kk = 0; kk < (nz+2*ng); ++kk) {
  for (int jj = 0; jj < (ny+2*ng); ++jj) {
  for (int ii = 0; ii < (nx+2*ng); ++ii) {
    field_t* ff = voxel(ii,jj,kk);
    ff->jfx0  = ff->jfx;
    ff->jfy0  = ff->jfy;
    ff->jfz0  = ff->jfz;
    ff->rhof0 = ff->rhof;
  }}}
}

// Copy E/B fields from live cells to ghost cells; overwrite old ghosts
void FieldArray::ghost_copy_eb() {
  int* ivghost = ivoxels_ghost;
  int* ivghsrc = ivoxels_ghsrc;
  for (int ii = 0; ii < nvg; ++ii) {
    field_t* gh  = &( f0[*ivghost] );
    field_t* src = &( f0[*ivghsrc] );
    gh->ex = src->ex;
    gh->ey = src->ey;
    gh->ez = src->ez;
    gh->bx = src->bx;
    gh->by = src->by;
    gh->bz = src->bz;
    ++ivghost;
    ++ivghsrc;
  }
}

// Copy E field from live cells to ghost cells; overwrite old ghosts
void FieldArray::ghost_copy_e() {
  int* ivghost = ivoxels_ghost;
  int* ivghsrc = ivoxels_ghsrc;
  for (int ii = 0; ii < nvg; ++ii) {
    field_t* gh  = &( f0[*ivghost] );
    field_t* src = &( f0[*ivghsrc] );
    gh->ex = src->ex;
    gh->ey = src->ey;
    gh->ez = src->ez;
    ++ivghost;
    ++ivghsrc;
  }
}

// Copy B field from live cells to ghost cells; overwrite old ghosts
void FieldArray::ghost_copy_b() {
  int* ivghost = ivoxels_ghost;
  int* ivghsrc = ivoxels_ghsrc;
  for (int ii = 0; ii < nvg; ++ii) {
    field_t* gh  = &( f0[*ivghost] );
    field_t* src = &( f0[*ivghsrc] );
    gh->bx = src->bx;
    gh->by = src->by;
    gh->bz = src->bz;
    ++ivghost;
    ++ivghsrc;
  }
}

// Copy j/rho fields from live cells to ghost cells; overwrite old ghosts
void FieldArray::ghost_copy_jrho() {
  int* ivghost = ivoxels_ghost;
  int* ivghsrc = ivoxels_ghsrc;
  for (int ii = 0; ii < nvg; ++ii) {
    field_t* gh  = &( f0[*ivghost] );
    field_t* src = &( f0[*ivghsrc] );
    gh->jfx  = src->jfx;
    gh->jfy  = src->jfy;
    gh->jfz  = src->jfz;
    gh->rhof = src->rhof;
    ++ivghost;
    ++ivghsrc;
  }
}

// Add all ghost cell values of (j, rho) to existing values on grid.
void FieldArray::ghost_deposit_jrho() {
  int* ivghost = ivoxels_ghost;
  int* ivghsrc = ivoxels_ghsrc;
  for (int ii = 0; ii < nvg; ++ii) {
    field_t* gh  = &( f0[*ivghost] );
    field_t* src = &( f0[*ivghsrc] );
    src->jfx  += gh-> jfx;
    src->jfy  += gh-> jfy;
    src->jfz  += gh-> jfz;
    src->rhof += gh-> rhof;
    ++ivghost;
    ++ivghsrc;
  }
}
