#include <assert.h>
#include <stdlib.h>  // for malloc, qsort
#include <stdio.h>  // for printf

#include "field.h"

// For field-init ghost cell testing. -ATr,2024feb10
// copy-pasted from https://en.cppreference.com/w/c/algorithm/qsort
int compare_ints(const void* a, const void* b) {
    int arg1 = *(const int*)a;
    int arg2 = *(const int*)b;
    if (arg1 < arg2) return -1;
    if (arg1 > arg2) return 1;
    return 0;
}

// Constructor
FieldArray::FieldArray(int nx_, int ny_, int nz_, int ng_,
                       float hx_, float hy_, float hz_, float dt_) {
  nx = nx_;
  ny = ny_;
  nz = nz_;
  ng = ng_;
  nvall = (nx+2*ng) * (ny+2*ng) * (nz+2*ng);
  nvg   = nvall - (nx*ny*nz);

  hx = hx_;
  hy = hy_;
  hz = hz_;
  dt = dt_;

  // default boundary condition is periodic unless user overrides
  particle_bc_x = 0;
  particle_bc_y = 0;
  particle_bc_z = 0;
  field_bc_x = 0;
  field_bc_y = 0;
  field_bc_z = 0;

  // assumes float = 4 bytes, not performance portable but it's such a common
  // standard that I think this check is worthwhile
  assert(sizeof(field_t) == 4*nfstruct);

  f0 = (field_t*) malloc( nvall*sizeof(field_t) );

  // ghost cell management
  //
  // The cell indexing scheme is:
  //     ghost = [    0,      ng-1] = [    0,      ng)
  //     live  = [   ng, nx+  ng-1] = [   ng, nx+  ng)
  //     ghost = [nx+ng, nx+2*ng-1] = [nx+ng, nx+2*ng)
  // where bracket=inclusive, parens=exclusive range bounds.
  // Example: for 10 live and 2 ghost cells,
  // indices 0-1 ghost, 2-11 live, 12-13 ghost.
  //
  // Construct two arrays of ghost cells' linear voxel indices, and the ghosts'
  // matched "source" cell voxel indices.  You must manage...
  // 6 slabs, total cell count nx*ny*(2*ng) + nx*(2*ng)*nz + (2*ng)*ny*nz
  // 8 corners, total cell count 8*ng^3
  // 12 edges, total cell count 4*nx*ng^2 + 4*ny*ng^2 + 4*nz*ng^2
  //
  // Beware "source" cells will appear multiple times (for example, a 2x2x2
  // domain with 1 ghost cell/side has 64 cells = 8 live, 56 ghost).
  // slab count 24, corner count 8, edge count 24 in this example.

  ivoxels_ghost = (int*) malloc(nvg*sizeof(int));
  ivoxels_ghsrc = (int*) malloc(nvg*sizeof(int));

  int count = 0;
  int* ivghost = ivoxels_ghost;
  int* ivghsrc = ivoxels_ghsrc;

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ghost slabs

  // Left x slab <- right; (0,y,z) <- (1,y,z)
  for (int kk=(   ng); kk < (nz+ng); ++kk) {
  for (int jj=(   ng); jj < (ny+ng); ++jj) {
  for (int ii=(    0); ii < (   ng); ++ii) {
    int ig0 = ii;
    int is1 = (nx+ng-1) - (ng-ii-1) % nx;  // resolves to nx+ii for ng <= nx
    //if ((jj==ng)&&(kk==ng)) printf("Left x slab %d <- %d\n", ig0, is1);
    *ivghost = ivoxel( ig0,  jj,  kk);
    *ivghsrc = ivoxel( is1,  jj,  kk);
    ++ivghost;
    ++ivghsrc;
    ++count;
  }}}
  // Left y slab <- right; (x,0,z) <- (x,1,z)
  for (int kk=(   ng); kk < (nz+ng); ++kk) {
  for (int jj=(    0); jj < (   ng); ++jj) {
  for (int ii=(   ng); ii < (nx+ng); ++ii) {
    int jg0 = jj;
    int js1 = (ny+ng-1) - (ng-jj-1) % ny;  // resolves to ny+jj for ng <= ny
    //if ((ii==ng)&&(kk==ng)) printf("Left y slab %d <- %d\n", jg0, js1);
    *ivghost = ivoxel(  ii, jg0,  kk);
    *ivghsrc = ivoxel(  ii, js1,  kk);
    ++ivghost;
    ++ivghsrc;
    ++count;
  }}}
  // Left z slab <- right; (x,y,0) <- (x,y,1)
  for (int kk=(    0); kk < (   ng); ++kk) {
  for (int jj=(   ng); jj < (ny+ng); ++jj) {
  for (int ii=(   ng); ii < (nx+ng); ++ii) {
    int kg0 = kk;
    int ks1 = (nz+ng-1) - (ng-kk-1) % nz;  // resolves to nz+kk for ng <= nz
    //if ((ii==ng)&&(jj==ng)) printf("Left z slab %d <- %d\n", kg0, ks1);
    *ivghost = ivoxel(  ii,  jj, kg0);
    *ivghsrc = ivoxel(  ii,  jj, ks1);
    ++ivghost;
    ++ivghsrc;
    ++count;
  }}}

  // Right x slab <- left; (1,y,z) <- (0,y,z)
  for (int kk=(   ng); kk < (nz+ng); ++kk) {
  for (int jj=(   ng); jj < (ny+ng); ++jj) {
  for (int ii=(    0); ii < (   ng); ++ii) {
    int ig1 = nx+ng+ii;
    int is0 = ng + ii % nx;
    //if ((jj==ng)&&(kk==ng)) printf("Right x slab %d <- %d\n", ig1, is0);
    *ivghost = ivoxel( ig1,  jj,  kk);
    *ivghsrc = ivoxel( is0,  jj,  kk);
    ++ivghost;
    ++ivghsrc;
    ++count;
  }}}
  // Right y slab <- left; (x,1,z) <- (x,0,z)
  for (int kk=(   ng); kk < (nz+ng); ++kk) {
  for (int jj=(    0); jj < (   ng); ++jj) {
  for (int ii=(   ng); ii < (nx+ng); ++ii) {
    int jg1 = ny+ng+jj;
    int js0 = ng + jj % ny;
    //if ((ii==ng)&&(kk==ng)) printf("Right y slab %d <- %d\n", jg1, js0);
    *ivghost = ivoxel(  ii, jg1,  kk);
    *ivghsrc = ivoxel(  ii, js0,  kk);
    ++ivghost;
    ++ivghsrc;
    ++count;
  }}}
  // Right z slab <- left; (x,y,1) <- (x,y,0)
  for (int kk=(    0); kk < (   ng); ++kk) {
  for (int jj=(   ng); jj < (ny+ng); ++jj) {
  for (int ii=(   ng); ii < (nx+ng); ++ii) {
    int kg1 = nz+ng+kk;
    int ks0 = ng + kk % nz;
    //if ((ii==ng)&&(jj==ng)) printf("Right z slab %d <- %d\n", kg1, ks0);
    *ivghost = ivoxel(  ii,  jj, kg1);
    *ivghsrc = ivoxel(  ii,  jj, ks0);
    ++ivghost;
    ++ivghsrc;
    ++count;
  }}}

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ghost corners
  // (x,y,z) = (0,0,0) = lower left
  // (x,y,z) = (1,1,1) = upper right
  // for now, don't worry about locality of voxel indices
  // ... b/c the cost of ghost cell ops should be relatively small
  // code is designed for block copy/paste during development

  for (int kk=0; kk < ng; ++kk) {
  for (int jj=0; jj < ng; ++jj) {
  for (int ii=0; ii < ng; ++ii) {

    // left(0) ghosts <- right(1) source cells
    int ig0 = ii;
    int jg0 = jj;
    int kg0 = kk;
    int is1 = (nx+ng-1) - (ng-ii-1) % nx;  // resolves to nx+ii for ng <= nx
    int js1 = (ny+ng-1) - (ng-jj-1) % ny;  // resolves to ny+jj for ng <= ny
    int ks1 = (nz+ng-1) - (ng-kk-1) % nz;  // resolves to nz+kk for ng <= nz
    // right(1) ghosts <- left(0) source cells
    int ig1 = nx+ng+ii;
    int jg1 = ny+ng+jj;
    int kg1 = nz+ng+kk;
    int is0 = ng + ii % nx;
    int js0 = ng + jj % ny;
    int ks0 = ng + kk % nz;

    // (0,0,0) <- (1,1,1)
    //printf("(0,0,0) corner %d,%d,%d <- %d,%d,%d\n", ig0,jg0,kg0, is1,js1,ks1);
    *ivghost = ivoxel( ig0, jg0, kg0 );
    *ivghsrc = ivoxel( is1, js1, ks1 );
    ++ivghost;
    ++ivghsrc;
    ++count;

    // (1,0,0) <- (0,1,1)
    //printf("(1,0,0) corner %d,%d,%d <- %d,%d,%d\n", ig1,jg0,kg0, is0,js1,ks1);
    *ivghost = ivoxel( ig1, jg0, kg0 );
    *ivghsrc = ivoxel( is0, js1, ks1 );
    ++ivghost;
    ++ivghsrc;
    ++count;

    // (0,1,0) <- (1,0,1)
    //printf("(0,1,0) corner %d,%d,%d <- %d,%d,%d\n", ig0,jg1,kg0, is1,js0,ks1);
    *ivghost = ivoxel( ig0, jg1, kg0 );
    *ivghsrc = ivoxel( is1, js0, ks1 );
    ++ivghost;
    ++ivghsrc;
    ++count;

    // (1,1,0) <- (0,0,1)
    //printf("(1,1,0) corner %d,%d,%d <- %d,%d,%d\n", ig1,jg1,kg0, is0,js0,ks1);
    *ivghost = ivoxel( ig1, jg1, kg0 );
    *ivghsrc = ivoxel( is0, js0, ks1 );
    ++ivghost;
    ++ivghsrc;
    ++count;

    // (0,0,1) <- (1,1,0)
    //printf("(0,0,1) corner %d,%d,%d <- %d,%d,%d\n", ig0,jg0,kg1, is1,js1,ks0);
    *ivghost = ivoxel( ig0, jg0, kg1 );
    *ivghsrc = ivoxel( is1, js1, ks0 );
    ++ivghost;
    ++ivghsrc;
    ++count;

    // (1,0,1) <- (0,1,0)
    //printf("(1,0,1) corner %d,%d,%d <- %d,%d,%d\n", ig1,jg0,kg1, is0,js1,ks0);
    *ivghost = ivoxel( ig1, jg0, kg1 );
    *ivghsrc = ivoxel( is0, js1, ks0 );
    ++ivghost;
    ++ivghsrc;
    ++count;

    // (0,1,1) <- (1,0,0)
    //printf("(0,1,1) corner %d,%d,%d <- %d,%d,%d\n", ig0,jg1,kg1, is1,js0,ks0);
    *ivghost = ivoxel( ig0, jg1, kg1 );
    *ivghsrc = ivoxel( is1, js0, ks0 );
    ++ivghost;
    ++ivghsrc;
    ++count;

    // (1,1,1) <- (0,0,0)
    //printf("(1,1,1) corner %d,%d,%d <- %d,%d,%d\n", ig1,jg1,kg1, is0,js0,ks0);
    *ivghost = ivoxel( ig1, jg1, kg1 );
    *ivghsrc = ivoxel( is0, js0, ks0 );
    ++ivghost;
    ++ivghsrc;
    ++count;
  }}}

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ghost edges

  // x-aligned edges
  for (int kk=(    0); kk < (   ng); ++kk) {
  for (int jj=(    0); jj < (   ng); ++jj) {
  for (int ii=(   ng); ii < (nx+ng); ++ii) {

    // left(0) ghosts <- right(1) source cells
    int jg0 = jj;
    int kg0 = kk;
    int js1 = (ny+ng-1) - (ng-jj-1) % ny;  // resolves to ny+jj for ng <= ny
    int ks1 = (nz+ng-1) - (ng-kk-1) % nz;  // resolves to nz+kk for ng <= nz
    // right(1) ghosts <- left(0) source cells
    int jg1 = ny+ng+jj;
    int kg1 = nz+ng+kk;
    int js0 = ng + jj % ny;
    int ks0 = ng + kk % nz;

    // (x,0,0) <- (x,1,1)
    *ivghost = ivoxel(  ii, jg0, kg0 );
    *ivghsrc = ivoxel(  ii, js1, ks1 );
    ++ivghost;
    ++ivghsrc;
    ++count;
    // (x,1,0) <- (x,0,1)
    *ivghost = ivoxel(  ii, jg1, kg0 );
    *ivghsrc = ivoxel(  ii, js0, ks1 );
    ++ivghost;
    ++ivghsrc;
    ++count;
    // (x,0,1) <- (x,1,0)
    *ivghost = ivoxel(  ii, jg0, kg1 );
    *ivghsrc = ivoxel(  ii, js1, ks0 );
    ++ivghost;
    ++ivghsrc;
    ++count;
    // (x,1,1) <- (x,0,0)
    *ivghost = ivoxel(  ii, jg1, kg1 );
    *ivghsrc = ivoxel(  ii, js0, ks0 );
    ++ivghost;
    ++ivghsrc;
    ++count;
  }}}

  // y-aligned edges
  for (int kk=(    0); kk < (   ng); ++kk) {
  for (int jj=(   ng); jj < (ny+ng); ++jj) {
  for (int ii=(    0); ii < (   ng); ++ii) {

    // left(0) ghosts <- right(1) source cells
    int ig0 = ii;
    int kg0 = kk;
    int is1 = (nx+ng-1) - (ng-ii-1) % nx;  // resolves to nx+ii for ng <= nx
    int ks1 = (nz+ng-1) - (ng-kk-1) % nz;  // resolves to nz+kk for ng <= nz
    // right(1) ghosts <- left(0) source cells
    int ig1 = nx+ng+ii;
    int kg1 = nz+ng+kk;
    int is0 = ng + ii % nx;
    int ks0 = ng + kk % nz;

    // (0,y,0) <- (1,y,1)
    *ivghost = ivoxel( ig0,  jj, kg0);
    *ivghsrc = ivoxel( is1,  jj, ks1);
    ++ivghost;
    ++ivghsrc;
    ++count;
    // (1,y,0) <- (0,y,1)
    *ivghost = ivoxel( ig1,  jj, kg0);
    *ivghsrc = ivoxel( is0,  jj, ks1);
    ++ivghost;
    ++ivghsrc;
    ++count;
    // (0,y,1) <- (1,y,0)
    *ivghost = ivoxel( ig0,  jj, kg1);
    *ivghsrc = ivoxel( is1,  jj, ks0);
    ++ivghost;
    ++ivghsrc;
    ++count;
    // (1,y,1) <- (0,y,0)
    *ivghost = ivoxel( ig1,  jj, kg1);
    *ivghsrc = ivoxel( is0,  jj, ks0);
    ++ivghost;
    ++ivghsrc;
    ++count;
  }}}

  // z-aligned edges
  for (int kk=(   ng); kk < (nz+ng); ++kk) {
  for (int jj=(    0); jj < (   ng); ++jj) {
  for (int ii=(    0); ii < (   ng); ++ii) {

    // left(0) ghosts <- right(1) source cells
    int ig0 = ii;
    int jg0 = jj;
    int is1 = (nx+ng-1) - (ng-ii-1) % nx;  // resolves to nx+ii for ng <= nx
    int js1 = (ny+ng-1) - (ng-jj-1) % ny;  // resolves to ny+jj for ng <= ny
    // right(1) ghosts <- left(0) source cells
    int ig1 = nx+ng+ii;
    int jg1 = ny+ng+jj;
    int is0 = ng + ii % nx;
    int js0 = ng + jj % ny;

    // (0,0,z) <- (1,1,z)
    *ivghost = ivoxel( ig0, jg0,  kk);
    *ivghsrc = ivoxel( is1, js1,  kk);
    ++ivghost;
    ++ivghsrc;
    ++count;
    // (1,0,z) <- (0,1,z)
    *ivghost = ivoxel( ig1, jg0,  kk);
    *ivghsrc = ivoxel( is0, js1,  kk);
    ++ivghost;
    ++ivghsrc;
    ++count;
    // (0,1,z) <- (1,0,z)
    *ivghost = ivoxel( ig0, jg1,  kk);
    *ivghsrc = ivoxel( is1, js0,  kk);
    ++ivghost;
    ++ivghsrc;
    ++count;
    // (1,1,z) <- (0,0,z)
    *ivghost = ivoxel( ig1, jg1,  kk);
    *ivghsrc = ivoxel( is0, js0,  kk);
    ++ivghost;
    ++ivghsrc;
    ++count;
  }}}

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // all done counting ghosts, do some sanity checks

  // Need sorted ivoxels for some tests.
  // DO NOT USE SORTED LISTS IN PRODUCTION CODE, because the relative order
  // of ivoxels_ghost, ivoxels_ghsrc must be maintained.  --ATr,2024feb10
  int* ivoxels_test_g = (int*) malloc(nvg*sizeof(int));
  int* ivoxels_test_s = (int*) malloc(nvg*sizeof(int));
  for (int ii =0; ii < count; ++ii) {
    ivoxels_test_g[ii] = ivoxels_ghost[ii];
    ivoxels_test_s[ii] = ivoxels_ghsrc[ii];
  }
  qsort(ivoxels_test_g, count, sizeof(int), compare_ints);
  qsort(ivoxels_test_s, count, sizeof(int), compare_ints);

  printf("Field initialization: testing ghost count %d expected %d ghosts.\n", count, nvg);
  assert(count == nvg);

  printf("Field initialization: testing ghost linear voxel indices unique.\n");
  //printf("linear voxel unsorted %d sorted %d\n", ivoxels_ghost[0], ivoxels_test_g[0]);
  for (int ii = 1; ii < count; ++ii) {
    //printf("linear voxel unsorted %d sorted %d\n", ivoxels_ghost[ii], ivoxels_test_g[ii]);
    assert(ivoxels_test_g[ii-1] != ivoxels_test_g[ii]);
  }

  // Useful test when any dimension < number of ghost cells, i.e.,
  // min(nx,ny,nz) < ng so the ghost cells have to alias around
  printf("Field initialization: testing ghosts only point to live cells.\n");
  {
    int iig = 0;
    int iis = 0;
    // Fail if any src voxels are also ghost voxels
    // Simultaneous loop requires both lists to be sorted
    while (iis < count) {
      int ivg = ivoxels_test_g[iig];
      int ivs = ivoxels_test_s[iis];
      while (ivg < ivs && iig < count) {
        // ghost iterator catches up to source when ivg >= ivs (note "=")
        ++iig;
        ivg = ivoxels_test_g[iig];
      }
      assert(ivs != ivg);
      ++iis;
    }
  }

  // Cleanup from testing
  free(ivoxels_test_g);
  free(ivoxels_test_s);

  printf("Field initialization: all tests passed.\n");

}
