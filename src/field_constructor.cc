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
    //if ((jj==ng)&&(kk==ng)) printf("Left x slab %d <- %d\n", ii, nx+ii);
    *ivghost = ivoxel(      ii,       jj,       kk);
    *ivghsrc = ivoxel(   nx+ii,       jj,       kk);
    ++ivghost;
    ++ivghsrc;
    ++count;
  }}}
  // Left y slab <- right; (x,0,z) <- (x,1,z)
  for (int kk=(   ng); kk < (nz+ng); ++kk) {
  for (int jj=(    0); jj < (   ng); ++jj) {
  for (int ii=(   ng); ii < (nx+ng); ++ii) {
    //if ((ii==ng)&&(kk==ng)) printf("Left y slab %d <- %d\n", jj, ny+jj);
    *ivghost = ivoxel(      ii,       jj,       kk);
    *ivghsrc = ivoxel(      ii,    ny+jj,       kk);
    ++ivghost;
    ++ivghsrc;
    ++count;
  }}}
  // Left z slab <- right; (x,y,0) <- (x,y,1)
  for (int kk=(    0); kk < (   ng); ++kk) {
  for (int jj=(   ng); jj < (ny+ng); ++jj) {
  for (int ii=(   ng); ii < (nx+ng); ++ii) {
    //if ((ii==ng)&&(jj==ng)) printf("Left z slab %d <- %d\n", kk, nz+kk);
    *ivghost = ivoxel(      ii,       jj,       kk);
    *ivghsrc = ivoxel(      ii,       jj,    nz+kk);
    ++ivghost;
    ++ivghsrc;
    ++count;
  }}}

  // Right x slab <- left; (1,y,z) <- (0,y,z)
  for (int kk=(   ng); kk < (nz+ng); ++kk) {
  for (int jj=(   ng); jj < (ny+ng); ++jj) {
  for (int ii=(    0); ii < (   ng); ++ii) {
    //if ((jj==ng)&&(kk==ng)) printf("Right x slab %d <- %d\n", nx+ng+ii, ng+ii);
    *ivghost = ivoxel(nx+ng+ii,       jj,       kk);
    *ivghsrc = ivoxel(   ng+ii,       jj,       kk);
    ++ivghost;
    ++ivghsrc;
    ++count;
  }}}
  // Right y slab <- left; (x,1,z) <- (x,0,z)
  for (int kk=(   ng); kk < (nz+ng); ++kk) {
  for (int jj=(    0); jj < (   ng); ++jj) {
  for (int ii=(   ng); ii < (nx+ng); ++ii) {
    //if ((ii==ng)&&(kk==ng)) printf("Right y slab %d <- %d\n", ny+ng+jj, ng+jj);
    *ivghost = ivoxel(      ii, ny+ng+jj,       kk);
    *ivghsrc = ivoxel(      ii,    ng+jj,       kk);
    ++ivghost;
    ++ivghsrc;
    ++count;
  }}}
  // Right z slab <- left; (x,y,1) <- (x,y,0)
  for (int kk=(    0); kk < (   ng); ++kk) {
  for (int jj=(   ng); jj < (ny+ng); ++jj) {
  for (int ii=(   ng); ii < (nx+ng); ++ii) {
    //if ((ii==ng)&&(jj==ng)) printf("Right z slab %d <- %d\n", nz+ng+kk, ng+kk);
    *ivghost = ivoxel(      ii,       jj, nz+ng+kk);
    *ivghsrc = ivoxel(      ii,       jj,    ng+kk);
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
  // printf statements are deisgned for block copy/paste during code
  // ... development

  for (int kk=0; kk < ng; ++kk) {
  for (int jj=0; jj < ng; ++jj) {
  for (int ii=0; ii < ng; ++ii) {
    // (0,0,0) <- (1,1,1)
    //printf("(0,0,0) corner %d,%d,%d <- %d,%d,%d\n",
    //                              ii,       jj,       kk,
    //                           nx+ii,    ny+jj,    nz+kk
    //);
    *ivghost = ivoxel(      ii,       jj,       kk);
    *ivghsrc = ivoxel(   nx+ii,    ny+jj,    nz+kk);
    ++ivghost;
    ++ivghsrc;
    ++count;

    // (1,0,0) <- (0,1,1)
    //printf("(1,0,0) corner %d,%d,%d <- %d,%d,%d\n",
    //                        nx+ng+ii,       jj,       kk,
    //                           ng+ii,    ny+jj,    nz+kk
    //);
    *ivghost = ivoxel(nx+ng+ii,       jj,       kk);
    *ivghsrc = ivoxel(   ng+ii,    ny+jj,    nz+kk);
    ++ivghost;
    ++ivghsrc;
    ++count;

    // (0,1,0) <- (1,0,1)
    //printf("(0,1,0) corner %d,%d,%d <- %d,%d,%d\n",
    //                              ii, ny+ng+jj,       kk,
    //                           nx+ii,    ng+jj,    nz+kk
    //);
    *ivghost = ivoxel(      ii, ny+ng+jj,       kk);
    *ivghsrc = ivoxel(   nx+ii,    ng+jj,    nz+kk);
    ++ivghost;
    ++ivghsrc;
    ++count;

    // (1,1,0) <- (0,0,1)
    //printf("(1,1,0) corner %d,%d,%d <- %d,%d,%d\n",
    //                        nx+ng+ii, ny+ng+jj,       kk,
    //                           ng+ii,    ng+jj,    nz+kk
    //);
    *ivghost = ivoxel(nx+ng+ii, ny+ng+jj,       kk);
    *ivghsrc = ivoxel(   ng+ii,    ng+jj,    nz+kk);
    ++ivghost;
    ++ivghsrc;
    ++count;

    // (0,0,1) <- (1,1,0)
    //printf("(0,0,1) corner %d,%d,%d <- %d,%d,%d\n",
    //                              ii,       jj, nz+ng+kk,
    //                           nx+ii,    ny+jj,    ng+kk
    //);
    *ivghost = ivoxel(      ii,       jj, nz+ng+kk);
    *ivghsrc = ivoxel(   nx+ii,    ny+jj,    ng+kk);
    ++ivghost;
    ++ivghsrc;
    ++count;

    // (1,0,1) <- (0,1,0)
    //printf("(1,0,1) corner %d,%d,%d <- %d,%d,%d\n",
    //                        nx+ng+ii,       jj, nz+ng+kk,
    //                           ng+ii,    ny+jj,    ng+kk
    //);
    *ivghost = ivoxel(nx+ng+ii,       jj, nz+ng+kk);
    *ivghsrc = ivoxel(   ng+ii,    ny+jj,    ng+kk);
    ++ivghost;
    ++ivghsrc;
    ++count;

    // (0,1,1) <- (1,0,0)
    //printf("(0,1,1) corner %d,%d,%d <- %d,%d,%d\n",
    //                              ii, ny+ng+jj, nz+ng+kk,
    //                           nx+ii,    ng+jj,    ng+kk
    //);
    *ivghost = ivoxel(      ii, ny+ng+jj, nz+ng+kk);
    *ivghsrc = ivoxel(   nx+ii,    ng+jj,    ng+kk);
    ++ivghost;
    ++ivghsrc;
    ++count;

    // (1,1,1) <- (0,0,0)
    //printf("(1,1,1) corner %d,%d,%d <- %d,%d,%d\n",
    //                        nx+ng+ii, ny+ng+jj, nz+ng+kk,
    //                           ng+ii,    ng+jj,    ng+kk
    //);
    *ivghost = ivoxel(nx+ng+ii, ny+ng+jj, nz+ng+kk);
    *ivghsrc = ivoxel(   ng+ii,    ng+jj,    ng+kk);
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
    // (x,0,0) <- (x,1,1)
    //if (ii==ng) printf("x-aligned edge (%d,%d) <- (%d,%d)\n",jj,kk,ny+jj,nz+kk);
    *ivghost = ivoxel(      ii,       jj,       kk);
    *ivghsrc = ivoxel(      ii,    ny+jj,    nz+kk);
    ++ivghost;
    ++ivghsrc;
    ++count;
    // (x,1,0) <- (x,0,1)
    //if (ii==ng) printf("x-aligned edge (%d,%d) <- (%d,%d)\n",ny+ng+jj,kk,ng+jj,nz+kk);
    *ivghost = ivoxel(      ii, ny+ng+jj,       kk);
    *ivghsrc = ivoxel(      ii,    ng+jj,    nz+kk);
    ++ivghost;
    ++ivghsrc;
    ++count;
    // (x,0,1) <- (x,1,0)
    //if (ii==ng) printf("x-aligned edge (%d,%d) <- (%d,%d)\n",jj,nz+ng+kk,ny+jj,ng+kk);
    *ivghost = ivoxel(      ii,       jj, nz+ng+kk);
    *ivghsrc = ivoxel(      ii,    ny+jj,    ng+kk);
    ++ivghost;
    ++ivghsrc;
    ++count;
    // (x,1,1) <- (x,0,0)
    //if (ii==ng) printf("x-aligned edge (%d,%d) <- (%d,%d)\n",ny+ng+jj,nz+ng+kk,ng+jj,ng+kk);
    *ivghost = ivoxel(      ii, ny+ng+jj, nz+ng+kk);
    *ivghsrc = ivoxel(      ii,    ng+jj,    ng+kk);
    ++ivghost;
    ++ivghsrc;
    ++count;
  }}}

  // y-aligned edges
  for (int kk=(    0); kk < (   ng); ++kk) {
  for (int jj=(   ng); jj < (ny+ng); ++jj) {
  for (int ii=(    0); ii < (   ng); ++ii) {
    // (0,y,0) <- (1,y,1)
    *ivghost = ivoxel(      ii,       jj,       kk);
    *ivghsrc = ivoxel(   nx+ii,       jj,    nz+kk);
    ++ivghost;
    ++ivghsrc;
    ++count;
    // (1,y,0) <- (0,y,1)
    *ivghost = ivoxel(nx+ng+ii,       jj,       kk);
    *ivghsrc = ivoxel(   ng+ii,       jj,    nz+kk);
    ++ivghost;
    ++ivghsrc;
    ++count;
    // (0,y,1) <- (1,y,0)
    *ivghost = ivoxel(      ii,       jj, nz+ng+kk);
    *ivghsrc = ivoxel(   nx+ii,       jj,    ng+kk);
    ++ivghost;
    ++ivghsrc;
    ++count;
    // (1,y,1) <- (0,y,0)
    *ivghost = ivoxel(nx+ng+ii,       jj, nz+ng+kk);
    *ivghsrc = ivoxel(   ng+ii,       jj,    ng+kk);
    ++ivghost;
    ++ivghsrc;
    ++count;
  }}}

  // z-aligned edges
  for (int kk=(   ng); kk < (nz+ng); ++kk) {
  for (int jj=(    0); jj < (   ng); ++jj) {
  for (int ii=(    0); ii < (   ng); ++ii) {
    // (0,0,z) <- (1,1,z)
    *ivghost = ivoxel(      ii,       jj,       kk);
    *ivghsrc = ivoxel(   nx+ii,    ny+jj,       kk);
    ++ivghost;
    ++ivghsrc;
    ++count;
    // (1,0,z) <- (0,1,z)
    *ivghost = ivoxel(nx+ng+ii,       jj,       kk);
    *ivghsrc = ivoxel(   ng+ii,    ny+jj,       kk);
    ++ivghost;
    ++ivghsrc;
    ++count;
    // (0,1,z) <- (1,0,z)
    *ivghost = ivoxel(      ii, ny+ng+jj,       kk);
    *ivghsrc = ivoxel(   nx+ii,    ng+jj,       kk);
    ++ivghost;
    ++ivghsrc;
    ++count;
    // (1,1,z) <- (0,0,z)
    *ivghost = ivoxel(nx+ng+ii, ny+ng+jj,       kk);
    *ivghsrc = ivoxel(   ng+ii,    ng+jj,       kk);
    ++ivghost;
    ++ivghsrc;
    ++count;
  }}}

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // all done counting cost, do some sanity checks

  printf("Field initialization: testing ghost count %d expected %d ghosts.\n", count, nvg);
  assert(count == nvg);

  printf("Field initialization: testing ghost linear voxel indices unique.\n");
  // Make a sorted copy and looking for adjacent duplicates.
  // DO NOT USE THE SORTED LIST IN PRODUCTION CODE, because the relative order
  // of ivoxels_ghost, ivoxels_ghsrc must be maintained.  --ATr,2024feb10
  int* ivoxels_test = (int*) malloc(nvg*sizeof(int));
  for (int ii =0; ii < count; ++ii) {
    ivoxels_test[ii] = ivoxels_ghost[ii];
  }
  qsort(ivoxels_test, count, sizeof(int), compare_ints);
  //printf("linear voxel unsorted %d sorted %d\n", ivoxels_ghost[0], ivoxels_test[0]);
  for (int ii = 1; ii < count; ++ii) {
    //printf("linear voxel unsorted %d sorted %d\n", ivoxels_ghost[ii], ivoxels_test[ii]);
    assert(ivoxels_test[ii-1] != ivoxels_test[ii]);
  }
  free(ivoxels_test);

  printf("Field initialization: all tests passed.\n");
  //assert(0);  // TODO while writing code

}
