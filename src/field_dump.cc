#include <assert.h>
//#include <stdlib.h>  // for malloc, qsort
#include <stdio.h>
#include <string.h>

#include "hdf5.h"

#include "field.h"

// Dump useful field arrays to disk, live cells only (omit ghost cells).
// step is only used to set the filename.
void FieldArray::dump(int step, const char* formatstr, bool with_ghost) {

  char fname[100];
  snprintf(fname, 100, formatstr, step);

  herr_t status;

  hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  if (with_ghost) {
    hputf_gh( file_id,   "ex" );
    hputf_gh( file_id,   "ey" );
    hputf_gh( file_id,   "ez" );
    hputf_gh( file_id,   "bx" );
    hputf_gh( file_id,   "by" );
    hputf_gh( file_id,   "bz" );
    hputf_gh( file_id,  "jfx" );
    hputf_gh( file_id,  "jfy" );
    hputf_gh( file_id,  "jfz" );
    hputf_gh( file_id, "rhof" );
  } else {
    // "normal" case
    hputf( file_id,   "ex" );
    hputf( file_id,   "ey" );
    hputf( file_id,   "ez" );
    hputf( file_id,   "bx" );
    hputf( file_id,   "by" );
    hputf( file_id,   "bz" );
    hputf( file_id,  "jfx" );
    hputf( file_id,  "jfy" );
    hputf( file_id,  "jfz" );
    hputf( file_id, "rhof" );
  }

  status = H5Fclose(file_id);
  assert(status >= 0);

  printf("Dumped %s\n", fname);

  return;
}

// hputf = put float-type field buffer into HDF5 file object
void FieldArray::hputf(hid_t file_id, const char* field_name) {

  // Output buffer omits ghost cells
  float buf[nx][ny][nz];
  for (int kk = ng; kk < (nz+ng); ++kk) {
  for (int jj = ng; jj < (ny+ng); ++jj) {
  for (int ii = ng; ii < (nx+ng); ++ii) {
    field_t* fv = voxel(ii,jj,kk);
    buf[ii-ng][jj-ng][kk-ng] = *(fseek_key(field_name, fv));
  }}}

  hsize_t dims[3] = { (hsize_t)nx, (hsize_t)ny, (hsize_t)nz };

  hid_t dspace_id = H5Screate_simple(3, dims, NULL);

  hid_t dset_id = H5Dcreate2(file_id, field_name, H5T_NATIVE_FLOAT,
                             dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  herr_t status;
  status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    &buf);
  assert(status >= 0);

  status = H5Dclose(dset_id);
  assert(status >= 0);

  status = H5Sclose(dspace_id);
  assert(status >= 0);
}

// hputf_gh = put float-type field buffer with ghost cells into HDF5 file object
void FieldArray::hputf_gh(hid_t file_id, const char* field_name) {

  // Output buffer with ghost cells
  float buf[nx+2*ng][ny+2*ng][nz+2*ng];
  for (int kk = 0; kk < (nz+2*ng); ++kk) {
  for (int jj = 0; jj < (ny+2*ng); ++jj) {
  for (int ii = 0; ii < (nx+2*ng); ++ii) {
    field_t* fv = voxel(ii,jj,kk);
    buf[ii][jj][kk] = *(fseek_key(field_name, fv));
  }}}

  hsize_t dims[3] = { (hsize_t)(nx+2*ng), (hsize_t)(ny+2*ng), (hsize_t)(nz+2*ng) };

  hid_t dspace_id = H5Screate_simple(3, dims, NULL);

  hid_t dset_id = H5Dcreate2(file_id, field_name, H5T_NATIVE_FLOAT,
                             dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  herr_t status;
  status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    &buf);
  assert(status >= 0);

  status = H5Dclose(dset_id);
  assert(status >= 0);

  status = H5Sclose(dspace_id);
  assert(status >= 0);
}

// hputd = put double-type field buffer into HDF5 file object
//void FieldArray::hputd(hid_t file_id, const char* field_name) {
//}
