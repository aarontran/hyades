#include <assert.h>
#include <stdlib.h> // need stdlib or cstdint to get int32_t at compile time???,
                    // don't understand why --ATr,2024feb09
#include <stdio.h>

#include "hdf5.h"

//#include "field.h"
#include "interp.h"
#include "particle.h"

void ParticleArray::dump(int step) {

  char fname[100];
  snprintf(fname, 100, "output/prtl.%d.hdf5", step);

  herr_t status;

  hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  // float (32-bit) outputs
  hputf( file_id,   "x" );
  hputf( file_id,   "y" );
  hputf( file_id,   "z" );
  hputf( file_id,  "ux" );
  hputf( file_id,  "uy" );
  hputf( file_id,  "uz" );

  // int (32-bit) outputs
  hputi( file_id, "ind" );

  status = H5Fclose(file_id);
  assert(status >= 0);

  printf("Dumped %s\n", fname);
  return;
}

// hputf = put float-type particle attribute buffer into HDF5 file object
void ParticleArray::hputf(hid_t file_id, const char* attr_name) {

  float buf[np];
  particle_t* p = p0;
  for (int ii=0; ii<np; ++ii) {
    buf[ii] = *(pseek_fkey(attr_name, p));
    ++p;
  }

  hsize_t dims[1] = { (hsize_t)np };

  hid_t dspace_id = H5Screate_simple(1, dims, NULL);

  hid_t dset_id = H5Dcreate2(file_id, attr_name, H5T_NATIVE_FLOAT,
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

// hputi = put int32_t particle attribute buffer into HDF5 file object
void ParticleArray::hputi(hid_t file_id, const char* attr_name) {

  // Note: HDF5 "C9x" integer types let us specify bit width while using the
  // native endianness (in contrast to H5T_NATIVE_INT or H5T_STD_U32LE).
  // Documentation is at:
  // https://docs.hdfgroup.org/hdf5/develop/group___p_d_t_c9x.html

  int32_t buf[np];
  particle_t* p = p0;
  for (int ii=0; ii<np; ++ii) {
    buf[ii] = *(pseek_ikey(attr_name, p));
    ++p;
  }

  hsize_t dims[1] = { (hsize_t)np };

  hid_t dspace_id = H5Screate_simple(1, dims, NULL);

  hid_t dset_id = H5Dcreate2(file_id, attr_name, H5T_NATIVE_INT32,
                             dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  herr_t status;
  status = H5Dwrite(dset_id, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    &buf);
  assert(status >= 0);

  status = H5Dclose(dset_id);
  assert(status >= 0);

  status = H5Sclose(dspace_id);
  assert(status >= 0);
}
