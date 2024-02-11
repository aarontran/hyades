#ifndef FIELDS_H
#define FIELDS_H

typedef struct field {
  // If you change field layout, also update
  // FieldArray.nfstruct, FieldArray.iget
  float ex,     ey,     ez,     epsx;
  float bx,     by,     bz,     epsy;
  float bx0,    by0,    bz0,    epsz;
  float jfx,    jfy,    jfz,    rhof;
  //float jfxold, jfyold, jfzold, rhofold;
  //float tempx,  tempy,  tempz,  tmpsm;
  //float smex,   smey,   smez,   pexx;
  //float smcbx,  smcby,  smcbz,  cmat;
} field_t;

class FieldArray {

  private:
    const static int nfstruct = 16;  // number of floats in field_t typedef

  public:
    FieldArray(int nx_, int ny_, int nz_, int ng_, float hx_, float hy_, float hz_, float dt_);

    int nx;  // x-axis cell count (no ghosts)
    int ny;  // y-axis ...
    int nz;  // z-axis ...
    int ng;  // ghost cell count, per coordinate axis, per side
    int nvall;  // total number of cells INCLUDING ghosts
    int nvg;  // total number of ghost cells (corners and slabs)
    float hx;  // cell size
    float hy;
    float hz;
    float dt;  // simulation timestep
    // the big kahuna
    field_t* f0;
    // ghost cell management
    // we might want to make (nx,ny,nz,ng) compile time constants so these
    // arrays can be pre-allocated and static... but that's a lot of work
    int* ivoxels_ghost;  // linear voxel indices for all ghost cells
    int* ivoxels_ghsrc;  // linear voxel indices for all ghosts' paired "source" cells

    // Low-level methods to iterate over and get/set field values
    field_t*      voxel(int ii, int jj, int kk);
    int          ivoxel(int ii, int jj, int kk);
    float*         fget(field_t* ff, int ii);
    void          fcopy(field_t* ff, field_t* gg);
    void         freset(float v);

    // High-level methods to set field values
    void      uniform_b(float bx, float by, float bz);
    void      uniform_e(float ex, float ey, float ez);
    void      ghost_sync_slab(); // TODO update my code...
    //void      ghost_sync_corner(); // TODO
    //void      ghost_cur_reset(); // TODO
    //void      ghost_cur_reduce(); // TODO
    // I need to be able to
    // * copy ALL source -> ghost, ALL fields
    // * add (aka reduce) ghost -> source, SELECTED fields
    // * zero ghost, SELECTED fields
    // * zero ghost, ALL fields
    //
    // TODO want a method to invert linear index into 3-tuple

};

#endif  // FIELDS_H
