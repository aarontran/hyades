#ifndef FIELD_H
#define FIELD_H

#include "hdf5.h"

typedef struct field {
  // If you change struct layout, also update FieldArray member data and
  // subroutines: nfstruct, dump(...), fseek_one(...), fseek_key(...), and any
  // other code that accesses struct members by explicit index or name
  float ex,     ey,     ez,     tmpx;
  float bx,     by,     bz,     tmpy;
  float bx0,    by0,    bz0,    tmpz;
  float jfx,    jfy,    jfz,    rhof;
  float jfx0,   jfy0,   jfz0,   rhof0;
  float smex,   smey,   smez,   pad;//pexx;
  //float smcbx,  smcby,  smcbz,  cmat;
} field_t;

class FieldArray {

  private:
    const static int nfstruct = 24;  // number of floats in field_t typedef

  public:
    FieldArray(int nx_, int ny_, int nz_, int ng_, float hx_, float hy_, float hz_, float dt_);

    // --------------------------------------------------
    // Data

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

    // how many subcycles per simulation timestep
    int nsubcycle_;

    // Hybrid algorithm parameters
    float hyb_te_ref_;
    float hyb_ne_ref_;
    float hyb_ne_floor_;
    float hyb_eta_;
    float hyb_hypereta_;

    // the big kahuna
    field_t* f0;
    // ghost cell management
    // we might want to make (nx,ny,nz,ng) compile time constants so these
    // arrays can be pre-allocated and static... but that's a lot of work
    int* ivoxels_ghost;  // linear voxel indices for all ghost cells
    int* ivoxels_ghsrc;  // linear voxel indices for all ghosts' paired "source" cells

    // --------------------------------------------------
    // Low-level methods to iterate over and get/set field values at single
    // voxels, not distinguishing live/ghost cells.
    // Argument ordering: we get/reset/copy the field member mm from voxel(s) ff, gg.
    // Dummy vars: ii,jj,kk for grid indices; mm for field struct members.
    field_t*      voxel(int ii, int jj, int kk);
    int          ivoxel(int ii, int jj, int kk);
    float*    fseek_key(const char* name,field_t* ff);
    float*    fseek_one(int mm,          field_t* ff);
    void       fset_one(int mm, float v, field_t* ff);
    void       fset_all(        float v, field_t* ff);
    void      fcopy_one(int mm, field_t* ff, field_t* gg);
    void      fcopy_all(        field_t* ff, field_t* gg);
    // TODO want a method to invert linear index into 3-tuple

    // --------------------------------------------------
    // Intermediate-level methods to get/set field values on live/ghost mesh
    // blocks.
    // Ghost methods operate on ALL of face, edge, and corner blocks.
    // Notes on naming?
    // * Synonyms for "reset" include: zero, null, reset, clear, empty, flush.
    // * Synonyms for "copy/reduce" include: fetch, sync, update, clone, yank,
    //   pull/push, gather/scatter/reduce, get/put, fill, mirror, replicate.
    // * I'd like a verb that connotes a one-way data movement for ghosts.
    // --ATr,2024feb10-11
    //void mesh_set_one     (int mm, float v);  // not sure if needed yet
    void mesh_set_all     (        float v);
    //void alive_set_one    (int mm, float v);  // not sure if needed yet
    //void alive_set_all    (        float v);  // not sure if needed yet
    //void ghost_set_one    (int mm, float v);  // not implemented
    //void ghost_set_all    (        float v);  // not sure if needed yet
    //void ghost_copy_one   (int mm);  // not sure if neeeded yet
    void ghost_copy_all   ();
    //void ghost_reduce_one (int mm);  // not sure if neeeded yet

    // --------------------------------------------------
    // High-level methods to set field values on the mesh
    void mesh_set_b   (float bx, float by, float bz);
    void mesh_set_e   (float ex, float ey, float ez);
    void mesh_set_jrho(float jx, float jy, float jz, float rho);
    void mesh_set_jrho0();

    //void ghost_set_eb      (float v);  // not implemented
    //void ghost_set_jrho    (float v);
    void ghost_copy_eb     ();
    void ghost_copy_e      ();
    void ghost_copy_b      ();
    void ghost_copy_jrho   ();
    void ghost_deposit_jrho();

    // --------------------------------------------------
    // High-level methods for top-level hybrid algorithm
    void advance_eb_rk4_ctrmesh   ();
    void advance_eb_rk4_ctrmesh   (int isub, int nsub);
    void advance_e_ctrmesh        (float frac);

    // High-level dump methods
    void dump    (int step, const char* formatstr, bool with_ghost);
    void hputf   (hid_t file_id, const char* field_name);
    void hputf_gh(hid_t file_id, const char* field_name);
    //void hputd(hid_t file_id, const char* field_name);

};

#endif  // FIELD_H
