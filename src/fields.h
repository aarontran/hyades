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
    field_t* f0;

    // Low-level methods to iterate over and get/set field values
    field_t*      voxel(int ii, int jj, int kk);
    int          ivoxel(int ii, int jj, int kk);
    float*         fget(field_t* ff, int ii);
    void          fcopy(field_t* ff, field_t* gg);
    void         freset(float v);

    // High-level methods to set field values
    void      uniform_b(float bx, float by, float bz);
    void      uniform_e(float ex, float ey, float ez);
    void      update_ghost();

};

#endif  // FIELDS_H
