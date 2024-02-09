#ifndef FIELDS_H
#define FIELDS_H

typedef struct field {
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

  public:
    FieldArray(int nx_, int ny_, int nz_, int ng_, float hx_, float hy_, float hz_, float dt_);

    int nx;  // x-axis cell count (no ghosts)
    int ny;  // y-axis ...
    int nz;  // z-axis ...
    int ng;  // ghost cell count, per coordinate axis, per side
    int nv;  // total number of cells INCLUDING ghosts
    float hx;  // cell size
    float hy;
    float hz;
    float dt;  // simulation timestep
    field_t* f;

    field_t*      voxel(int ii, int jj, int kk);
    int          ivoxel(int ii, int jj, int kk);
    void      uniform_b(float bx, float by, float bz);
    void      uniform_e(float ex, float ey, float ez);

};

#endif  // FIELDS_H
