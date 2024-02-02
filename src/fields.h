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

  private:
    int nx;
    int ny;
    int nz;
    int ng;
    int nv;

  public:
    FieldArray(int nx_, int ny_, int nz_, int ng_);

    field_t* f;

    field_t*      voxel(int ii, int jj, int kk);
    int          ivoxel(int ii, int jj, int kk);
    void      uniform_b(float bx, float by, float bz);
    void      uniform_e(float ex, float ey, float ez);

};

#endif  // FIELDS_H
