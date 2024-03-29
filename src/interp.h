#ifndef INTERP_H
#define INTERP_H

#include "field.h"

typedef struct interp {
  // Second-order accurate interpolation coefficients
  // TODO are the mixed partials needed? Disabled for now. -ATr,2024feb09
  float ex, dexdx, dexdy, dexdz, d2exdx, d2exdy, d2exdz; //d2exdydz;
  float ey, deydx, deydy, deydz, d2eydx, d2eydy, d2eydz; //d2eydzdx;
  float ez, dezdx, dezdy, dezdz, d2ezdx, d2ezdy, d2ezdz; //d2ezdxdy;
  float bx, dbxdx, dbxdy, dbxdz, d2bxdx, d2bxdy, d2bxdz;
  float by, dbydx, dbydy, dbydz, d2bydx, d2bydy, d2bydz;
  float bz, dbzdx, dbzdy, dbzdz, d2bzdx, d2bzdy, d2bzdz;
  float _pad[2];                       // 16-byte align    
} interp_t;

class InterpArray {

  public:
    InterpArray(FieldArray fa_);

    FieldArray  fa;
    interp_t*  ic0;  // ic = interpolation coefficients

    interp_t* voxel(int ii, int jj, int kk);
    void      update();

    // Interpolate field(s) to (x,y,z) location given by voxel + offsets.
    // Voxel offsets dx,dy,dz are within interval [-1,1].
    float exloc(interp_t* ic, float dx, float dy, float dz);
    float eyloc(interp_t* ic, float dx, float dy, float dz);
    float ezloc(interp_t* ic, float dx, float dy, float dz);
    float bxloc(interp_t* ic, float dx, float dy, float dz);
    float byloc(interp_t* ic, float dx, float dy, float dz);
    float bzloc(interp_t* ic, float dx, float dy, float dz);

};

#endif  // INTERP_H
