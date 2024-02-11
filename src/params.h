#ifndef PARAMS_H
#define PARAMS_H

typedef struct param {

  int idump;
  int isort;
  int ilast;
  int seed;

  float Lx;
  float Ly;
  float Lz;

  int nx;
  int ny;
  int nz;

  int npmax;
  int nppc;

} param_t;

#endif  //PARAMS_H
