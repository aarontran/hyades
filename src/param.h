#ifndef PARAM_H
#define PARAM_H

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
  int ng;

  int npmax;
  int nppc;

} param_t;

#endif  //PARAM_H
