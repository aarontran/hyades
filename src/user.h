#ifndef USER_H
#define USER_H

#include "field.h"
#include "param.h"
#include "particle.h"

void user_param(param_t* par);

void user_initialize(FieldArray* fa, ParticleArray* ions, param_t* par);

#endif  // USER_H
