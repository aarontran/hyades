#include <stdlib.h> // need stdlib or cstdint to get int32_t at compile time???,
                    // don't understand why --ATr,2024feb09
#include <stdio.h>
//#include "field.h"
#include "interp.h"
#include "particle.h"

void ParticleArray::dump(int step) {
  printf("prtls would dump at step %d\n", step);
  return;
}
