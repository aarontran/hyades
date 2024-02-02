#include <cmath>  // for log, cos, sin, pi
#include <stdio.h>  // for printf

#include "random.h"

// Constructor
Random::Random(int seed_) {
  if (seed_ <= 0) {
    seed = 1; // RNG scheme doesn't work for zero seed
  } else {
    seed = (long long)seed_;
  }
}

// Uniform random numbers on [0,1) except for a small interval near 1.
// https://en.wikipedia.org/wiki/Lehmer_random_number_generator
// See also std::minstd_rand
double Random::drand() {
  long long seed0 = seed;
  seed = (amult*seed) % m31;
  return seed0 / ((double)m31p1);
}

// Actual range is 2^31
long long Random::rand() {
  long long seed0 = seed;
  seed = (amult*seed) % m31;
  return seed0;
}

// Normal distribution using Box-Muller transform
// Not very fast for time-sensitive work
double Random::normal(float mu, float sigma) {
  float u1 = drand();
  float u2 = drand();
  // discard the unused sin(...) sample
  float x = sqrt(-2*log(u1))*cos(M_PI*u2/2);
  return (x-mu)*sigma;
}
