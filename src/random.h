#ifndef RANDOM_H
#define RANDOM_H

// https://en.wikipedia.org/wiki/Lehmer_random_number_generator
// See also std::minstd_rand
class Random {

  private:
    // Uniform RNG
    long long seed;
    const static long long m31   = 2147483647;  // Mersenne prime 2^31 - 1
    const static long long m31p1 = 2147483648;  // 2^31
    const static long long amult = 48271;

  public:
    Random(int seed_);

    long long rand();
    double drand();
    double normal(float mu, float sigma);

};

#endif  // RANDOM_H
