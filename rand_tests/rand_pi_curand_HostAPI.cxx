#include <iostream>
#include <ctime>
#include <cstdlib>
#include "fill_rand.h"


int main(int argc, char **argv)
{
  const int N = (argc > 1)?std::atoi(argv[1]):8192;
  const size_t vsize = 2*N*sizeof(float);

  const unsigned int seed = static_cast<unsigned int>(std::time(0));
  CURNG gen = makeCURNG_MT(seed);

  float *vec;
  cudaMalloc((void **)&vec, vsize);
  curandGenerateUniform(gen.generator(), vec, 2*N);

  int inside = 0;
#pragma acc parallel loop \
  deviceptr(vec) \
  reduction(+:inside)
  for(int i=0; i<N; ++i){
    float x = vec[2*i];
    float y = vec[2*i+1];
    if(x*x + y*y < 1.0f)
      ++inside;
  }

  std::cout << N << " " << (4.0*inside)/N << std::endl;
}
