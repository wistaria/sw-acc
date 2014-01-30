/*
 * pgc++ cannot generate GPGPU code :(
 */

#include <iostream>
#include <cstdlib>
#include <ctime>

#include <cuda_runtime.h>
#include <curand_kernel.h>

int main(int argc, char **argv)
{
  const int N = (argc > 1) ? std::atoi(argv[1]) : 8192;
  const int M = (argc > 2) ? std::atoi(argv[2]) : 1024;

  const unsigned int seed = static_cast<unsigned int>(std::time(0));

  int inside=0;
#pragma acc parallel loop \
  copyin(seed) \
  reduction(+:inside)
  for(int i=0; i<N; ++i){
    curandState *state;
    curand_init(seed, i, 0, state);
    for(int j=0; j<M; ++j){
      float x = curand_uniform(state);
      float y = curand_uniform(state);
      if(x*x + y*y < 1.0f){
        ++inside;
      }
    }
  }

  std::cout << (4.0*inside)/(N*M) << std::endl;
}
