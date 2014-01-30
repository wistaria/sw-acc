
#include <iostream>
#include <cstdlib>
#include <ctime>

#include <cuda_runtime.h>
#include <curand_kernel.h>

__device__ void calc_pi(int N, curandState *state, double *pi)
{
  int insided = 0;
  for(int i=0; i<N; ++i){
    float x = curand_uniform(state);
    float y = curand_uniform(state);
    insided += (x*x+y*y < 1.0f)? 1 : 0;
  }
  *pi = (4.0*insided)/N;
}

int main(int argc, char **argv)
{
  const int N = (argc > 1) ? std::atoi(argv[1]) : 8192;
  const int M = (argc > 2) ? std::atoi(argv[2]) : 1024;

  const unsigned int seed = static_cast<unsigned int>(std::time(0));

  double pi = 0.0;
#pragma acc parallel loop \
  copyin(seed) \
  reduction(+:pi)
  for(int i=0; i<N; ++i){
    curandState *state;
    curand_init(seed, i, 0, state);
    calc_pi(M, state, &pi);
  }

  std::cout << pi/N << std::endl;
}
