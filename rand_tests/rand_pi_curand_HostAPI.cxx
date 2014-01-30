#include <iostream>
#include <ctime>
#include <cstdlib>
#include <boost/timer.hpp>
#include "../fill_rand.h"


int main(int argc, char **argv)
{
  const int N = (argc > 1)?std::atoi(argv[1]):8192;
  const size_t vsize = 2*N*sizeof(float);

  const int MCS = (argc > 2)?std::atoi(argv[2]):1000;

  const unsigned int seed = static_cast<unsigned int>(std::time(0));
  CURNG gen = makeCURNG_MT(seed);

  double pisum = 0.0;
  double pisum2 = 0.0;
  boost::timer tm;
  for(int mcs = 0; mcs < MCS; ++mcs){
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
    const double pi = 4.0*inside/N;
    pisum += pi;
    pisum2 += pi*pi;
    cudaFree(vec);
  }
  std::cout << "# N nproc pi time(ms)" << std::endl;
  std::cout << N << " " << " " << pisum/MCS << " " << 1000*tm.elapsed()/MCS  << std::endl;
}
