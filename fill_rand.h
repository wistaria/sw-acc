#include <stdlib.h>
#include <vector>

#include <cuda_runtime.h>
#include <curand.h>

class CURNG{
public:
  CURNG(unsigned int seed = 12345U, curandRngType_t type=CURAND_RNG_PSEUDO_DEFAULT)
  {
    curandCreateGenerator(&gen_, type);
    set_seed(seed);
  }
  ~CURNG()
  {
    curandDestroyGenerator(gen_);
  }
  void set_seed(unsigned int seed)
  {
    curandSetPseudoRandomGeneratorSeed(gen_, seed);
  }
  curandGenerator_t generator(){return gen_;}
private:
  curandGenerator_t gen_;
};

inline 
CURNG makeCURNG_MT(unsigned int seed = 12345U)
{
  return CURNG(seed, CURAND_RNG_PSEUDO_MTGP32); 
}

void fill_rand( float *vec, size_t N, CURNG &gen)
{
  float *devData;
  const size_t vsize = N*sizeof(float);
  cudaMalloc((void **)&devData, vsize);
  curandGenerateUniform(gen.generator(), devData, N);
  cudaMemcpy(vec, devData, vsize, cudaMemcpyDeviceToHost);
  cudaFree(devData);
}

inline
void fill_rand( std::vector<float> &vec, CURNG &gen)
{
  const size_t N = vec.size();
  float *inner = &(vec[0]);
  fill_rand(inner, N, gen);
}

