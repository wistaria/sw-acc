#include <cuda.h>
#include <curand.h>

#include <stdlib.h>
#include <vector>

class CURNG{
public:
  CURNG(unsigned int seed = 12345U)
  {
    init(seed);
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
  void init(unsigned int seed)
  {
    curandCreateGenerator(&gen_,CURAND_RNG_PSEUDO_DEFAULT);
    set_seed(seed);
  }
  curandGenerator_t gen_;
};

template <typename T>
void fill_rand( std::vector<T> &vec, CURNG &gen)
{
  T *devData;
  const size_t N = vec.size();
  const size_t vsize = N*sizeof(T);
  cudaMalloc((void **)&devData, vsize);
  curandGenerateUniform(gen.generator(), devData, N);
  cudaMemcpy(&(vec[0]), devData, vsize, cudaMemcpyDeviceToHost);
  cudaFree(devData);
}

