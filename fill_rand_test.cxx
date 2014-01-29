#include <iostream>
#include "fill_rand.h"
#include <cstdlib>
#include <ctime>

int main(int argc, char **argv)
{
  using namespace std;
  const int N = (argc > 1)?std::atoi(argv[1]):8192;
  cout << "N = " << N << endl;


  const int seed = static_cast<unsigned int>(std::time(0));
  CURNG rng(seed);

  vector<float> vec(2*N,0);
  fill_rand(vec,rng);
  /*
  for(int i=0; i<vec.size(); ++i){
    cout << vec[i] << endl;
  }
  */

  int inside = 0;
  for(int i=0; i<N; ++i){
    float x = vec[2*i];
    float y = vec[2*i+1];
    float r2 = x*x + y*y;
    if(r2 < 1.0f){
      ++inside;
    }
  }
  cout << (4.0*inside)/N << endl;
}
