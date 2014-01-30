#include <iostream>
#include <cstdlib>
#include <ctime>
#include <boost/random.hpp>
#include <boost/timer.hpp>
#include <omp.h>

int main(int argc, char **argv)
{
  using namespace std;

  const int N = (argc > 1)?std::atoi(argv[1]):8192;

  const int MCS = (argc > 2)?std::atoi(argv[2]):1000;

  const int seed = static_cast<unsigned int>(std::time(0));
  boost::uniform_real<float> dist(0.f,1.0f);

  const int nth = omp_get_max_threads();

  std::vector<boost::mt19937> rngs(nth);
  for(int i=0; i<nth; ++i){
    rngs[i] = boost::mt19937(seed + 137*i);
  }

  double pisum = 0.0;
  double pisum2 = 0.0;
  //boost::timer tm;
  double tm1 = omp_get_wtime();
  for(int mcs = 0; mcs < MCS; ++mcs){

    int inside = 0;
#pragma omp parallel for reduction(+:inside)
    for(int i=0; i<N; ++i){
      int id = omp_get_thread_num();
      //cout << id << " " << endl;
      float x = dist(rngs[id]);
      float y = dist(rngs[id]);
      float r2 = x*x + y*y;
      if(r2 < 1.0f){
        ++inside;
      }
    }

    const double pi =  (4.0*inside)/N;
    pisum += pi;
    pisum2 += pi*pi;

  }
  //cout << "# N proc pi milisec" << endl;
  cout << N << " " << " " << nth  << " "  << pisum/MCS << " " << 1000*(omp_get_wtime() - tm1)/MCS << endl;
}
