#include "observable.h"
#include "options.h"
#include "square_lattice.h"

#include <boost/foreach.hpp>
#include <boost/random.hpp>
#include <boost/timer.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

int main(int argc, char* argv[]) {
  // parameters
  options p(argc, argv, 8, 2.27);
  if (!p.valid) std::exit(-1);
  const unsigned int sweeps = p.sweeps;
  const unsigned int therm = p.therm;
  const double beta = 1 / p.temperature;
  const double prob = 1 - std::exp(- 2 * beta);

  // sqaure lattice
  square_lattice lattice(p.length);

  // random number generators
  typedef boost::mt19937 engine_t;
  typedef boost::uniform_01<engine_t&> random01_t;
  engine_t engine(29833u);
  random01_t random01(engine);

  // spin configuration
  std::vector<int> spins(lattice.num_sites());
  std::fill(spins.begin(), spins.end(), 0 /* all up */);
  std::vector<double> rn(std::max(lattice.num_sites(), lattice.num_bonds()));

  // cluster information
  std::vector<int> label(lattice.num_sites());

  // oservables
  observable num_clusters;
  observable usus; // uniform susceptibility
  int nc;
  double mag2;

  //
  // Monte Carlo steps
  //

  boost::timer tm;

  for (unsigned int mcs = 0; mcs < therm + sweeps; ++mcs) {

    //
    // cluster construction
    //

    // initialize cluster information
    for (int s = 0; s < lattice.num_sites(); ++s) label[s] = s;

    // scan
    for (int b = 0; b < lattice.num_bonds(); ++b) rn[b] = random01();
    for (int b = 0; b < lattice.num_bonds(); ++b) {
      int s = lattice.source(b);
      int t = lattice.target(b);
      if (spins[s] == spins[t] && rn[b] < prob) label[s] = label[t] = std::min(s, t);
    }

    // analysis
    int rem;
    do {
      rem = 0;
      for (int s = 0; s < lattice.num_sites(); ++s) {
        if (label[label[s]] != label[s]) {
          label[s] = label[label[s]];
          ++rem;
        }
      }
    } while (rem > 0);

    //
    // cluster flip
    //

    // assign cluster id & accumulate cluster properties
    nc = 0;
    for (int s = 0; s < lattice.num_sites(); ++s) {
      if (label[s] == s) ++nc;
    }
    mag2 = 0;
    for (int s = 0; s < lattice.num_sites(); ++s) mag2 += 2 * spins[s] - 1;
    mag2 = mag2 * mag2;

    // determine whether clusters are flipped or not
    for (int s = 0; s < lattice.num_sites(); ++s) rn[s] = random01();

    // flip spins
    for (int s = 0; s < lattice.num_sites(); ++s) if (rn[label[s]] < 0.5) spins[s] ^= 1;

    //
    // measurements
    //

    if (mcs >= therm) {
      num_clusters << nc;
      usus << beta * mag2 / lattice.num_sites();
    }
  }

  double elapsed = tm.elapsed();
  std::clog << "Elapsed time = " << elapsed << " sec\n"
            << "Speed = " << (therm + sweeps) / elapsed << " MCS/sec\n";
  std::cout << "Number of Clusters        = "
            << num_clusters.mean() << " +- " << num_clusters.error() << std::endl
            << "Uniform Susceptibility    = "
            << usus.mean() << " +- " << usus.error() << std::endl;
}
