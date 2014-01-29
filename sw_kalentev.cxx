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
  std::vector<boost::array<int, 4> > neighbors(lattice.num_sites());
  std::vector<boost::array<int, 4> > neighbor_bonds(lattice.num_sites());
  {
    std::vector<int> num_neighbors(lattice.num_sites(), 0);
    for (int b = 0; b < lattice.num_bonds(); ++b) {
      int s = lattice.source(b);
      int t = lattice.target(b);
      neighbors[s][num_neighbors[s]] = t;
      neighbor_bonds[s][num_neighbors[s]++] = b;
      neighbors[t][num_neighbors[t]] = s;
      neighbor_bonds[t][num_neighbors[t]++] = b;
    }
  }

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
  std::vector<int> active(lattice.num_bonds());
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
    for (int b = 0; b < lattice.num_bonds(); ++b) rn[b] = random01();
    for (int b = 0; b < lattice.num_bonds(); ++b)
      active[b] =
        ((spins[lattice.source(b)] == spins[lattice.target(b)]) && (rn[b] < prob)) ? 1 : 0;
    for (int s = 0; s < lattice.num_sites(); ++s) label[s] = s;

    while (true) {
      // scan
      bool converge = true;
      for (int s = 0; s < lattice.num_sites(); ++s) {
        for (int i = 0; i < 4; ++i) {
          if (active[neighbor_bonds[s][i]]) {
            converge &= (label[s] == label[neighbors[s][i]]);
            label[s] = std::min(label[s], label[neighbors[s][i]]);
          }
        }
      }
      if (converge) break;

      // analysis
      while (true) {
        bool changed = false;
        for (int s = 0; s < lattice.num_sites(); ++s) {
          changed |= (label[label[s]] != label[s]);
          label[label[s]] = label[s];
        }
        if (!changed) break;
      }

      for (int b = 0; b < lattice.num_bonds(); ++b) {
      }
      if (converge) break;
    }

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
