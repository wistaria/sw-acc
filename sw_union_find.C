#include "observable.h"
#include "options.h"
#include "square_lattice.h"
#include "union_find.h"

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
  typedef looper::union_find::node fragment_t;
  std::vector<fragment_t> fragments(lattice.num_sites());
  std::vector<bool> to_flip(lattice.num_sites());

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
    for (int s = 0; s < lattice.num_sites(); ++s) fragments[s] = fragment_t();


    for (int b = 0; b < lattice.num_bonds(); ++b) rn[b] = random01();
    for (int b = 0; b < lattice.num_bonds(); ++b) {
      if (spins[lattice.source(b)] == spins[lattice.target(b)] && rn[b] < prob) {
        unify(fragments, lattice.source(b), lattice.target(b));
      }
    }

    //
    // cluster flip
    //

    // assign cluster id & accumulate cluster properties
    nc = 0;
    mag2 = 0;
    for (int f = 0; f < fragments.size(); ++f) {
      if (fragments[f].is_root()) {
        fragments[f].set_id(nc++);
        mag2 += fragments[f].weight() * fragments[f].weight();
      }
    }
    for (int f = 0; f < fragments.size(); ++f) fragments[f].set_id(cluster_id(fragments, f));

    // determine whether clusters are flipped or not
    for (int s = 0; s < lattice.num_sites(); ++s) rn[s] = random01();
    for (int c = 0; c < nc; ++c) to_flip[c] = false;
    for (int s = 0; s < lattice.num_sites(); ++s) {
      int c = cluster_id(fragments, s);
      to_flip[c] = to_flip[c] ^ (rn[s] < 0.5);
    }

    // flip spins
    for (int s = 0; s < lattice.num_sites(); ++s) if (to_flip[fragments[s].id()]) spins[s] ^= 1;

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
