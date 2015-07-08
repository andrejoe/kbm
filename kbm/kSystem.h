#ifndef KSYSTEM_H
#define KSYSTEM_H

#include <vector>
#include "kCluster.h"

// kSystem holds a single kCluster. It is to define all system wide parameters.
// For example external forces such as gravity or boundary containers which
// affect the kCluster
// as a whole.

class kSystem {
private:
public:
  double r = 0.5;
  //  double r0 = 1.1123630599639879; // for 4 spheres.
  double r0 = 1.1; // for
  double h = 1.0;

  double min = 0.0001 * r0;
  double hmin = 0.01;
  double relpemin = 0.002;

  int cnt = 0;

  std::vector<kCluster> Cluster;

  // This is what is described as the A0 algorithm for the symmetric relocation
  // schemea

  // When the correct value for r0 is known the solutions are found very fast
  // and very acurately.
  // TODO: works well but can be tightened up alot and adjusted to better work
  // with spheres (it was designed for 2d).
  void srl_A0() {
    while ((Cluster[0].pe(r, r0) > min) || (h >= hmin) ||
           relpe(Cluster[0], Cluster[1]) >= relpemin) {
      Cluster[1].update(r, r0, h);
      if (Cluster[1].pe(r, r0) > Cluster[0].pe(r, r0)) {
        h *= 0.8;
      }
      Cluster[0] = Cluster[1];
      cnt++;

      std::cout << cnt << std::endl;
      std::cout << h << " " << Cluster[0].pe(r, r0) << " " << Cluster[1].pe(r, r0) << " "
                << relpe(Cluster[0], Cluster[1]) << std::endl;
      Cluster[1].print();
    }
  }

  // relative potential measures the relative difference between the potentials
  // of two Clusters.
  // The primary purpose is to identify local minima when searching for
  // solutions using steepest decent.
  // (eq 8) Wenqi and Yan 2004.

  double relpe(kCluster a, kCluster b) {
    double tmp1 = a.pe(r, r0);
    double tmp2 = b.pe(r, r0);

    if (tmp2 > tmp1) {
      return (tmp2 - tmp1) / tmp1;
    } else
      return (tmp1 - tmp2) / tmp1;
  }

  kSystem() {
    // We only need two clusters. One for the present and one for the next time
    // frame

    Cluster.push_back(kCluster(5, 0.5));
    // Clone the initial. Giving two identical clusters.

    // Cluster.push_back(kCluster(Cluster[0]));

    /*
    Cluster[1].update(r, r0, h);
    Cluster[1].update(r, r0, h);
    Cluster[1].update(r, r0, h);

    Cluster[0].print();
    Cluster[1].print();

    std::cout << "Made it here ... where is the rest?" << std::endl;

    std::cout << "mesum0: " << Cluster[0].mesum() << std::endl;
    std::cout << "mesum1: " << Cluster[1].mesum() << std::endl;
     */

    /*
      srl_A0();

    std::cout << cnt << std::endl;
    std::cout << h << " " << Cluster[0].pe(r, r0) << " " << Cluster[1].pe(r, r0)
              << " " << relpe(Cluster[0], Cluster[1]) << std::endl;
    Cluster[1].print();
      */
  }

  ~kSystem() {}
};

#endif
