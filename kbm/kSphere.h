#ifndef kbm_kSphere_h
#define kbm_kSphere_h

#include "kMath.h"

class kSphere {
private:
  // switch to ensure we only seed the random number generator once.
  bool rseeded = false;

public:
  v3 pos;

  //  double r; // radius of the sphere // Since we will be using sphere so equal size, this is
  // moved to the cluster class.
  double R; // boundary radius in which all spheres are created.

  // double dst;

  kSphere();
  kSphere(double _R);

  ~kSphere();

  int initialise_rseed();
  double create_random(double radius);
};

#endif
