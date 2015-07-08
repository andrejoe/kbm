#include <cstdlib>
#include <cmath>
#include <ctime>
#include "kSphere.h"

#include "kMath.h"

// Let the default constructor simply create a unit sphere at the origin.
kSphere::kSphere() {
  pos.x = 0.0;
  pos.y = 0.0;
  pos.z = 0.0;
}

/*
   kSphere::kSphere() {
    initialise_rseed();
    create_random(4.0);
    kSphere::r = 1.0;
   }
 */

kSphere::kSphere(double _R) {
  initialise_rseed();
  create_random(_R);
}

kSphere::~kSphere() {}

int kSphere::initialise_rseed() {
  if (!rseeded) {
    srand48((long)std::clock());
    rseeded = true;
  }
  return 0;
}

double kSphere::create_random(double R) {
  double phi, theta, rho;

  phi = M_PI * drand48();
  theta = 2.0 * M_PI * drand48();

  rho = R * drand48();

  kSphere::pos.x = rho * cos(theta) * sin(phi);
  kSphere::pos.y = rho * sin(theta) * sin(phi);
  kSphere::pos.z = rho * cos(phi);

  return kMath::v3abs(pos);
}
