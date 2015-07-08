#ifndef KCLUSTER_H
#define KCLUSTER_H

#include <iostream>
#include <vector>

#include <gsl/gsl_combination.h>

#include "kSphere.h"
#include "kMath.h"

class kCluster {

private:
  double r = 1.0;   // the radii of the spheres
  double r0 = 10.0; // the radii of the bourndary. Use a setter (set_boundary) to adjust.

public:
  std::vector<kSphere> Sphere;

  void setBoundary(double _r0) { r0 = _r0; }
  double getBoundary() { return r0; }

  double mesum() {
    double tmp;
    double result = 0.0;
    // sum of the sqaured distances between all spheres

    for (std::size_t j = 0; j < Sphere.size(); j++) {
      for (std::size_t i = 0; i < j; i++) {

        tmp = kMath::v3abs(kMath::v3diff(Sphere[i].pos, Sphere[j].pos));

        result += (tmp * tmp);
      }
    }
    return result;
  }

  // inverts the position of the spheres in the cluster which have there indices in the list of
  // indices.
  int invert(std::vector<size_t> Indice) {
    for (size_t i = 0; i < Indice.size(); i++) {
      Sphere[Indice[i]].pos.x = -Sphere[Indice[i]].pos.x;
      Sphere[Indice[i]].pos.y = -Sphere[Indice[i]].pos.y;
      Sphere[Indice[i]].pos.z = -Sphere[Indice[i]].pos.z;
    }
    return 0;
  }

  // returns a list of k indices (in a vector) which indicate the k spheres from the main
  // cluster
  // which represent the smallest potential energy among all posible groups of k spheres.
  std::vector<size_t> minU(size_t k, double r, double r0) {
    double pe0;
    double pe1 = std::numeric_limits<double>::max();
    std::vector<kSphere> S0;

    std::vector<size_t> tmp;
    std::vector<size_t> result;
    size_t n = Sphere.size();
    tmp.resize(k);
    result.resize(k);
    S0.resize(k);
    gsl_combination *c0 = gsl_combination_calloc(n, k);

    do {
      for (size_t i = 0; i < k; i++) {
        S0[i] = Sphere[c0->data[i]];
        tmp[i] = c0->data[i];
      }
      pe0 = pe(tmp, r, r0);
      // keep hold of the smallest potential on each iteration of the combinations.
      if (pe0 < pe1) {
        pe1 = pe0;
        result = tmp;
      }
    } while (gsl_combination_next(c0) == GSL_SUCCESS);

    gsl_combination_free(c0);
    return result;
  }

  // returns a list of k indices (in a vector) which indicate the k spheres from the main cluster
  // which represent the largest potential energy among all posible groups of k spheres.
  std::vector<size_t> maxU(size_t k, double r, double r0) {
    double pe0;
    double pe1 = std::numeric_limits<double>::max();
    std::vector<kSphere> S0;
    std::vector<size_t> tmp;
    std::vector<size_t> result;
    size_t n = Sphere.size();
    tmp.resize(k);
    result.resize(k);
    S0.resize(k);
    gsl_combination *c0 = gsl_combination_calloc(n, k);

    do {
      std::cout << "-> ";
      for (size_t i = 0; i < k; i++) {
        S0[i] = Sphere[c0->data[i]];
        tmp[i] = c0->data[i];
      }
      pe0 = pe(tmp, r, r0);
      // keep hold of the smallest potential on each iteration of the combinations.
      if (pe0 > pe1) {
        pe1 = pe0;
        result = tmp;
      }
    } while (gsl_combination_next(c0) == GSL_SUCCESS);

    gsl_combination_free(c0);
    return result;
  }

  // TODO: write comment!
  double d0_mag(kSphere a, double r, double r0) {
    double tmp = kMath::v3abs(a.pos) + r;
    if (tmp > r0) {
      return tmp - r0;
    } else {
      return 0.0;
    }
  }

  // The d0 vector is a vector pointing from its position to the origin so just the negative of the
  // postion with the magnitude found in d0_mag to use the d0_mag we will need to first unitize the
  // d0 "deformation" vector.
  v3 d0_vec(kSphere a, double r, double r0) {
    return kMath::v3scale(kMath::v3unit(a.pos), d0_mag(a, r, r0));
  }

  // TODO: write comment!
  static v3 d2_vec(kSphere a, kSphere b, double r) {
    v3 v;
    v = kMath::v3scale(kMath::v3unit(kMath::v3diff(a.pos, b.pos)), d2_mag(a, b, r));
    return v;
  }

  // TODO: write comment!
  static double d2_mag(kSphere a, kSphere b, double r) {
    double tmp = kMath::v3abs(kMath::v3diff(a.pos, b.pos)) / 2.0; // tmp will equal 0 if a == b.
    if (tmp < r && tmp != 0.0) { // second condition to ensure that for a = b, 0
                                 // is returned.
      return r - tmp;
    } else {
      return 0.0;
    }
  }

  // TODO: write comment!
  double pe(double r, double r0) {
    double sum0 = 0.0;
    double sum2 = 0.0;

    for (std::size_t i = 0; i < Sphere.size(); i++) {
      sum0 += d0_mag(Sphere[i], r, r0);
    }
    for (std::size_t j = 0; j < Sphere.size(); j++) {
      for (std::size_t i = 0; i < Sphere.size(); i++) {
        if (i != j) {
          sum2 += d0_mag(Sphere[i], r, r0);
        }
      }
    }
    return sum0 + sum2;
  }

  // Potential Energy function which takes indices for the Sphere list to evaluate over. In stead of
  // the entire list.
  double pe(std::vector<size_t> Index, double r, double r0) {
    double sum0 = 0.0;
    double sum2 = 0.0;

    for (std::size_t i = 0; i < Index.size(); i++) {
      sum0 += d0_mag(Sphere[Index[i]], r, r0);
    }

    for (std::size_t j = 0; j < Index.size(); j++) {
      for (std::size_t i = 0; i < Index.size(); i++) {
        if (i != j) {
          sum2 += d2_mag(Sphere[Index[i]], Sphere[Index[j]], r);
        }
      }
    }
    return sum0 + sum2;
  }

  // The gradient of the potential for a single sphere (against all others and the boundary).
  v3 gradpe(size_t i, double r, double r0) {
    v3 u;
    v3 v;

    u = d0_vec(Sphere[i], r, r0);

    for (std::size_t j = 0; j < Sphere.size(); j++) {
      // TODO: test to check if we need the following condition which should
      // probalby be covered by d2_vec.
      if (i != j) { // I dont think we need to test here if d2_vec is working
                    // correctly ... test
        v = d2_vec(Sphere[i], Sphere[j], r);
      } else {
        v.x = 0.0;
        v.y = 0.0;
        v.z = 0.0;
      }
      u.x += v.x;
      u.y += v.y;
      u.z += v.z;
    }
    return u;
  }

  // move the position to the next frame.
  void update(double r, double r0, double h) {
    v3 tmp;
    for (size_t i = 0; i < Sphere.size(); i++) {
      tmp = gradpe(i, r, r0);
      Sphere[i].pos.x -= h * tmp.x;
      Sphere[i].pos.y -= h * tmp.y;
      Sphere[i].pos.z -= h * tmp.z;
    }
    // std::cout << "size: " << tmp.x << std::endl;
  }

  // generic print function for the cluster.
  // TODO: override the << operator instead.
  void print() {
    for (size_t i = 0; i < Sphere.size(); i++) {
      std::cout << Sphere[i].pos.x << " " << Sphere[i].pos.y << " " << Sphere[i].pos.z << std::endl;
    }
    std::cout << std::endl;
  }

  // TODO: write comment!
  kCluster(size_t size, double _r) {
    r = _r; // I moved the radius from kSphere to here (kCluster) since we are using equal size
            // spheres. r0 is given an large default value but will need to be adjusted for
            // each configuration size and sphere radius.
    for (size_t i = 0; i < size; i++) {
      Sphere.push_back(kSphere());
    }
    //  minU(4, 0.5, 1.2);
  }

  // Copy Constructor.
  kCluster(const kCluster &obj) {
    // All we need to do (for the moment at least) is copy the Sphere vector.
    // Im hoping that this copies the values and not the addresses. Check!
    Sphere = obj.Sphere;
  }

  // TODO: write comment!
  void operator=(const kCluster &obj) { Sphere = obj.Sphere; }

  // TODO: write comment!
  ~kCluster() {}
};

#endif