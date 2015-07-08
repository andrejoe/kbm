//
//  kMath.h
//  kbm
//
//  Created by Erwin on 22/06/2015.
//  Copyright (c) 2015 Erwin. All rights reserved.
//

#ifndef kbm_kMath_h
#define kbm_kMath_h

#include <cmath>

typedef struct {
  double x;
  double y;
  double z;
} v3;

class kMath {
private:
public:
  static double v3abs(const v3 &v) { return sqrt(v.x * v.x + v.y * v.y + v.z * v.z); }

  // return a unit vector created from the vector argument.
  static v3 v3unit(const v3 &u) {
    v3 v;
    double mag = v3abs(u);
    v.x = u.x / mag;
    v.y = u.y / mag;
    v.z = u.z / mag;
    return v;
  }

  static v3 v3scale(const v3 &u, double scale) {
    v3 v;
    v.x = u.x * scale;
    v.y = u.y * scale;
    v.z = u.z * scale;
    return v;
  }

  static v3 v3diff(const v3 &u, const v3 &v) {
    v3 w;
    w.x = v.x - u.x;
    w.y = v.y - u.y;
    w.z = v.z - u.z;
    return w;
  }
};

#endif
