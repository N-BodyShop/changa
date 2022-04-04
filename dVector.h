#ifndef DVECTOR_H
#define DVECTOR_H

#ifndef TESTCODE
#include <pup_stl.h>
#endif


#define NDIM 3

class dVector {
  double data[NDIM];

public:

  dVector() {}

  dVector( const dVector &a)
  {
    for( int i = 0; i < NDIM; i++)
      data[i] = a.data[i];
  }

  dVector( int size) {}

  dVector(const double &x, const double &y, const double &z) {
    data[0] = x; data[1] = y; data[2] = z;
  }

  inline int size() {
    return NDIM;
  }

  inline dVector& operator=(const dVector& v) {
    for (int i = 0; i < NDIM; i++) data[i] = v.data[i];

    return *this;
  }

  inline dVector& operator=(const double& d) {
    for (int i = 0; i < NDIM; i++) data[i] = d;

    return *this;
  }

  inline double sum() {
    double sum = 0.;

    for (int i = 0; i < NDIM; i++) sum += data[i];
    return sum;
  }

  inline double& operator[](const int& i) {
    return data[i];
  }

  inline dVector operator*(const dVector& v) {
    dVector a;

    for (int i = 0; i < NDIM; i++) a.data[i] = data[i] * v.data[i];
    return a;
  }

  inline dVector operator*(const double d) {
    dVector a;

    for (int i = 0; i < NDIM; i++)
      a.data[i] = data[i] * d;
    return a;
  }

  inline dVector operator/(const double d) {
    dVector a;

    for (int i = 0; i < NDIM; i++) a.data[i] = data[i] / d;
    return a;
  }

  inline dVector& operator/=(const double d) {
    for (int i = 0; i < NDIM; i++) data[i] /= d;
    return *this;
  }

  inline dVector& operator*=(const double d) {
    for (int i = 0; i < NDIM; i++) data[i] *= d;
    return *this;
  }

  inline dVector operator-(const dVector& v) {
    dVector a;

    for (int i = 0; i < NDIM; i++) a.data[i] = data[i] - v.data[i];
    return a;
  }

  inline dVector operator+(const dVector& v) {
    dVector a;

    for (int i = 0; i < NDIM; i++) a.data[i] = data[i] + v.data[i];
    return a;
  }

  inline void     resize(int size) {}

  inline dVector& operator+=(const dVector& v) {
    for (int i = 0; i < NDIM; i++) data[i] += v.data[i];

    return *this;
  }

  inline dVector& operator-=(const dVector& v) {
    for (int i = 0; i < NDIM; i++) data[i] -= v.data[i];

    return *this;
  }

#ifndef TESTCODE
  void pup(PUP::er& p)
  {
    for (int i = 0; i < NDIM; i++) p | data[i];
  }
#endif
};

inline dVector operator*( const double d, dVector& a) {
  //double x = a[0];
  return dVector( a[0]*d, a[1]*d, a[2]*d);
}

void   crossProduct(dVector& a,
                    dVector& b,
                    dVector& c);
double dotProduct(dVector& a,
                  dVector& b);
double length(dVector& a);
#endif
