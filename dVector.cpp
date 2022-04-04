#include "dVector.h"

void crossProduct(dVector& a, dVector& b, dVector& c)
{
  double x1 = a[0], x2 = b[0];
  double y1 = a[1], y2 = b[1];
  double z1 = a[2], z2 = b[2];

  c[0] = y1 * z2 - z1 * y2;
  c[1] = z1 * x2 - x1 * z2;
  c[2] = x1 * y2 - y1 * x2;
}

double dotProduct(dVector& a, dVector& b)
{
  return (a * b).sum();
}

double length(dVector& a)
{
  return sqrt((a * a).sum());
}
