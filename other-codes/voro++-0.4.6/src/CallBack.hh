#ifndef VORO_CALLBACK_HH
#define VORO_CALLBACK_HH

#include "cell.hh"
using namespace voro;

namespace voro {
class voronoiCallBack
{
public :
  virtual void invoke( int id, voronoicell_neighbor &cell) = 0;
  virtual void invoke( int id, voronoicell &cell) = 0;
};
}
#endif
