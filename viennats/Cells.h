#ifndef DEF_CELLS_H
#define DEF_CELLS_H

/* =========================================================================
   Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.

                            -----------------
                 ViennaTS - The Vienna Topography Simulator
                            -----------------

   Contact:         viennats@iue.tuwien.ac.at

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


#include <vector>

namespace geom {

  template <int D> class cell {
  public:
    unsigned int Points[1 << D];
  };

  template <int D> class cells:public std::vector<cell<D> >  {};

}


#endif //DEF_CELLS_H
