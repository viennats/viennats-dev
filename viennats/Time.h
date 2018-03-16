#ifndef MYTIMER_H_
#define MYTIMER_H_

/* =========================================================================
   Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.

                            -----------------
                 ViennaTS - The Vienna Topography Simulator
                            -----------------

   Contact:         viennats@iue.tuwien.ac.at

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


#include <omp.h>

namespace my {
  ///Namespace for time or timing related methods.
  namespace time {

    double GetTime() {
      return omp_get_wtime();
    }

  }
}


#endif /*MYTIMER_H_*/
