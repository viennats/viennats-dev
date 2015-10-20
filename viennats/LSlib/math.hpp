#ifndef MATH_HPP_
#define MATH_HPP_

/* =========================================================================
   Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.

                            -----------------
                 ViennaTS - The Vienna Topography Simulator
                            -----------------

   Contact:         viennats@iue.tuwien.ac.at

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


#include "config.hpp"
#include <cmath>

#ifndef LVLSET_HAVE_COPYSIGN
    #ifdef LVLSET_HAVE__COPYSIGN
        #define copysign _copysign
        #define LVLSET_HAVE_COPYSIGN
    #endif
#endif

namespace lvlset {

    namespace math {

        template<class T> inline int sign(T x) {
            if (x>T(0)) return 1;
            if (x<T(0)) return -1;
            return 0;
        }

        template<class T> inline T pow2(T x) {
            return x*x;
        }

        template<class T> inline bool signbit(T x) {
            #ifdef LVLSET_HAVE_STD_SIGNBIT
                return std::signbit(x);
            #else
                return copysign(static_cast<T>(1),x)<static_cast<T>(0);
            #endif
        }

        template< class T> inline T abs(T x) {
            if (x>=T(0)) return x; else return -x;
        }

        inline int abs(int x) {
            return std::abs(x);
        }

        inline long int abs(long int x) {
            return std::abs(x);
        }

        inline double abs(double x) {
            return std::fabs(x);
        }

        inline long double abs(long double x) {
            return std::fabs(x);
        }

        inline float abs(float x) {
            return std::fabs(x);
        }


    }
}
#endif /*MATH_HPP_*/
