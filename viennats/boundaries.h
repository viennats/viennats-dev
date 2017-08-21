#ifndef BOUNDARIES_H_
#define BOUNDARIES_H_

/* =========================================================================
   Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.

                            -----------------
                 ViennaTS - The Vienna Topography Simulator
                            -----------------

   Contact:         viennats@iue.tuwien.ac.at

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


#include <cassert>
#include <cmath>

///Namespace for boundary condition related objects
namespace bnc {

    //#############################################
    //Boundary conditions definitions
    //#############################################

    enum boundary_condition_type {
        REFLECTIVE_BOUNDARY,
        INFINITE_BOUNDARY,
        PERIODIC_BOUNDARY,
        EXTENDED_BOUNDARY
    };


    ///Stores and hangles boundary conditions of simulation domain
    class boundary_conditions_type {
    public:
        boundary_condition_type min;
        boundary_condition_type max;

        ///maps coordinates into simulation domain so they do not go outside
        template <class T, class T2>
        T2 map_coordinate(T minCoord, T maxCoord, T2 Coord, bool &ReverseSign) const {

            const T ext=maxCoord-minCoord;

            bool cycles=true;

            //assert(min!=INFINITE_BOUNDARY);         //TODO
            if (min!=EXTENDED_BOUNDARY) {
                while (Coord<minCoord) {
                    cycles=!cycles;
                    Coord+=ext;
                }
            }

            //assert(max!=INFINITE_BOUNDARY);         //TODO
            if (max!=EXTENDED_BOUNDARY) {
                while (Coord>maxCoord) {
                    cycles=!cycles;
                    Coord-=ext;
                }
            }

            if (cycles || (min==PERIODIC_BOUNDARY)) {
                ReverseSign=false;
                return Coord-minCoord;
            } else {
                ReverseSign=true;
                return maxCoord-Coord;
            }
        }

    };






}


#endif /* BOUNDARIES_H_ */
