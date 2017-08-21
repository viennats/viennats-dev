#ifndef BOX_HPP_
#define BOX_HPP_

/* =========================================================================
   Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.

                            -----------------
                 ViennaTS - The Vienna Topography Simulator
                            -----------------

   Contact:         viennats@iue.tuwien.ac.at

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


#include <limits>
#include <algorithm>
#include "vector.hpp"

namespace lvlset {

    /// Class defining a box used in ray tracing optimisation
    template<class T,int D> class box {
        vec<T,D> xMin,xMax;
    public:

        explicit box(const vec<T,D>&,const vec<T,D>&);

        box(const box<T,D>&);

        explicit box(const vec<T,D>&);

        box() {};

        box<T,D>& operator=(const box<T,D>&);

        ///Checks whether there are points in the box and returns result.
        bool is_empty() const {
            for (int i=0;i<D;i++) if (xMax[i]<xMin[i]) return true;
            return false;
        }

        ///Iterator over all grid points, contained by a box.
        class iterator {
            vec<T,D> pos;
            const box<T,D>& b;

        public:
            iterator(const box<T,D>& bx) : pos(bx.min()), b(bx) {}

            iterator& operator++() {
                int i;
                for (i=0;i<D-1;i++) {
                    pos[i]++;
                    if (pos[i]<=b.xMax[i]) {
                        break;
                    } else {
                        pos[i]=b.xMin[i];
                    }
                }
                if (i==D-1) pos[i]++;
                return *this;
            }

            iterator operator++(int) {
                iterator tmp=*this;
                ++(*this);
                return tmp;
            }

            bool is_finished() const {
                return (pos[D-1]>b.xMax[D-1]);
            }

            const vec<T,D>& operator*() const {
                return pos;
            }


        };

        ///Returns vector of lowest point of box in each dimension.
        const vec<T,D>& min() const {
            return xMin;
        }
        ///Returns vector of highest point of box in each dimension.
        const vec<T,D>& max() const {
            return xMax;
        }
    };

    ///Sets xMin, xMax to the lowest and highest value of the two vectors for each dimension respectively.
    template<class T,int D> inline box<T,D>::box(const vec<T,D>& idx0,const vec<T,D>& idx1) {
        xMin=Min(idx0,idx1);
        xMax=Max(idx0,idx1);
    }

    ///Creates a copy of the passed box.
    template<class T,int D> inline box<T,D>::box(const box<T,D>& box) {
        xMin=box.xMin;
        xMax=box.xMax;
    }

    ///Sets both xMin and Xmax to the passed vector.
    template<class T,int D> inline box<T,D>::box(const vec<T,D>& idx) {
        xMin=idx;
        xMax=idx;
    }

    template<class T,int D> inline box<T,D>& box<T,D>::operator=(const box<T,D>& b)  {
        xMin=b.xMin;
        xMax=b.xMax;
        return *this;
    }
}


#endif /*BOX_HPP_*/
