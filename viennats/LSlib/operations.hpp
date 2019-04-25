#ifndef OPERATIONS_HPP_
#define OPERATIONS_HPP_

/* =========================================================================
   Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.

                            -----------------
                 ViennaTS - The Vienna Topography Simulator
                            -----------------

   Contact:         viennats@iue.tuwien.ac.at

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


#include "kernel.hpp"
#include "math.hpp"
#include <algorithm>
#include <cassert>

namespace lvlset {

    namespace {
        template <class GridTraitsType, class LevelSetTraitsType, class BinaryType>
        levelset<GridTraitsType, LevelSetTraitsType> MinMax(const levelset<GridTraitsType, LevelSetTraitsType>& lA, const levelset<GridTraitsType, LevelSetTraitsType>& lB, const BinaryType& min_or_max) {
            //this function is called by the min and the max function defined below
            //this function should not be called directly by the user,
            //      and is therefore within an anonymous namespace
            //depending on which binary functor is passed to the library (min or max)
            //      this function calculates the minimum or maximum of two level set functions
            //      with an optimal complexity O(lA.num_pts()+lB.num_pts())

            assert(lA.number_of_layers()>0);
            assert(lB.number_of_layers()>0);

            typedef levelset<GridTraitsType, LevelSetTraitsType> LevelSetType;

            LevelSetType tmp(lA.grid());

            tmp.initialize();

            typename LevelSetType::const_iterator_runs itA(lA);
            typename LevelSetType::const_iterator_runs itB(lB);

            while (!itA.is_finished() || !itB.is_finished()) {

                typename LevelSetType::value_type d=min_or_max(itA.value(),itB.value());
                vec<typename LevelSetType::index_type,GridTraitsType::dimensions> pos=std::max(itA.start_indices(),itB.start_indices());

                if (math::abs(d)<std::numeric_limits<typename LevelSetType::value_type>::max()) {
                    tmp.push_back(0,pos, d);      //TODO
                } else {
                    tmp.push_back_undefined(0, pos, (tmp.sign(d)==POS_SIGN)?LevelSetType::POS_PT:LevelSetType::NEG_PT);     //TODO
                }

                switch(compare(itA.end_indices(), itB.end_indices())) {
                    case -1:
                        itA.next();
                        break;
                    case 0:
                        itA.next();
                    default:
                        itB.next();
                }

            }

            tmp.finalize(std::min(lA.number_of_layers(), lB.number_of_layers()));

            return tmp;
        }

        /*template <class GridTraitsType, class LevelSetTraitsType, class BinaryType>
        void MinMax2(levelset<GridTraitsType, LevelSetTraitsType>& lA, const levelset<GridTraitsType, LevelSetTraitsType>& lB, const BinaryType& min_or_max) {
            //this function is called by the min and the max function defined below
            //this function should not be called directly by the user,
            //      and is therefore within an anonymous namespace
            //depending on which binary functor is passed to the library (min or max)
            //      this function calculates the minimum or maximum of two level set functions
            //      with an optimal complexity O(lA.num_pts()+lB.num_pts())

            assert(lA.number_of_layers()>0);
            assert(lB.number_of_layers()>0);

            typedef levelset<GridTraitsType, LevelSetTraitsType> LevelSetType;

            LevelSetType tmp(lA.grid().grid_traits());

            tmp.initialize(lA.get_new_segmentation());

            #pragma omp parallel num_threads(lA.sub_levelsets.size())
            {
                #pragma omp for schedule(static, 1)
                for (int p=0;p<static_cast<int>(lA.sub_levelsets.size());++p) {

                    typename LevelSetType::const_iterator_runs itA(lA,(p==0)?grid().min_point_index():segmentation[p-1]);
                    typename LevelSetType::const_iterator_runs itB(lB,(p==0)?grid().min_point_index():segmentation[p-1]);

                    vec<typename LevelSetType::index_type, D> end_v=(p!=static_cast<int>(lA.segmentation.size()))?segmentation[p]:lA.grid().increment_indices(grid().max_point_index());

                    while ((itA.start_indices()<end_v) || (itB.start_indices()<end_v)) {

                        typename LevelSetType::value_type d=min_or_max(itA.value(),itB.value());
                        vec<typename LevelSetType::index_type,GridTraitsType::dimensions> pos=std::max(itA.start_indices(),itB.start_indices());

                        if (math::abs(d)<std::numeric_limits<typename LevelSetType::value_type>::max()) {
                            tmp.push_back(0,pos, d);      //TODO
                        } else {
                            tmp.push_back_undefined(0, pos, (tmp.sign(d)==POS_SIGN)?LevelSetType::POS_PT:LevelSetType::NEG_PT);     //TODO
                        }

                        switch(compare(itA.end_indices(), itB.end_indices())) {
                            case -1:
                                itA.next();
                                break;
                            case 0:
                                itA.next();
                            default:
                                itB.next();
                        }
                    }
                }
            }

            tmp.finalize(std::min(lA.number_of_layers(), lB.number_of_layers()));

            tmp.swap(lA);
        }*/

        class max_wrapper {
        public:
            template <class G> G operator()(G a, G b) const {
                return std::max(a,b);
            }
        };

        class min_wrapper {
        public:
            template <class G> G operator()(G a, G b) const {
                return std::min(a,b);
            }
        };
    }

    template <class GridTraitsType, class LevelSetTraitsType>
    levelset<GridTraitsType, LevelSetTraitsType> max(const levelset<GridTraitsType, LevelSetTraitsType>& lA, const levelset<GridTraitsType, LevelSetTraitsType>& lB) {
        //this function returns the maximum of two level set functions
        return MinMax(lA, lB, max_wrapper());
    }

    template <class GridTraitsType, class LevelSetTraitsType>
    levelset<GridTraitsType, LevelSetTraitsType> min(const levelset<GridTraitsType, LevelSetTraitsType>& lA, const levelset<GridTraitsType, LevelSetTraitsType>& lB) {
        //this function returns the minimum of two level set functions
        return MinMax(lA, lB, min_wrapper());
    }




    /*template <class GridTraitsType, class LevelSetTraitsType>
    void max2(levelset<GridTraitsType, LevelSetTraitsType>& lA, const levelset<GridTraitsType, LevelSetTraitsType>& lB) {
        //this function returns the maximum of two level set functions
        MinMax2(lA, lB, max_wrapper());
    }

    template <class GridTraitsType, class LevelSetTraitsType>
    void min2(levelset<GridTraitsType, LevelSetTraitsType>& lA, const levelset<GridTraitsType, LevelSetTraitsType>& lB) {
        //this function returns the minimum of two level set functions
        MinMax2(lA, lB, min_wrapper());
    }*/


    template <class GridTraitsType, class LevelSetTraitsType>
    levelset<GridTraitsType, LevelSetTraitsType> invert(const levelset<GridTraitsType, LevelSetTraitsType>& l) {
        //this function returns an inverted level set function (the signs of all level set values are reversed)
        levelset<GridTraitsType, LevelSetTraitsType> tmp=l;
        tmp.invert();
        return tmp;
    }

    template <class GridTraitsType, class LevelSetTraitsType> bool have_equal_signs(const levelset<GridTraitsType, LevelSetTraitsType>& lA, const levelset<GridTraitsType, LevelSetTraitsType>& lB) {
        //this function returns true if the signs of all grid points are equal for the level set functions lA and lB

        typename levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_runs itA(lA);
        typename levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_runs itB(lB);

        vec<typename levelset<GridTraitsType, LevelSetTraitsType>::index_type,GridTraitsType::dimensions>  startA, startB;

        while (!itA.is_finished() || !itB.is_finished()) {

            if (itA.sign()!=itB.sign()) return false;

            if (itA.end_indices()<itB.end_indices()) {
                itA.next();
            } else if (itA.end_indices()>itB.end_indices()) {
                itB.next();
            } else {
                itA.next();
                itB.next();
            }
        }
        return true;
    }


    template <class GridTraitsType, class LevelSetTraitsType> bool operator<=(const levelset<GridTraitsType, LevelSetTraitsType>& lA, const levelset<GridTraitsType, LevelSetTraitsType>& lB) {
        //this function returns true if for all grid points the level set values of lA are smaller or equal than those of lB

        typename levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_runs itA(lA);
        typename levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_runs itB(lB);

        vec<typename levelset<GridTraitsType, LevelSetTraitsType>::index_type,GridTraitsType::dimensions> startA, startB;

        while (!itA.is_finished() || !itB.is_finished()) {

            if (itA.value()>itB.value()) {
              itA.print();              //TODO
              itB.print();
              std::cout << "diff value = " << (itA.value()-itB.value()) << std::endl;
              return false;
            }

            if (itA.end_indices()<itB.end_indices()) {
                itA.next();
            } else if (itA.end_indices()>itB.end_indices()) {
                itB.next();
            } else {
                itA.next();
                itB.next();
            }
        }
        return true;
    }

    //this function returns a*lA + b*lB
    template <class GridTraitsType, class LevelSetTraitsType> levelset<GridTraitsType, LevelSetTraitsType> numerical_linear_combination(const levelset<GridTraitsType, LevelSetTraitsType>& lA, const levelset<GridTraitsType, LevelSetTraitsType>& lB, double a, double b ) {


        typedef levelset<GridTraitsType, LevelSetTraitsType> LevelSetType;

        typename levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_runs itA(lA);
        typename levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_runs itB(lB);

        LevelSetType tmp(lA.grid());
        tmp.initialize();

        while (!itA.is_finished() || !itB.is_finished()) {

            typename LevelSetType::value_type d =  a*itA.value() + b*itB.value();
            vec<typename LevelSetType::index_type,GridTraitsType::dimensions> pos=std::max(itA.start_indices(),itB.start_indices());


            if (math::abs(d)<std::numeric_limits<typename LevelSetType::value_type>::max()) { // TODO: this is dangerous, value might overflow and give wrong results
                tmp.push_back(0,pos, d);
            } else {
                tmp.push_back_undefined(0, pos, (tmp.sign(d)==POS_SIGN)?LevelSetType::POS_PT:LevelSetType::NEG_PT);
            }

            switch(compare(itA.end_indices(), itB.end_indices())) {
                case -1:
                    itA.next();
                    break;
                case 0:
                    itA.next();
                default:
                    itB.next();
            }
        }

        tmp.finalize(std::min(lA.number_of_layers(), lB.number_of_layers()));


        return tmp;
    }

    template <class GridTraitsType, class LevelSetTraitsType> bool operator>=(const levelset<GridTraitsType, LevelSetTraitsType>& lA, const levelset<GridTraitsType, LevelSetTraitsType>& lB) {
        //this function returns true if for all grid points the level set values of lA are greater or equal than those of lB
        return (lB<=lA);
    }

}

#endif /*OPERATIONS_HPP_*/
