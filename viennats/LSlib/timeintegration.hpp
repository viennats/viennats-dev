#ifndef TIMEINTEGRATION_HPP_
#define TIMEINTEGRATION_HPP_

/* =========================================================================
   Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.

                            -----------------
                 ViennaTS - The Vienna Topography Simulator
                            -----------------

   Contact:         viennats@iue.tuwien.ac.at

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


#include <cmath>
#include <algorithm>
#include <cassert>
#include <vector>
#include "math.hpp"
#include "integration_schemes.hpp"
#include "message.h"

#include "operations.hpp"

namespace lvlset {

    using namespace math;

     namespace {
        template <class A> class PointerAdapter {
        public:
            typedef A result;
            static const A& deref(const A& x ) {return x;}
            static A& deref(A& x ) {return x;}
        };

        template <class A> class PointerAdapter<A*> {
        public:
            typedef A result;
            static const A& deref(const A* x ) {return *x;}
            static A& deref(A* x ) {return *x;}
        };
     }  //end anonymous namespace





    template <class LevelSetsType, class IntegrationSchemeType, class TimeStepRatioType, class TempRatesStopsType, class VelocityClassType, class SegmentationType>
    typename PointerAdapter<typename LevelSetsType::value_type>::result::value_type get_max_time_step(
            LevelSetsType& LevelSets,
            TempRatesStopsType& TempRatesStops,
            TimeStepRatioType TimeStepRatio,
            const VelocityClassType& Velocities,
            const IntegrationSchemeType& IntegrationScheme,
            const SegmentationType& seg

    ) {

        if (TimeStepRatio>0.4999) TimeStepRatio=0.4999;        //restriction of the CFL-Condition "TimeStepRatio" to values less than 0.5
                                                               //that means that the maximum distance which the surface moves within one time step
                                                               //is less than 0.499*GridSpacing

        typedef PointerAdapter<typename LevelSetsType::value_type> ptr;
        typedef typename ptr::result LevelSetType;
        typedef typename LevelSetType::value_type value_type;
//        typedef typename LevelSetType::size_type size_type;

        assert(std::numeric_limits<value_type>::is_iec559);
        assert(std::numeric_limits<value_type>::has_infinity);


        LevelSetType & LevelSet=ptr::deref(LevelSets.back());         //top level set



        for (typename LevelSetsType::const_iterator it=LevelSets.begin();&(*it)!=&(LevelSets.back());++it) assert(ptr::deref(*it).number_of_layers()>=2);

        //return value maximum possible time increment initialized with infinite value
        value_type MaxTimeStep2=std::numeric_limits<value_type>::max();

        //typename LevelSetType::points_type
        //seg=LevelSet.get_new_segmentation();      //TODO

        TempRatesStops.resize(seg.size()+1, 0);

        //#pragma omp parallel for schedule(static,1)
        //for (int p=0;p<=static_cast<int>(seg.size());++p) {

        #pragma omp parallel num_threads(seg.size()+1)  //use num_threads(seg.size()+1) threads
        {
            int p=0;
            #ifdef _OPENMP
            p=omp_get_thread_num();
            #endif

            TempRatesStops[p]=new std::vector<std::pair<value_type, value_type> >;
            TempRatesStops[p]->reserve(LevelSet.num_active_pts()*allocation_factor/(seg.size()+1));   //TODO

            value_type MaxTimeStep=std::numeric_limits<value_type>::max();

            typename lvlset::IntegrationScheme<LevelSetType, VelocityClassType, IntegrationSchemeType> scheme(LevelSet, Velocities, IntegrationScheme);


            typename LevelSetType::point_type start_v=(p==0)?LevelSet.grid().min_point_index():seg[p-1];
            typename LevelSetType::point_type end_v=(p!=static_cast<int>(seg.size()))?seg[p]:LevelSet.grid().increment_indices(LevelSet.grid().max_point_index());

            //iterators which iterate simultaneously over level sets
      std::vector< typename LevelSetType::const_iterator_runs> ITs;
      for (typename LevelSetsType::const_iterator it=LevelSets.begin();&(*it)!=&(LevelSets.back());++it)  ITs.push_back(typename LevelSetType::const_iterator_runs(ptr::deref(*it),start_v));


            //iteration
            for (typename LevelSetType::const_iterator_runs srfIT(LevelSet, start_v);srfIT.start_indices()<end_v;srfIT.next()) {
                if (!srfIT.is_active()) continue;

                value_type Phi_im1=srfIT.value();

                //remaining distance
                value_type q=TimeStepRatio;

                value_type tmp_tmax=0;

                typename LevelSetsType::size_type i=LevelSets.size()-1;
                while(i!=0) {

                    //calculate increment rate for actual Material i
                    value_type v=scheme(srfIT, LevelSets.size()-1-i);

                    //if rate is positive
                    if (v>0.) {
                        tmp_tmax+=q/v;
                        TempRatesStops[p]->push_back(std::make_pair(v,-std::numeric_limits<value_type>::max()));
                        break;
                    }

                    ITs[i-1].go_to_indices_sequential(srfIT.start_indices());

                    value_type Phi_i=Phi_im1;
                    Phi_im1=ITs[i-1].value();

                    if (v==0.) {
                        if (Phi_i<Phi_im1) {
                            tmp_tmax=std::numeric_limits<value_type>::max();
                            TempRatesStops[p]->push_back(std::make_pair(v,std::numeric_limits<value_type>::max()));
                            break;
                        }
                    } else {

                        if (Phi_i<Phi_im1) {
                            value_type tmp=Phi_im1-Phi_i;
                            if (tmp>=q) {
                                tmp_tmax-=q/v;
                                TempRatesStops[p]->push_back(std::make_pair(v,std::numeric_limits<value_type>::max()));
                                break;
                            } else {
                                tmp_tmax-=tmp/v;
                                TempRatesStops[p]->push_back(std::make_pair(v,Phi_im1));
                                q-=tmp;
                            }
                        }
                    }
                    i--;
                }

                if (i==0) {
                    value_type v=scheme(srfIT,LevelSets.size()-1);
                    if (v>0.) {
                        tmp_tmax+=q/v;
                    } else if (v==0.) {
                        tmp_tmax=std::numeric_limits<value_type>::max();
                    } else {
                        tmp_tmax-=q/v;
                    }
                    TempRatesStops[p]->push_back(std::make_pair(v,std::numeric_limits<value_type>::max()));
                }

                if (tmp_tmax<MaxTimeStep) MaxTimeStep=tmp_tmax;

            }

            #pragma omp critical  //execute as single thread
            {
                if (MaxTimeStep<MaxTimeStep2) MaxTimeStep2=MaxTimeStep;
            }

        }
        assert(MaxTimeStep2>0.);

        return MaxTimeStep2;

    }



    template <class GridTraitsType, class LevelSetTraitsType, class TempRatesStopsType, class TimeStepType, class SegmentationType>
    void time_integrate_active_grid_points(
                    levelset<GridTraitsType, LevelSetTraitsType>& LS,       //single level set function
                    TimeStepType  TimeStep,
                    const TempRatesStopsType& TempRatesStops,
                    const SegmentationType& seg
                    ) {


        //this function integrates all active grid points of the top most level set function over the time given by "TimeStep"

        typedef levelset<GridTraitsType, LevelSetTraitsType> LevelSetType;
        typedef typename LevelSetType::value_type value_type;
        typedef typename LevelSetType::size_type size_type;

        assert(std::numeric_limits<value_type>::has_infinity);
        assert(std::numeric_limits<value_type>::is_iec559);

        //typename LevelSetType::points_type seg=LS.get_new_segmentation();      //TODO

        //#pragma omp parallel for schedule(static,1)
        //for (int p=0;p<=static_cast<int>(seg.size());++p) {

        #pragma omp parallel num_threads(seg.size()+1) //use num_threads(seg.size()+1) threads
        {

            int p=0;
            #ifdef _OPENMP
            p=omp_get_thread_num();
            #endif

            //counter to access right entries of gradient arrays
            size_type counter=0;
             //TempRatesStopsType::value_type::const_iterator
            typename std::vector<std::pair<value_type, value_type> >::const_iterator itRS=TempRatesStops[p]->begin();

            typename LevelSetType::point_type start_v=(p==0)?LS.grid().min_point_index():seg[p-1];
            typename LevelSetType::point_type end_v=(p!=static_cast<int>(seg.size()))?seg[p]:LS.grid().increment_indices(LS.grid().max_point_index());

            //iteration
            for (typename LevelSetType::const_iterator_runs srfIT(LS, start_v);srfIT.start_indices()<end_v;srfIT.next()) {
                if (!srfIT.is_active()) continue;

                value_type Phi=srfIT.value();

                TimeStepType t=TimeStep;

                while (math::abs(itRS->second-Phi)<math::abs(t*itRS->first)) {
                    t-=math::abs((itRS->second-Phi)/itRS->first);
                    Phi=itRS->second;
                    ++itRS;
                }

                Phi-=t*itRS->first;

               LS.set_value(srfIT.pt_id(),Phi);

                while(math::abs(itRS->second)!=std::numeric_limits<value_type>::max()) ++itRS;

                ++itRS;

                ++counter;
            }

            assert(itRS==TempRatesStops[p]->end());
            delete TempRatesStops[p];

        }


    }

/*
    template <class LevelSetsType, class IntegrationSchemeType, class TimeStepRatioType, class VelocityClassType,class TimeStepType>
    typename PointerAdapter<typename LevelSetsType::value_type>::result::value_type time_integrate_active_grid_points2(
            LevelSetsType& LevelSets,
            TimeStepRatioType TimeStepRatio,
            const VelocityClassType& Velocities,
            const IntegrationSchemeType& IntegrationScheme,
            const bool separate_materials,
            //const SegmentationType& seg,
            TimeStepType MaxTimeStep2=std::numeric_limits<TimeStepType>::max()         //TODO change to value_type

    ) {
        if (TimeStepRatio>0.4999) TimeStepRatio=0.4999;        //restriction of the CFL-Condition "TimeStepRatio" to values less than 0.5
                                                               //that means that the maximum distance which the surface moves within one time step
                                                               //is less than 0.499*GridSpacing

        typedef PointerAdapter<typename LevelSetsType::value_type> ptr;
        typedef typename ptr::result LevelSetType;
        typedef typename LevelSetType::value_type value_type;
        typedef typename LevelSetType::size_type size_type;

        assert(std::numeric_limits<value_type>::is_iec559);
        assert(std::numeric_limits<value_type>::has_infinity);


        LevelSetType & LevelSet=ptr::deref(LevelSets.back());         //top level set

        typename LevelSetType::points_type seg=LevelSet.get_new_segmentation();

        for (typename LevelSetsType::const_iterator it=LevelSets.begin();&(*it)!=&(LevelSets.back());++it) assert(ptr::deref(*it).number_of_layers()>=2);

        LevelSetType new_lvlset(LevelSet.grid());

        new_lvlset.initialize(seg, LevelSet.get_allocation()*(double(1)/LevelSet.number_of_layers()));

        std::vector<size_type> active_pt_count(seg.size()+1, LevelSetType::INACTIVE);

        #pragma omp parallel num_threads(seg.size()+1)
        {
            int p=0;
            #ifdef _OPENMP
            p=omp_get_thread_num();
            #endif

            std::vector<std::pair<value_type, value_type> > TempRatesStops;

            TempRatesStops.reserve(LevelSet.num_active_pts()*allocation_factor/(seg.size()+1));   //TODO

            value_type MaxTimeStep=std::numeric_limits<value_type>::max();
            //std::cout << "1.\n";
            typename lvlset::IntegrationScheme<LevelSetType, VelocityClassType, IntegrationSchemeType> scheme(LevelSet, Velocities, IntegrationScheme);

            //iterators which iterate simultaneously over level sets
            std::vector<typename LevelSetType::const_iterator_runs> ITs;
            for (typename LevelSetsType::const_iterator it=LevelSets.begin();&(*it)!=&(LevelSets.back());++it)  ITs.push_back(typename LevelSetType::const_iterator_runs(ptr::deref(*it)));

            typename LevelSetType::point_type start_v=(p==0)?LevelSet.grid().min_point_index():seg[p-1];
            typename LevelSetType::point_type end_v=(p!=static_cast<int>(seg.size()))?seg[p]:LevelSet.grid().increment_indices(LevelSet.grid().max_point_index());

            //iteration
            int iterator_int=0.;
            for (typename LevelSetType::const_iterator_runs srfIT(LevelSet, start_v);srfIT.start_indices()<end_v;srfIT.next()) {
                if (!srfIT.is_active()) {
                    assert(math::abs(srfIT.value())>0.5);
                    new_lvlset.push_back_undefined(p,srfIT.start_indices(),(srfIT.sign()==POS_SIGN)?LevelSetType::POS_PT:LevelSetType::NEG_PT);
                    continue;
                }
                if (active_pt_count[p]==LevelSetType::INACTIVE) active_pt_count[p]=srfIT.active_pt_id();


                value_type Phi_im1=srfIT.value();

                assert(math::abs(srfIT.value())<=0.5);

                new_lvlset.push_back(p,srfIT.start_indices(),Phi_im1);

                //remaining distance
                value_type q=TimeStepRatio;

                value_type tmp_tmax=0;

                typename LevelSetsType::size_type i=LevelSets.size()-1;
                while(i!=0) {

                  //need to deal with two types of velocities. One modifies the surface
                  //and one moves the entire LS surface in a given direction

                    //calculate increment rate for actual Material i
                    value_type v=scheme(srfIT, LevelSets.size()-1-i, 0.);

                    ITs[i-1].go_to_indices_sequential(srfIT.start_indices());

                    value_type Phi_i=Phi_im1;
                    Phi_im1=ITs[i-1].value();

                    //if rate is positive
                    if (v>0.) {
                      if (separate_materials) {
                        if (Phi_i>-Phi_im1) {
                          value_type tmp=Phi_im1+Phi_i;
                          if (tmp>=q) {
                            tmp_tmax+=q/v;
                            TempRatesStops.push_back(std::make_pair(v,-std::numeric_limits<value_type>::max()));
                            break;
                          } else {
                            tmp_tmax+=tmp/v;
                            TempRatesStops.push_back(std::make_pair(v,-Phi_im1));
                            q-=tmp;
                          }
                        }
                      } else {
                            if (Phi_i<Phi_im1) {
                                tmp_tmax=std::numeric_limits<value_type>::max();
                                TempRatesStops.push_back(std::make_pair(v,std::numeric_limits<value_type>::max()));
                                break;
                            }
                      }
                    } else if (v==0.) {
                        if (Phi_i<Phi_im1) {
                            tmp_tmax=std::numeric_limits<value_type>::max();
                            TempRatesStops.push_back(std::make_pair(v,std::numeric_limits<value_type>::max()));
                            break;
                        }
                    } else {

                        if (Phi_i<Phi_im1) {
                            value_type tmp=Phi_im1-Phi_i;

                            if (tmp>=q) {
                                tmp_tmax-=q/v;
                                TempRatesStops.push_back(std::make_pair(v,std::numeric_limits<value_type>::max()));
                                break;
                            } else {
                                tmp_tmax-=tmp/v;
                                TempRatesStops.push_back(std::make_pair(v,Phi_im1));
                                q-=tmp;
                            }
                        }
                    }
                    i--;
                }

                if (i==0) {
                  //std::cout << "iterator_int: " << iterator_int << std::endl;
                  iterator_int++;
                    value_type v=scheme(srfIT,LevelSets.size()-1,iterator_int);
                    //if (v>-0.9){
                    //  std::cout << "v(" << i << "): " << v << std::endl;
                    //  std::cout << "q: " << q << std::endl;
                    //}
                    //std::cout << "v(-)ti: " << v << std::endl;

                    if (v>0.) {
                        tmp_tmax+=q/v;
                    } else if (v==0.) {
                        tmp_tmax=std::numeric_limits<value_type>::max();
                    } else {
                        tmp_tmax-=q/v;
                    }
                    TempRatesStops.push_back(std::make_pair(v,std::numeric_limits<value_type>::max()));
                }

                if (tmp_tmax<MaxTimeStep) MaxTimeStep=tmp_tmax;

            }

            //determine max time step
            #pragma omp critical
            {
                if (MaxTimeStep<MaxTimeStep2) MaxTimeStep2=MaxTimeStep;
            }
            #pragma omp barrier
            #pragma omp single
            {
                new_lvlset.finalize(1);
                assert(new_lvlset.num_active_pts()==LevelSet.num_active_pts());
                assert(new_lvlset.num_pts()==new_lvlset.num_active_pts());
                LevelSet.swap(new_lvlset);
            }

            //update level set values
            typename std::vector<std::pair<value_type, value_type> >::const_iterator itRS=TempRatesStops.begin();
            for (size_type local_pt_id=0;local_pt_id<LevelSet.num_pts(p);++local_pt_id) {

                value_type Phi=LevelSet.value(p, local_pt_id);

                TimeStepType t=MaxTimeStep2;

                while (math::abs(itRS->second-Phi)<math::abs(t*itRS->first)) {
                    t-=math::abs((itRS->second-Phi)/itRS->first);
                    Phi=itRS->second;
                    ++itRS;
                }

                Phi-=t*itRS->first;

                LevelSet.set_value(p, local_pt_id, Phi);

                while(math::abs(itRS->second)!=std::numeric_limits<value_type>::max()) ++itRS;
                ++itRS;
            }


            assert(itRS==TempRatesStops.end());

        }

        return MaxTimeStep2;

    }
*/

// SFINAE (Substitution Failure Is Not An Error): If IntegrationScheme is STENCIL_LOCAL_LAX_FRIEDRICHS_SCALAR_TYPE the time step is reduced depending on dissipation coefficients.
template<class LevelSetType,class VelocityClassType,class IntegrationSchemeType, class TimeStepType,
         typename std::enable_if< std::is_same< lvlset::STENCIL_LOCAL_LAX_FRIEDRICHS_SCALAR_TYPE, IntegrationSchemeType>::value>::type* = nullptr>
         void reduce_timestep_hamilton_jacobi( lvlset::IntegrationScheme<LevelSetType, VelocityClassType, IntegrationSchemeType>& scheme, TimeStepType& MaxTimeStep) {

  typedef typename LevelSetType::value_type value_type;

  const double alpha_maxCFL = 1.0; //TODO Can be potentially smaller than 1 (user input???)
  //second time step test, based on alphas
  vec<value_type,3> alphas = scheme.getFinalAlphas();
  vec<value_type,3> dxs = scheme.getDx();

  MaxTimeStep=0;
  for(int i = 0; i < 3; ++i){
    if(math::abs(dxs[i]) > 1e-6 ){
      MaxTimeStep += alphas[i] / dxs[i];
    }
  }

  MaxTimeStep = alpha_maxCFL / MaxTimeStep;
}

// SFINAE (Substitution Failure Is Not An Error): IntegrationScheme != STENCIL_LOCAL_LAX_FRIEDRICHS_SCALAR_TYPE
template<class LevelSetType,class VelocityClassType,class IntegrationSchemeType, class TimeStepType,
         typename std::enable_if< !std::is_same< lvlset::STENCIL_LOCAL_LAX_FRIEDRICHS_SCALAR_TYPE, IntegrationSchemeType>::value>::type* = nullptr>
         void reduce_timestep_hamilton_jacobi( lvlset::IntegrationScheme<LevelSetType, VelocityClassType, IntegrationSchemeType>& scheme, TimeStepType& MaxTimeStep2) {}


    template <class LevelSetsType, class IntegrationSchemeType, class TimeStepRatioType, class VelocityClassType,class TimeStepType>
    typename PointerAdapter<typename LevelSetsType::value_type>::result::value_type time_integrate_active_grid_points2(
            LevelSetsType& LevelSets,
            TimeStepRatioType TimeStepRatio,
            const VelocityClassType& Velocities,
            const IntegrationSchemeType& IntegrationScheme,
//            const bool separate_materials,
            //const SegmentationType& seg,
            TimeStepType MaxTimeStep2=std::numeric_limits<TimeStepType>::max()        //TODO change to value_type
    ) {
        if (TimeStepRatio>0.4999) TimeStepRatio=0.4999;        //restriction of the CFL-Condition "TimeStepRatio" to values less than 0.5
                                                               //that means that the maximum distance which the surface moves within one time step
                                                               //is less than 0.499*GridSpacing

        typedef PointerAdapter<typename LevelSetsType::value_type> ptr;
        typedef typename ptr::result LevelSetType;
        typedef typename LevelSetType::value_type value_type;
        typedef typename LevelSetType::size_type size_type;

        assert(std::numeric_limits<value_type>::is_iec559);
        assert(std::numeric_limits<value_type>::has_infinity);


        LevelSetType & LevelSet=ptr::deref(LevelSets.back());         //top level set

        typename LevelSetType::points_type seg=LevelSet.get_new_segmentation();

        for (typename LevelSetsType::const_iterator it=LevelSets.begin();&(*it)!=&(LevelSets.back());++it) assert(ptr::deref(*it).number_of_layers()>=2);

        // the levelset filled with the new values
        LevelSetType new_lvlset(LevelSet.grid());

        new_lvlset.initialize(seg, LevelSet.get_allocation()*(double(1)/LevelSet.number_of_layers()));

        std::vector<size_type> active_pt_count(seg.size()+1, LevelSetType::INACTIVE);

        // split over threads via segmentations
        #pragma omp parallel num_threads(seg.size()+1) //use num_threads(seg.size()+1) threads
        {
            int p=0;
            #ifdef _OPENMP
            p=omp_get_thread_num();
            #endif

            // holds the advection velocities which will be applied later
            std::vector<std::pair<value_type, value_type> > TempRatesStops;
            // allocation factor is ~1.2 which is the maximum amount of new points one would expect
            TempRatesStops.reserve(LevelSet.num_active_pts()*allocation_factor/(seg.size()+1));   //TODO

            value_type MaxTimeStep=std::numeric_limits<value_type>::max();

            // initialise integration by passing all found velocities to the scheme class
            typename lvlset::IntegrationScheme<LevelSetType, VelocityClassType, IntegrationSchemeType> scheme(LevelSet, Velocities, IntegrationScheme);

            //iterators which iterate simultaneously over level sets
            std::vector<typename LevelSetType::const_iterator_runs> ITs;
            for (typename LevelSetsType::const_iterator it=LevelSets.begin();&(*it)!=&(LevelSets.back());++it){
              ITs.push_back(typename LevelSetType::const_iterator_runs(ptr::deref(*it)));
            }

            // set start and end indices in the grid for each thread
            typename LevelSetType::point_type start_v=(p==0)?LevelSet.grid().min_point_index():seg[p-1];
            typename LevelSetType::point_type end_v=(p!=static_cast<int>(seg.size()))?seg[p]:LevelSet.grid().increment_indices(LevelSet.grid().max_point_index());

            // initialise top levelset iterator
            ITs.push_back(typename LevelSetType::const_iterator_runs(LevelSet, start_v));

            // iterate over the LS points in the segmentation of each thread
            for (typename LevelSetType::const_iterator_runs& srfIT=ITs.back(); srfIT.start_indices()<end_v; srfIT.next()) {

                // if LS point is not an active point, push an undefined point to new levelset
                if (!srfIT.is_active()) {
                    assert(math::abs(srfIT.value())>0.5);
                    new_lvlset.push_back_undefined(p,srfIT.start_indices(),(srfIT.sign()==POS_SIGN)?LevelSetType::POS_PT:LevelSetType::NEG_PT);
                    continue;
                }


                if (active_pt_count[p]==LevelSetType::INACTIVE) active_pt_count[p]=srfIT.active_pt_id();

                // Phi_im1 holds the LS value of the current grid point
                value_type Phi_im1=srfIT.value();

                assert(math::abs(srfIT.value())<=0.5);

                // push current active point to new_lvlset
                new_lvlset.push_back(p,srfIT.start_indices(),Phi_im1);

                // q = cfl condition
                value_type q=TimeStepRatio;
                value_type tmp_tmax=0;

                for(typename LevelSetsType::size_type i=LevelSets.size()-1; i>=0; --i){

                  //value_type v=scheme(srfIT, 0);  // rate of topmost levelset
                  value_type v = 0;

                  // check if there is any other levelset at the same point:
                  // if yes, take the velocity of the lowest levelset
                  for(unsigned level_num=0; level_num<LevelSets.size(); ++level_num){
                    // put iterator to same position as the top levelset
                    ITs[level_num].go_to_indices_sequential(srfIT.start_indices());

                    // if the lower surface is actually outside, i.e. its LS value is lower or equal
                    if(ITs[level_num].value() <= Phi_im1){
                      v = scheme(srfIT, LevelSets.size()-1-level_num); //NOTE AT numerical hamiltonian is calculated
                      break;
                    }
                  }

                  // Phi_i LS value of material i
                  value_type Phi_i=Phi_im1;

                  // find LS point of material i-1, which is below i
                  if(i>0){
                    ITs[i-1].go_to_indices_sequential(srfIT.start_indices());
                    Phi_im1=ITs[i-1].value();
                  }else Phi_im1 = std::numeric_limits<value_type>::max();


                  // if velocity is positive, set maximum time step possible without violating the cfl condition
                  if (v>0.) {
                    tmp_tmax+=q/v;
                    TempRatesStops.push_back(std::make_pair(v,-std::numeric_limits<value_type>::max()));
                    break;
                    // if velocity is 0, maximum time step is infinite
                  } else if (v==0.) {
                    tmp_tmax=std::numeric_limits<value_type>::max();
                    TempRatesStops.push_back(std::make_pair(v,std::numeric_limits<value_type>::max()));
                    break;
                    // if the velocity is negative apply the velocity for as long as possible without infringing on material below
                  } else {
                    value_type tmp= math::abs(Phi_im1-Phi_i);

                    if (tmp>=q) {
                      tmp_tmax-=q/v;
                      TempRatesStops.push_back(std::make_pair(v,std::numeric_limits<value_type>::max()));
                      break;
                    } else {
                      tmp_tmax-=tmp/v;
                      // the second part of the pair indicates how far we can move in this time step until the end of the material is reached
                      TempRatesStops.push_back(std::make_pair(v,Phi_im1));
                      q-=tmp;
                    }
                  }
                }

                if (tmp_tmax<MaxTimeStep) MaxTimeStep=tmp_tmax;

            }

            //determine max time step
            #pragma omp critical  //execute as single thread
            {
                if (MaxTimeStep<MaxTimeStep2) MaxTimeStep2=MaxTimeStep;

                //If scheme is STENCIL_LOCAL_LAX_FRIEDRICHS the time step is reduced depending on the dissipation coefficients
                //For all remaining schemes this function is empty.
                reduce_timestep_hamilton_jacobi(scheme, MaxTimeStep);

                if (MaxTimeStep<MaxTimeStep2){
                   #ifdef VERBOSE
                    std::cout << "HJ: Reduced time step from " << MaxTimeStep2;
                    std::cout << " to " << MaxTimeStep << std::endl;
                   #endif
                   MaxTimeStep2=MaxTimeStep;
                }
            }
            #pragma omp barrier //wait until all other threads in section reach the same point.
            #pragma omp single //section of code that must be run by a single available thread.
            {



                new_lvlset.finalize(1);
                assert(new_lvlset.num_active_pts()==LevelSet.num_active_pts());
                assert(new_lvlset.num_pts()==new_lvlset.num_active_pts());
                LevelSet.swap(new_lvlset);
            }

            /*
                This is where the magic happens
                Here the velocites are applied to the LS values
            */

            // iterator over all velocities in TempStopRates to apply them
            typename std::vector<std::pair<value_type, value_type> >::const_iterator itRS=TempRatesStops.begin();

            // iterate over all points in the segmentation of new topmost levelset
            for (size_type local_pt_id=0;local_pt_id<LevelSet.num_pts(p);++local_pt_id) {
                // phi is the LS value at the current grid point
                value_type Phi=LevelSet.value(p, local_pt_id);

                TimeStepType t=MaxTimeStep2;  // maximum time step we can take

                // if there is a change in materials during one time step, deduct the time taken to advect up to the end of the top material and set the LS value to the one below
                while (math::abs(itRS->second-Phi)<math::abs(t*itRS->first)) {
                    t-=math::abs((itRS->second-Phi)/itRS->first);
                    Phi=itRS->second;
                    ++itRS;
                }


                //NOTE AT now timestep is define, can be compared with SLF dt


                // now deduct the velocity times the time step we take
                Phi-=t*itRS->first;

                // set the new LS value for the new layer
                LevelSet.set_value(p, local_pt_id, Phi);

                // this is run when two materials are close but the velocity is too slow to actually reach the second material, to get read of the extra entry in the TempRatesStop
                while(math::abs(itRS->second)!=std::numeric_limits<value_type>::max()) ++itRS;

                // advance the TempStopRates iterator by one
                ++itRS;
            }

            assert(itRS==TempRatesStops.end());

        }

        return MaxTimeStep2;

    }

    namespace {

        template <class VelocityType, class LevelSetType> class VelocityAdapter {       //this velocity adapter is used
            const VelocityType& vc;                                                     //to make a binary velocity-functor from an
        public:                                                                         //unary one, in case of time integration of
                                                                                        //a single level set function
            VelocityAdapter(const VelocityType& v) : vc(v) {}

            template <class A, class B>
            typename LevelSetType::value_type operator()(A a, B b) const {
                return vc(a);
            }
        };


        template <class PointDataType, class PointDataSizeType>
        class DataCopyType {
            PointDataType & data;
            PointDataType old_data;
            const PointDataSizeType data_size;
        public:

            DataCopyType(PointDataType & d, PointDataSizeType ds):data(d), data_size(ds) {
                std::swap(old_data, data);
            }

            template <class I>
            void set_size(I i) const {
                data.resize(data_size*i);
            }

            template <class I>
            void copy_data(I i, I j) const {
                //assert((j+1)*data_size<=old_data.size());
                //assert((i+1)*data_size<=data.size());
                std::copy(old_data.begin()+j*data_size, old_data.begin()+(j+1)*data_size, data.begin()+i*data_size);
            }
        };
    }

    template <class GridTraitsType, class LevelSetTraitsType, class TempRatesStopsType, class IntegrationSchemeType, class TimeStepRatioType, class VelocityClassType, class SegmentationType>
    typename levelset<GridTraitsType, LevelSetTraitsType>::value_type  get_max_time_step(
                    levelset<GridTraitsType, LevelSetTraitsType>& LS,           //single level set function
                    TempRatesStopsType& TempRatesStop,
                    TimeStepRatioType TimeStepRatio,
                    const VelocityClassType& v,
                    const IntegrationSchemeType& IntegrationScheme,
                    const SegmentationType& seg) {

         //specialization of get_max_time_step for a single level set function

        typedef levelset<GridTraitsType, LevelSetTraitsType> LevelSetType;

        std::vector<LevelSetType*> tmp(1,&LS);

        return get_max_time_step(tmp, TempRatesStop, TimeStepRatio, v, IntegrationScheme, seg); //VelocityAdapter<VelocityType, LevelSetType>(Velocities));

    }


    template <class LevelSetsType, class TempRatesStopsType, class TimeStepType, class SegmentationType>
    void time_integrate_active_grid_points(
            LevelSetsType& LS,
            TimeStepType TimeStep,
            const TempRatesStopsType& TempRatesStops,
            const SegmentationType& seg
    ) {

        //specialization of time_integrate_active_grid_points for multiple level set functions

        typedef PointerAdapter<typename LevelSetsType::value_type> ptr;

        time_integrate_active_grid_points(ptr::deref(LS.back()), TimeStep, TempRatesStops, seg);

    }



    template <  class LevelSetsType,
                class VelocitiesClass,
                class IntegrationSchemeType,
                class CFLType,
                class TimeType,
                class PointDataType,
                class PointDataSizeType>
    TimeType time_integrate(    LevelSetsType& LevelSets,
                                const VelocitiesClass& Velocities,
                                const IntegrationSchemeType& integration_scheme,
                                CFLType CFL,
                                TimeType MaxTimeStep,
                                PointDataType& PointData,
                                PointDataSizeType PointDataSize=1,
                                const bool is_selective_depo=false
                            ) {

//    msg::print_message("Velocity: ", Velocities.getVelocity());
        //this function integrates a the level set functions "LevelSets" over time
        //these level set functions define different regions of materials
        //only the top most level set function is time integarated
        //however, the different surface velocities in different material regions
        //are consistently taken into account.
        //see also the definition of "time_integrate_active_grid_points"-function for more details

        typedef PointerAdapter<typename LevelSetsType::value_type> ptr;
        typedef typename ptr::result LevelSetType;
//        typedef typename LevelSetType::value_type value_type;
//        typedef typename LevelSetType::size_type size_type;

        //std::vector<std::vector<std::pair<value_type, value_type> > > TempRS;        //vector to store rates and stops for multi level set
        //TempRS.reserve(ptr::deref(LevelSets.back()).num_active_pts());
        //std::vector<std::vector<std::pair<value_type, value_type> >* > TempRS;

        //double tmp;

        //tmp=-my::time::GetTime();
        IntegrationScheme<LevelSetType, VelocitiesClass, IntegrationSchemeType>::prepare_surface_levelset(ptr::deref(LevelSets.back()));
        //tmp+=my::time::GetTime();
        //std::cout << tmp;




        /*typename LevelSetType::points_type seg=ptr::deref(LevelSets.back()).get_new_segmentation();

        tmp=-my::time::GetTime();
        time_step=get_max_time_step(           //determine maximum possible time step
                LevelSets,
                TempRS,
                CFL,
                Velocities,
                integration_scheme,
                seg);
        tmp+=my::time::GetTime();
        std::cout << " " << tmp;

        tmp=-my::time::GetTime();
        ptr::deref(LevelSets.back()).reduce(1);         //reduce the top most level set function to one layer
        tmp+=my::time::GetTime();
        std::cout << " " << tmp;


        if (time_step>MaxTimeStep) time_step=MaxTimeStep;

        tmp=-my::time::GetTime();
        time_integrate_active_grid_points(  LevelSets,
                                            time_step,
                                            TempRS,
                                            seg);

        tmp+=my::time::GetTime();
        std::cout << " " << tmp;*/

        //tmp=-my::time::GetTime();
        TimeType time_step=time_integrate_active_grid_points2(           //determine maximum possible time step
                                                LevelSets,
                                                CFL,
                                                Velocities,
                                                integration_scheme,
                                                MaxTimeStep);


            //tmp+=my::time::GetTime();
        //std::cout << " " << tmp;

        /*tmp=-my::time::GetTime();
        ptr::deref(LevelSets.back()).reduce(1);         //reduce the top most level set function to one layer
        tmp+=my::time::GetTime();
        std::cout << " " << tmp;*/


        //tmp=-my::time::GetTime();
        ptr::deref(LevelSets.back()).rebuild(DataCopyType<PointDataType, PointDataSizeType>(PointData, PointDataSize));     //rebuild the top most level set function,
                                                                            //and also map data assigned to each active grid point
                                                                            //to the new level set function
        //tmp+=my::time::GetTime();
        //std::cout << " " << tmp;

        //tmp=-my::time::GetTime();

        for (typename LevelSetsType::iterator it=LevelSets.begin();&(*it)!=&(LevelSets.back());++it) {
            //adjust all level set functions below the top most level set function
            //For selective deposition this is not necessary.
            if(is_selective_depo == false){
                ptr::deref(*it).max(ptr::deref(LevelSets.back()));
            }

            ptr::deref(*it).prune();            //remove grid points which do not have at least one opposite signed neighbor
            ptr::deref(*it).segment();
        }

        //tmp+=my::time::GetTime();
        //std::cout << " " << tmp << std::endl;

        return time_step;       //return the time_step which was used for time integration
    }

    template <  class GridTraitsType,
                class LevelSetTraitsType,
                class VelocitiesClass,
                class IntegrationSchemeType,
                class CFLType,
                class TimeType,
                class PointDataType,
                class PointDataSizeType>
    TimeType time_integrate(    levelset<GridTraitsType, LevelSetTraitsType>& LevelSet,
                                const VelocitiesClass& Velocities,
                                const IntegrationSchemeType& integration_scheme,
                                CFLType CFL,
                                TimeType MaxTimeStep,
                                PointDataType& PointData,
                                PointDataSizeType PointDataSize=1,
                                const bool is_selective_depo=false
                            ) {
      //std::cout << "time_integrate!\n";
        //this specialization of "time_integrate" is for the time integration of just one level set function
        //it initilaizes a new container and adds a reference to the given level set functions
        //then the general "time-integrate"-function is called

        typedef levelset<GridTraitsType, LevelSetTraitsType> LevelSetType;
        std::vector<LevelSetType*> tmp(1,&LevelSet);    //setup a container with a reference to the level set function

        return time_integrate( tmp,
                                //VelocityAdapter<VelocitiesClass, LevelSetType>(Velocities)
                                Velocities,
                                integration_scheme,
                                CFL,
                                MaxTimeStep,
                                PointData,
                                PointDataSize
                        );
    }

    template <  class LevelSetsType,
                class VelocitiesClass,
                class IntegrationSchemeType,
                class CFLType,
                class TimeType>
    TimeType time_integrate(    LevelSetsType& LevelSet,
                                const VelocitiesClass& Velocities,
                                const IntegrationSchemeType& integration_scheme,
                                CFLType CFL,
                                TimeType MaxTimeStep
                                ) {
        //this version of the "time_integrate" function is used if
        //no data assigned to the active grid points should be mapped to the new level set function
        //after time integration

        std::vector<int> dummy;             //dummy data-vector
        return time_integrate(  LevelSet,
                                 Velocities,
                                 integration_scheme,
                                 CFL,
                                 MaxTimeStep,
                                 dummy,
                                 0);
    }
}



#endif /*TIMEINTEGRATION_HPP_*/
