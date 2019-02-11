#ifndef INTEGRATION_SCHEMES_HPP_
#define INTEGRATION_SCHEMES_HPP_

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

#include "../Math.h"

namespace lvlset {




    class ENGQUIST_OSHER_SCALAR_1ST_ORDER_TYPE {};
    const ENGQUIST_OSHER_SCALAR_1ST_ORDER_TYPE ENGQUIST_OSHER_SCALAR_1ST_ORDER=ENGQUIST_OSHER_SCALAR_1ST_ORDER_TYPE();

    class ENGQUIST_OSHER_SCALAR_2ND_ORDER_TYPE {};
    const ENGQUIST_OSHER_SCALAR_2ND_ORDER_TYPE ENGQUIST_OSHER_SCALAR_2ND_ORDER=ENGQUIST_OSHER_SCALAR_2ND_ORDER_TYPE();

    class ENGQUIST_OSHER_VECTOR_1ST_ORDER_TYPE {};
    const ENGQUIST_OSHER_VECTOR_1ST_ORDER_TYPE ENGQUIST_OSHER_VECTOR_1ST_ORDER=ENGQUIST_OSHER_VECTOR_1ST_ORDER_TYPE();

    class ENGQUIST_OSHER_VECTOR_2ND_ORDER_TYPE {};
    const ENGQUIST_OSHER_VECTOR_2ND_ORDER_TYPE ENGQUIST_OSHER_VECTOR_2ND_ORDER=ENGQUIST_OSHER_VECTOR_2ND_ORDER_TYPE();

    class ENGQUIST_OSHER_SV_1ST_ORDER_TYPE {};
    const ENGQUIST_OSHER_SV_1ST_ORDER_TYPE ENGQUIST_OSHER_SV_1ST_ORDER=ENGQUIST_OSHER_SV_1ST_ORDER_TYPE();

    class ENGQUIST_OSHER_SV_2ND_ORDER_TYPE {};
    const ENGQUIST_OSHER_SV_2ND_ORDER_TYPE ENGQUIST_OSHER_SV_2ND_ORDER=ENGQUIST_OSHER_SV_2ND_ORDER_TYPE();

    class LAX_FRIEDRICHS_SCALAR_1ST_ORDER_TYPE {
    public:
        const double alpha;
        LAX_FRIEDRICHS_SCALAR_1ST_ORDER_TYPE(double a) : alpha(a) {}
    };

    LAX_FRIEDRICHS_SCALAR_1ST_ORDER_TYPE LAX_FRIEDRICHS_SCALAR_1ST_ORDER(double alpha) {
        return LAX_FRIEDRICHS_SCALAR_1ST_ORDER_TYPE(alpha);
    }

    //at
    class LAX_FRIEDRICHS_SCALAR_2ND_ORDER_TYPE {
    public:
        const double alpha;
        LAX_FRIEDRICHS_SCALAR_2ND_ORDER_TYPE(double a) : alpha(a) {}
    };

    LAX_FRIEDRICHS_SCALAR_2ND_ORDER_TYPE LAX_FRIEDRICHS_SCALAR_2ND_ORDER(double alpha) {
        return LAX_FRIEDRICHS_SCALAR_2ND_ORDER_TYPE(alpha);
    }

    //at

    class STENCIL_LOCAL_LAX_FRIEDRICHS_SCALAR_TYPE {
    public:
        const double gamma;
        STENCIL_LOCAL_LAX_FRIEDRICHS_SCALAR_TYPE(double a) : gamma(a) {}
    };

    STENCIL_LOCAL_LAX_FRIEDRICHS_SCALAR_TYPE STENCIL_LOCAL_LAX_FRIEDRICHS(double gamma) {
        return STENCIL_LOCAL_LAX_FRIEDRICHS_SCALAR_TYPE(gamma);
    }


    class SMOOTHING_SCHEME_TYPE {
    public:
        const int material_level;
        const double max_curvature;
        const double min_curvature;

        SMOOTHING_SCHEME_TYPE(int mat, double max, double min) : material_level(mat), max_curvature(max), min_curvature(min) {}
    };

    SMOOTHING_SCHEME_TYPE SMOOTHING_SCHEME(int material_level, double max_curvature, double min_curvature) {
        return SMOOTHING_SCHEME_TYPE(material_level, max_curvature, min_curvature);
    }


    template <class LevelSetType, class VelocityType, int order>
    class EngquistOsherScalar;

    template <class LevelSetType, class VelocityType, int order>
    class EngquistOsherVector;

    template <class LevelSetType, class VelocityType, int order>
    class EngquistOsherSV;

    template <class LevelSetType, class VelocityType, int order>
    class LaxFriedrichsScalar;

    template <class LevelSetType, class VelocityType, int order>
    class StencilLocalLaxFriedrichsScalar;



    template <class LevelSetType>
    class SmoothingScheme;


    template<class LevelSetType, class VelocityType, class SchemeType> class IntegrationScheme {};



    template <class LevelSetType, class VelocityType>
    class IntegrationScheme<LevelSetType, VelocityType, ENGQUIST_OSHER_SCALAR_1ST_ORDER_TYPE>:public EngquistOsherScalar<LevelSetType, VelocityType, 1> {
    public:
        IntegrationScheme(LevelSetType& l, const VelocityType& v, const ENGQUIST_OSHER_SCALAR_1ST_ORDER_TYPE& s):EngquistOsherScalar<LevelSetType, VelocityType, 1>(l,v) {}
    };

    template <class LevelSetType, class VelocityType>
    class IntegrationScheme<LevelSetType, VelocityType, ENGQUIST_OSHER_SCALAR_2ND_ORDER_TYPE>:public EngquistOsherScalar<LevelSetType, VelocityType, 2> {
    public:
        IntegrationScheme(LevelSetType& l, const VelocityType& v, const ENGQUIST_OSHER_SCALAR_2ND_ORDER_TYPE& s):EngquistOsherScalar<LevelSetType, VelocityType, 2>(l,v) {}
    };

    template <class LevelSetType, class VelocityType>
    class IntegrationScheme<LevelSetType, VelocityType, ENGQUIST_OSHER_VECTOR_1ST_ORDER_TYPE>:public EngquistOsherVector<LevelSetType, VelocityType, 1> {
    public:
        IntegrationScheme(LevelSetType& l, const VelocityType& v, const ENGQUIST_OSHER_VECTOR_1ST_ORDER_TYPE& s):EngquistOsherVector<LevelSetType, VelocityType, 1>(l,v) {}
    };

    template <class LevelSetType, class VelocityType>
    class IntegrationScheme<LevelSetType, VelocityType, ENGQUIST_OSHER_VECTOR_2ND_ORDER_TYPE>:public EngquistOsherVector<LevelSetType, VelocityType, 2> {
    public:
        IntegrationScheme(LevelSetType& l, const VelocityType& v, const ENGQUIST_OSHER_VECTOR_2ND_ORDER_TYPE& s):EngquistOsherVector<LevelSetType, VelocityType, 2>(l,v) {}
    };

    template <class LevelSetType, class VelocityType>
    class IntegrationScheme<LevelSetType, VelocityType, ENGQUIST_OSHER_SV_1ST_ORDER_TYPE>:public EngquistOsherSV<LevelSetType, VelocityType, 1> {
    public:
        IntegrationScheme(LevelSetType& l, const VelocityType& v, const ENGQUIST_OSHER_SV_1ST_ORDER_TYPE& s):EngquistOsherSV<LevelSetType, VelocityType, 1>(l,v) {}
    };

    template <class LevelSetType, class VelocityType>
    class IntegrationScheme<LevelSetType, VelocityType, ENGQUIST_OSHER_SV_2ND_ORDER_TYPE>:public EngquistOsherSV<LevelSetType, VelocityType, 2> {
    public:
        IntegrationScheme(LevelSetType& l, const VelocityType& v, const ENGQUIST_OSHER_SV_2ND_ORDER_TYPE& s):EngquistOsherSV<LevelSetType, VelocityType, 2>(l,v) {}
    };

    template <class LevelSetType, class VelocityType>
    class IntegrationScheme<LevelSetType, VelocityType, LAX_FRIEDRICHS_SCALAR_1ST_ORDER_TYPE>:public LaxFriedrichsScalar<LevelSetType, VelocityType, 1> {
    public:
        IntegrationScheme(LevelSetType& l, const VelocityType& v, const LAX_FRIEDRICHS_SCALAR_1ST_ORDER_TYPE& s):LaxFriedrichsScalar<LevelSetType, VelocityType, 1>(l,v,s.alpha) {}
    };

    //at
    template <class LevelSetType, class VelocityType>
    class IntegrationScheme<LevelSetType, VelocityType, LAX_FRIEDRICHS_SCALAR_2ND_ORDER_TYPE>:public LaxFriedrichsScalar<LevelSetType, VelocityType, 2> {
    public:
        IntegrationScheme(LevelSetType& l, const VelocityType& v, const LAX_FRIEDRICHS_SCALAR_2ND_ORDER_TYPE& s):LaxFriedrichsScalar<LevelSetType, VelocityType, 2>(l,v,s.alpha) {}
    };

    //at
    template <class LevelSetType, class VelocityType>
    class IntegrationScheme<LevelSetType, VelocityType, STENCIL_LOCAL_LAX_FRIEDRICHS_SCALAR_TYPE>:public StencilLocalLaxFriedrichsScalar<LevelSetType, VelocityType, 1> {
    public:
        IntegrationScheme(LevelSetType& l, const VelocityType& v, const STENCIL_LOCAL_LAX_FRIEDRICHS_SCALAR_TYPE& s):StencilLocalLaxFriedrichsScalar<LevelSetType, VelocityType, 1>(l,v,s.gamma) {}
    };


    template <class LevelSetType, class VelocityType>
    class IntegrationScheme<LevelSetType, VelocityType, SMOOTHING_SCHEME_TYPE>:public SmoothingScheme<LevelSetType> {
    public:
        IntegrationScheme(LevelSetType& l, const VelocityType& v, const SMOOTHING_SCHEME_TYPE& s):SmoothingScheme<LevelSetType>(l,s.material_level, s.max_curvature, s.min_curvature) {}
    };




    template <class LevelSetType, class VelocityType, int order>
    class EngquistOsherScalar {

        std::vector<typename LevelSetType::const_iterator_runs_offset> it_neighbors;       //the neighbor iterators

        const LevelSetType & LS;

        const VelocityType& velocities;

        typedef typename LevelSetType::value_type value_type;
        typedef typename LevelSetType::index_type index_type;

        bool initialized;

    public:

        static void prepare_surface_levelset(LevelSetType& l) {
            assert((order==1) || (order==2));                   //the user in the level-set-traits-class

            l.expand(order*2+1);                         //expand the level set function to ensure that for all active grid points
                                                        //the level set values of the neighbor grid points,
                                                        //which are necessary to calculate the derivatives are also defined

        }

        EngquistOsherScalar(const LevelSetType& l, const VelocityType& v): LS(l), velocities(v),initialized(false) {}


        template <class IteratorType>
        value_type operator()(const IteratorType& it, unsigned int material) {

            assert(it.is_active());

            const int D=LevelSetType::dimensions;

            if (initialized) {

                for (int i=0;i<2*D*order;i++) it_neighbors[i].go_to_indices_sequential(it.start_indices());

            } else {
                for (int i=0;i<2*D;i++) {
          vec<index_type,D> tv(index_type(0));
          for (int j=0;j<order;j++) {
            if (i<D) tv[i]++; else tv[i-D]--;
            it_neighbors.push_back(typename LevelSetType::const_iterator_runs_offset(LS, tv,it.start_indices()));
          }
          initialized=true;
        }
            }

            value_type grad_pos=0.;
            value_type grad_neg=0.;

            for (int i=0;i<D;i++) {

                const value_type pos    =   LS.grid().grid_position_of_local_index(i,it.start_indices(i));

                const value_type d_p    =   math::abs(LS.grid().grid_position_of_global_index(i,it.start_indices(i)+1)-pos);
                const value_type d_n    =   -math::abs(LS.grid().grid_position_of_global_index(i,it.start_indices(i)-1)-pos);

                const value_type phi_0=it.value();
                const value_type phi_p=it_neighbors[i*order].value();
                const value_type phi_n=it_neighbors[(i+D)*order].value();

                value_type f_p= (phi_p-phi_0)/d_p;
                value_type f_n= (phi_n-phi_0)/d_n;


                if (order==2) {         //if second order time integration scheme is used

                    const value_type d_pp   =   math::abs(LS.grid().grid_position_of_global_index(i,it.start_indices(i)+2)-pos);
                    const value_type d_nn   =   -math::abs(LS.grid().grid_position_of_global_index(i,it.start_indices(i)-2)-pos);

                    const value_type phi_pp=it_neighbors[i*order+1].value();
                    const value_type phi_nn=it_neighbors[(i+D)*order+1].value();

                    const value_type f__0 = (((d_n*phi_p    -d_p    *phi_n) /   (d_p    -d_n)   +phi_0))    /(d_p   *d_n);
                    const value_type f__n = (((d_n*phi_nn   -d_nn   *phi_n) /   (d_nn   -d_n)   +phi_0))    /(d_nn  *d_n);
                    const value_type f__p = (((d_p*phi_pp   -d_pp   *phi_p) /   (d_pp   -d_p)   +phi_0))    /(d_pp  *d_p);

                    if (math::sign(f__0)==math::sign(f__p)) {
                        if (math::abs(f__p*d_p)<math::abs(f__0*d_n)) {
                            f_p-=d_p*f__p;
                        } else {
                            f_p+=d_n*f__0;
                        }
                    }

                    if (math::sign(f__0)==math::sign(f__n)) {
                        if (math::abs(f__n*d_n)<math::abs(f__0*d_p)) {
                            f_n-=d_n*f__n;
                        } else {
                            f_n+=d_p*f__0;
                        }
                    }

                }

                grad_pos+=math::pow2(std::max(f_n,value_type(0)))+math::pow2(std::min(f_p,value_type(0)));
                grad_neg+=math::pow2(std::min(f_n,value_type(0)))+math::pow2(std::max(f_p,value_type(0)));

            }

            value_type v=velocities(it.active_pt_id(), material);
            if (v>0) {
                return std::sqrt(grad_pos)*v;
            } else {
                return std::sqrt(grad_neg)*v;
            }
        }
    };





    template <class LevelSetType, class VelocityType, int order>
    class EngquistOsherVector {

        std::vector<typename LevelSetType::const_iterator_runs_offset> it_neighbors;       //the neighbor iterators

        const LevelSetType & LS;

        const VelocityType& velocities;

        typedef typename LevelSetType::value_type value_type;
        typedef typename LevelSetType::index_type index_type;

        bool initialized;

    public:

        static void prepare_surface_levelset(LevelSetType& l) {
            assert((order==1) || (order==2));                   //the user in the level-set-traits-class

            l.expand(order*2+1);                         //expand the level set function to ensure that for all active grid points
                                                        //the level set values of the neighbor grid points,
                                                        //which are necessary to calculate the derivatives are also defined

        }

        EngquistOsherVector(const LevelSetType& l, const VelocityType& v): LS(l), velocities(v),initialized(false) {}

        template <class IteratorType>
        value_type operator()(const IteratorType& it, unsigned int material) {

            assert(it.is_active());

            const int D=LevelSetType::dimensions;

            if (initialized) {
         for (int i=0;i<2*D*order;i++) it_neighbors[i].go_to_indices_sequential(it.start_indices());
            } else {
              for (int i=0;i<2*D;i++) {
          vec<index_type,D> tv(index_type(0));
          for (int j=0;j<order;j++) {
            if (i<D) tv[i]++; else tv[i-D]--;
            it_neighbors.push_back(typename LevelSetType::const_iterator_runs_offset(LS, tv,it.start_indices()));
          }
          initialized=true;
         }
       }

            value_type grad_pos[D]={0};
            value_type grad_neg[D]={0};

            for (int i=0;i<D;i++) {

                const value_type pos    =   LS.grid().grid_position_of_local_index(i,it.start_indices(i));

                const value_type d_p    =   math::abs(LS.grid().grid_position_of_global_index(i,it.start_indices(i)+1)-pos);
                const value_type d_n    =   -math::abs(LS.grid().grid_position_of_global_index(i,it.start_indices(i)-1)-pos);

                const value_type phi_0=it.value();
                const value_type phi_p=it_neighbors[i*order].value();
                const value_type phi_n=it_neighbors[(i+D)*order].value();

                value_type f_p= (phi_p-phi_0)/d_p;
                value_type f_n= (phi_n-phi_0)/d_n;

                if (order==2) {         //if second order time integration scheme is used

                    const value_type d_pp   =   math::abs(LS.grid().grid_position_of_global_index(i,it.start_indices(i)+2)-pos);
                    const value_type d_nn   =   -math::abs(LS.grid().grid_position_of_global_index(i,it.start_indices(i)-2)-pos);

                    const value_type phi_pp=it_neighbors[i*order+1].value();
                    const value_type phi_nn=it_neighbors[(i+D)*order+1].value();

                    const value_type f__0 = (((d_n*phi_p    -d_p    *phi_n) /   (d_p    -d_n)   +phi_0))    /(d_p   *d_n);
                    const value_type f__n = (((d_n*phi_nn   -d_nn   *phi_n) /   (d_nn   -d_n)   +phi_0))    /(d_nn  *d_n);
                    const value_type f__p = (((d_p*phi_pp   -d_pp   *phi_p) /   (d_pp   -d_p)   +phi_0))    /(d_pp  *d_p);

                    /*if (math::sign(f__0)==math::sign(f__p)) {
                        if (math::abs(f__p)<math::abs(f__0)) {
                            f_p-=d_p*f__p;
                        } else {
                            f_p+=d_n*f__0;
                        }
                    }

                    if (math::sign(f__0)==math::sign(f__n)) {
                        if (math::abs(f__n)<math::abs(f__0)) {
                            f_n-=d_n*f__n;
                        } else {
                            f_n+=d_p*f__0;
                        }
                    }*/

                    if (math::sign(f__0)==math::sign(f__p)) {
                        if (math::abs(f__p*d_p)<math::abs(f__0*d_n)) {
                            f_p-=d_p*f__p;
                        } else {
                            f_p+=d_n*f__0;
                        }
                    }

                    if (math::sign(f__0)==math::sign(f__n)) {
                        if (math::abs(f__n*d_n)<math::abs(f__0*d_p)) {
                            f_n-=d_n*f__n;
                        } else {
                            f_n+=d_p*f__0;
                        }
                    }
                }

                grad_pos[i]=f_n;
                grad_neg[i]=f_p;

            }

            vec<value_type,D> v=velocities(it.active_pt_id(), material);
            value_type vel_grad=0;
            for (int w=0;w<D;w++) {
                if (v[w]>0.) {
                    vel_grad+=v[w]*grad_pos[w];
                } else {
                    vel_grad+=v[w]*grad_neg[w];
                }
            }

            return vel_grad;
        }


    };



    template <class LevelSetType, class VelocityType, int order>
    class EngquistOsherSV {

        std::vector<typename LevelSetType::const_iterator_runs_offset> it_neighbors;       //the neighbor iterators

        const LevelSetType & LS;

        const VelocityType& velocities;

        typedef typename LevelSetType::value_type value_type;
        typedef typename LevelSetType::index_type index_type;

        bool initialized;

    public:

        static void prepare_surface_levelset(LevelSetType& l) {
            assert((order==1) || (order==2));                   //the user in the level-set-traits-class

            l.expand(order*2+1);                         //expand the level set function to ensure that for all active grid points
                                                        //the level set values of the neighbor grid points,
                                                        //which are necessary to calculate the derivatives are also defined
        }

        EngquistOsherSV(const LevelSetType& l, const VelocityType& v): LS(l), velocities(v),initialized(false){}


        template <class IteratorType>
//        value_type operator()(const IteratorType& it, unsigned int material, int iterator_int) {
        value_type operator()(const IteratorType& it, unsigned int material) {

            assert(it.is_active());

            const int D=LevelSetType::dimensions;


            if (initialized) {
         for (int i=0;i<2*D*order;i++) it_neighbors[i].go_to_indices_sequential(it.start_indices());
            } else {
              for (int i=0;i<2*D;i++) {
          vec<index_type,D> tv(index_type(0));
          for (int j=0;j<order;j++) {
            if (i<D) tv[i]++; else tv[i-D]--;
            it_neighbors.push_back(typename LevelSetType::const_iterator_runs_offset(LS, tv,it.start_indices()));
          }
          initialized=true;
         }
       }


            value_type grad_pos[D];
            value_type grad_neg[D];

            value_type grad_pos_total=0;
            value_type grad_neg_total=0;


            for (int i=0;i<D;i++) {

                const value_type pos    =   LS.grid().grid_position_of_local_index(i,it.start_indices(i));

                const value_type d_p    =   math::abs(LS.grid().grid_position_of_global_index(i,it.start_indices(i)+1)-pos);
                const value_type d_n    =   -math::abs(LS.grid().grid_position_of_global_index(i,it.start_indices(i)-1)-pos);

                const value_type phi_0=it.value();
                const value_type phi_p=it_neighbors[i*order].value();
                const value_type phi_n=it_neighbors[(i+D)*order].value();

                value_type f_p= (phi_p-phi_0)/d_p;
                value_type f_n= (phi_n-phi_0)/d_n;

                if (order==2) {         //if second order time integration scheme is used

                    const value_type d_pp   =   math::abs(LS.grid().grid_position_of_global_index(i,it.start_indices(i)+2)-pos);
                    const value_type d_nn   =   -math::abs(LS.grid().grid_position_of_global_index(i,it.start_indices(i)-2)-pos);

                    const value_type phi_pp=it_neighbors[i*order+1].value();
                    const value_type phi_nn=it_neighbors[(i+D)*order+1].value();

                    const value_type f__0 = (((d_n*phi_p    -d_p    *phi_n) /   (d_p    -d_n)   +phi_0))    /(d_p   *d_n);
                    const value_type f__n = (((d_n*phi_nn   -d_nn   *phi_n) /   (d_nn   -d_n)   +phi_0))    /(d_nn  *d_n);
                    const value_type f__p = (((d_p*phi_pp   -d_pp   *phi_p) /   (d_pp   -d_p)   +phi_0))    /(d_pp  *d_p);

                    if (math::sign(f__0)==math::sign(f__p)) {
                        if (math::abs(f__p*d_p)<math::abs(f__0*d_n)) {
                            f_p-=d_p*f__p;
                        } else {
                            f_p+=d_n*f__0;
                        }
                    }

                    if (math::sign(f__0)==math::sign(f__n)) {
                        if (math::abs(f__n*d_n)<math::abs(f__0*d_p)) {
                            f_n-=d_n*f__n;
                        } else {
                            f_n+=d_p*f__0;
                        }
                    }

                }

                grad_pos[i]=f_n;
                grad_neg[i]=f_p;

                grad_pos_total+=math::pow2(std::max(f_n,value_type(0)))+math::pow2(std::min(f_p,value_type(0)));
                grad_neg_total+=math::pow2(std::min(f_n,value_type(0)))+math::pow2(std::max(f_p,value_type(0)));

            }

            value_type vel_grad=0.;

            value_type v_sca=0;
            velocities.scalar_velocity(v_sca, it.active_pt_id(), material);

            if (v_sca>0) {
              vel_grad+=std::sqrt(grad_pos_total)*v_sca;
            } else {
              vel_grad+=std::sqrt(grad_neg_total)*v_sca;
            }

            value_type v_vec[D]={0};

            velocities.vector_velocity(v_vec,it.active_pt_id(), it.start_indices(0), material);

            for (int w=0;w<D;w++) {
        if (v_vec[w]>0.) {
          vel_grad+=v_vec[w]*grad_pos[w];
        } else {
          vel_grad+=v_vec[w]*grad_neg[w];
        }
      }

      return vel_grad;

        }
    };



    template <class LevelSetType, class VelocityType, int order>
    class LaxFriedrichsScalar {


        std::vector<typename LevelSetType::const_iterator_runs_offset> it_neighbors;       // derivative neighbor iterators, relative to central point

        const LevelSetType & LS;

        const VelocityType& velocities;

        typedef typename LevelSetType::value_type value_type;
        typedef typename LevelSetType::index_type index_type;

        const double alpha;

        bool initialized;

    public:

        static void prepare_surface_levelset(LevelSetType& l) {
            assert((order==1) || (order==2));                   //the user in the level-set-traits-class

            l.expand(order*2+1);                         //expand the level set function to ensure that for all active grid points
                                                        //the level set values of the neighbor grid points,
                                                        //which are necessary to calculate the derivatives are also defined
        }


        LaxFriedrichsScalar(const LevelSetType& l, const VelocityType& v, double a): LS(l), velocities(v), alpha(a),initialized(false) {}


        template <class IteratorType>
        value_type operator()(const IteratorType& it, unsigned int material) {

          assert(it.is_active());

          const int D=LevelSetType::dimensions;

          if (initialized) {
                for (int i=0;i<2*D*order;i++)
                  it_neighbors[i].go_to_indices_sequential(it.start_indices());
          } else {
            for (int i=0;i<2*D;i++) {
              vec<index_type,D> tv(index_type(0));

              for (int j=0;j<order;j++) {
                if ( i<D ){

                  tv[i]++;
                }
                else {
                  //tv[i-D]=index_type(1);
                  tv[i-D]--;
                }
                it_neighbors.push_back(typename LevelSetType::const_iterator_runs_offset(LS, tv,it.start_indices()));

              }
            initialized=true;
            }
          }

            value_type grad=0.;
            value_type dissipation=0.;

            for (int i=0;i<D;i++) { //iterate over dimensions

                const value_type pos    =   LS.grid().grid_position_of_local_index(i,it.start_indices(i));

                const value_type d_p    =   math::abs( LS.grid().grid_position_of_global_index(i,it.start_indices(i)+1)  - pos );
                const value_type d_n    =   -math::abs( LS.grid().grid_position_of_global_index(i,it.start_indices(i)-1) - pos );

                const value_type phi_0=it.value();
                const value_type phi_p=it_neighbors[i*order].value();
                const value_type phi_n=it_neighbors[(i+D)*order].value();

                value_type f_p= (phi_p-phi_0)/d_p;
                value_type f_n= (phi_n-phi_0)/d_n;

                const value_type f__0 = (((d_n*phi_p    -d_p    *phi_n) /   (d_p    -d_n)   +phi_0))    /(d_p   *d_n);

                if (order==2) {         //if second order time integration scheme is used


                    const value_type d_pp   =   math::abs(LS.grid().grid_position_of_global_index(i,it.start_indices(i)+2)-pos);
                    const value_type d_nn   =   -math::abs(LS.grid().grid_position_of_global_index(i,it.start_indices(i)-2)-pos);

                  //  //std::cout << "d_p, d_n, d_pp, d_nn = " << d_p << ", " << d_n << ", " << d_pp << ", " <<d_nn << ", " << std::endl;

                    const value_type phi_pp=it_neighbors[i*order+1].value();
                    const value_type phi_nn=it_neighbors[(i+D)*order+1].value();


                    const value_type f__n = (((d_n*phi_nn   -d_nn   *phi_n) /   (d_nn   -d_n)   +phi_0))    /(d_nn  *d_n);
                    const value_type f__p = (((d_p*phi_pp   -d_pp   *phi_p) /   (d_pp   -d_p)   +phi_0))    /(d_pp  *d_p);

                    if (math::sign(f__0)==math::sign(f__p)) {
                        if (math::abs(f__p*d_p)<math::abs(f__0*d_n)) {
                            f_p-=d_p*f__p;
                        } else {
                            f_p+=d_n*f__0;
                        }
                    }

                    if (math::sign(f__0)==math::sign(f__n)) {
                        if (math::abs(f__n*d_n)<math::abs(f__0*d_p)) {
                            f_n-=d_n*f__n;
                        } else {
                            f_n+=d_p*f__0;
                        }
                    }

                }

                grad+=math::pow2((f_n+f_p)*0.5);
                dissipation+=(f_p-f_n)*0.5*alpha;
            }

            value_type v=velocities(it.active_pt_id(), material);
            return v*std::sqrt(grad)-((v==0.)?0:dissipation);                   //TODO

        }
    };

    //at
    template <class LevelSetType, class VelocityType, int order>
    class StencilLocalLaxFriedrichsScalar {

        //std::vector<typename LevelSetType::const_iterator_runs_offset> it_slf_stencil_points;
        std::vector<typename LevelSetType::const_iterator_runs_offset> it_neighbors;       //the neighbor iterators
        const LevelSetType & LS;
        const VelocityType& velocities;

        typedef typename LevelSetType::value_type value_type;
        typedef typename LevelSetType::index_type index_type;

        const double gamma;



        //TODO at: hard coded just for testing
        vec<value_type,3> direction100{0,1,0};
        vec<value_type,3> direction010{1,0,-1};
        const value_type r100=0.0166;
        const value_type r110=0.0309;
        const value_type r111=0.00012;
        const value_type r311=0.030;


        bool initialized;

    public:

        static void prepare_surface_levelset(LevelSetType& l) {
            assert((order==1) || (order==2) || (order==3));                   //the user in the level-set-traits-class

            //TODO sparse field expansion must depend on slf stencil order!
            l.expand(11);
            // l.expand(order*2+1);                         //expand the level set function to ensure that for all active grid points
                                                        //the level set values of the neighbor grid points,
                                                        //which are necessary to calculate the derivatives are also defined
        }


        StencilLocalLaxFriedrichsScalar(const LevelSetType& l, const VelocityType& v, double a): LS(l), velocities(v), gamma(a),initialized(false) {}


        template <class IteratorType>
        value_type operator()(const IteratorType& it, unsigned int material) {

          assert(it.is_active());




          value_type v = velocities(it.active_pt_id(), material);

          if( v == value_type(0) ){
            return 0;
          }
          else{
            const int slf_order = 1;
            const int D=LevelSetType::dimensions;

            typename LevelSetType::neighbor_stencil n(LS,it, slf_order);

            std::vector< typename LevelSetType::star_stencil> stars = n.star_stencils(order);
            typename LevelSetType::star_stencil center = stars[n.get_center_index() ];


            value_type numhamiltonian=0.; //numerical Hamiltonian
            value_type dissipation=0.;

            numhamiltonian = NormL2(center.gradient()); //|grad(phi)|
            numhamiltonian *= velocities(it.active_pt_id(), material); //V |grad(phi)|

            std::cout << numhamiltonian << std::endl;

            std::vector< vec<value_type,D>> alphas( stars.size() );
            exit(0);



          }









            /*
            bool showLog=false;


            vec<value_type,3> NormalVector;

                  slf_gradient += math::pow2( (slf_dphi_n[i] + slf_dphi_p[i])*0.5); //math::pow2( slf_dphi_n[i]);
                  slf_dphi[i] = (slf_phi_n - slf_phi_p) / (2.0 * dx); //central difference for normal

                  slf_coord[i] = slf_pos;
                  slf_phi[i] = slf_phi_0;

                  ////std::cout << LS.grid().grid_position_of_global_index(i,it_slf_stencil_points[i_slf].start_indices(i) + slf_offset[i_slf][i]) << std::endl;
              }

              if(showLog){
                std::cout << "indices = " << it.start_indices() << ", offset = " << slf_offset[i_slf]<< ", dx = " << dx << ", pos = " << slf_coord << ", phi = " << slf_phi[0];// <<", slf_normal = " << slf_normal <<", slf_gradient = " << slf_gradient;

              }


              slf_normal=slf_dphi;
              slf_gradient = std::sqrt(slf_gradient);

              //normalize normal vectors
              for(int i = 0; i < D; ++i){
                if(std::fabs(slf_gradient) < 1e-4){
                  std::cout << slf_gradient;
                  slf_normal[i] /= math::abs(slf_normal.element_max());
                }else
                  slf_normal[i] /= slf_gradient;
              }


              //determine  velocity derivatives
              for(int i=0; i < D; ++i){

                vec<value_type,3> normal_p,normal_n;

                //now it gets really dirty: force 3D vector for fourRateInterpolation
                if(D == 3){
                  normal_p = slf_normal;
                  normal_n = slf_normal;
                } else if(D==2){
                  normal_p[0] = slf_normal[0];
                  normal_n[0] = slf_normal[0];
                  normal_p[1] = slf_normal[1];
                  normal_n[1] = slf_normal[1];
                  normal_p[2] = value_type(0);
                  normal_n[2] = value_type(0);
                } else{
                  std::cerr << "Are you really simulating a 1D problem????\n";
                  return 0;
                }

                normal_p[i] -= DN;
                normal_n[i] += DN;

                value_type vp = my::math::fourRateInterpolation<value_type>(normal_p, direction100, direction010, r100, r110, r111, r311);
                value_type vn = my::math::fourRateInterpolation<value_type>(normal_n, direction100, direction010, r100, r110, r111, r311);
                if(showLog)
                  if(std::fabs(slf_normal[1] - 1.0) > 1e-4  )
                    std::cout << "vp, vn = " << vp <<", " << vn << std::endl;

                //central difference
                slf_dv[i] = (vn - vp) / (2.0 * DN);

                normal_p[i] += DN;
                normal_n[i] -= DN;
              }

              //determine \partial H / \partial phi_l
              //NOTE  check d phi!!!!
              for (int i = 0 ; i < D; ++i) { //iterate over dimensions
                //Monti term

                value_type monti = 0;
                if(1){
                  for(int j = 0; j < D - 1; ++j ){ //phi_p**2 + phi_q**2
                       int idx = (i + 1) % D;
                        monti += slf_dphi[idx] * slf_dphi[idx];
                  }
                  monti  *= slf_dv[i] / (slf_gradient * slf_gradient);
                }
                //Toifl Quell term

                value_type toifl=0;
                if(1){
                  for(int j= 0; j < D - 1; ++j ){
                     int idx = (i + 1) % D;
                     toifl += slf_dphi[idx] * slf_dv[idx];
                  }
                toifl *= -slf_dphi[i] / (slf_gradient * slf_gradient);
                }
                //Osher (constant V) term
                value_type osher=0;
                if(1)
                  osher=velocities(it.active_pt_id(), material) * slf_normal[i];

                //Total derivative is sum of terms given above
                alphas[i_slf][i] = gamma * (std::fabs(monti + toifl + osher) );

              }
              if(showLog)
                std::cout  << ", slf_dv = " << slf_dv << ", alpha_local = " << alphas[i_slf] << std::endl;

            }
              if(showLog)
                std::cout << "max alphas = [";
            for(int d=0; d < D; ++d)
            {
                std::vector<value_type> alpha_comp(num_stencil_points);

                for(int i = 0 ; i < num_stencil_points ; ++i){
                  alpha_comp[i] = alphas[i][d];
                }
                value_type maxal = *std::max_element(alpha_comp.begin(),alpha_comp.end());
                if(showLog)
                  std::cout << maxal << ", ";

                dissipation += maxal * (dphi_n[d]-dphi_p[d]) * 0.5;
            }
            if(showLog)
              std::cout << "]\n\n";

            return numhamiltonian - dissipation;

        }*/

      }
    };


    template <class LevelSetType>
    class SmoothingScheme {

        std::vector<typename LevelSetType::const_iterator_runs_offset> it_neighbors;       //the neighbor iterators

        const LevelSetType & LS;


        typedef typename LevelSetType::value_type value_type;
        typedef typename LevelSetType::index_type index_type;

        const int material_level;
        const double max_curvature;
        const double min_curvature;


        static const int order=1;

        bool initialized;

    public:

        SmoothingScheme(LevelSetType& l, int mat, double max, double min): LS(l), material_level(mat), max_curvature(max), min_curvature(min),initialized(false) {}

        static void prepare_surface_levelset(LevelSetType& l) {

            l.expand(5);                        //expand the level set function to ensure that for all active grid points
                                                //the level set values of the neighbor grid points,
                                                //which are necessary to calculate the derivatives are also defined
        }


        template <class IteratorType>
        value_type operator()(const IteratorType& it, int material) {

            assert(it.is_active());

            const int D=LevelSetType::dimensions;

            assert(D==3);

            if (initialized) {
              for (unsigned int i=0;i<it_neighbors.size();i++) it_neighbors[i].go_to_indices_sequential(it.start_indices());
            } else {
              for (int i=-1;i<=1;++i) {
          for (int j=-1;j<=1;++j) {
            for (int k=-1;k<=1;++k) {
              if (((i!=0) || (j!=0) || (k!=0)) && ((i==0) || (j==0) || (k==0))) {
                vec<index_type,D> v(i,j,k);
                it_neighbors.push_back(typename LevelSetType::const_iterator_runs_offset(LS, v,it.start_indices()));
              }
            }
          }
        }
              initialized=true;
            }

            const int XmYmZ0=0;
            const int XmY0Zm=1;
            const int XmY0Z0=2;
            const int XmY0Zp=3;
            const int XmYpZ0=4;
            const int X0YmZm=5;
            const int X0YmZ0=6;
            const int X0YmZp=7;
            const int X0Y0Zm=8;
            const int X0Y0Zp=9;
            const int X0YpZm=10;
            const int X0YpZ0=11;
            const int X0YpZp=12;
            const int XpYmZ0=13;
            const int XpY0Zm=14;
            const int XpY0Z0=15;
            const int XpY0Zp=16;
            const int XpYpZ0=17;

            double PhiX=(it_neighbors[XpY0Z0].value()-it_neighbors[XmY0Z0].value())*0.5;
            double PhiY=(it_neighbors[X0YpZ0].value()-it_neighbors[X0YmZ0].value())*0.5;
            double PhiZ=(it_neighbors[X0Y0Zp].value()-it_neighbors[X0Y0Zm].value())*0.5;

            double PhiXX=it_neighbors[XpY0Z0].value()+it_neighbors[XmY0Z0].value()-2*it.value();
            double PhiYY=it_neighbors[X0YpZ0].value()+it_neighbors[X0YmZ0].value()-2*it.value();
            double PhiZZ=it_neighbors[X0Y0Zp].value()+it_neighbors[X0Y0Zm].value()-2*it.value();

            double PhiXY=(it_neighbors[XpYpZ0].value()+it_neighbors[XmYmZ0].value()-it_neighbors[XpYmZ0].value()-it_neighbors[XmYpZ0].value())*0.25;
            double PhiXZ=(it_neighbors[XpY0Zp].value()+it_neighbors[XmY0Zm].value()-it_neighbors[XpY0Zm].value()-it_neighbors[XmY0Zp].value())*0.25;
            double PhiYZ=(it_neighbors[X0YpZp].value()+it_neighbors[X0YmZm].value()-it_neighbors[X0YmZp].value()-it_neighbors[X0YpZm].value())*0.25;

            //const int mode=0;

            double denom=PhiX*PhiX+PhiY*PhiY+PhiZ*PhiZ;

            double num=     0.5*PhiX*PhiX*(PhiYY+PhiZZ)-PhiY*PhiZ*PhiYZ+        //mean curvature
                            0.5*PhiY*PhiY*(PhiXX+PhiZZ)-PhiX*PhiZ*PhiXZ+
                            0.5*PhiZ*PhiZ*(PhiXX+PhiYY)-PhiX*PhiY*PhiXY;

            double s=0.;

            if (material<material_level) {

                if (denom!=0) {

                    double k=num/(denom*std::sqrt(denom));

                    if ((k>max_curvature) || (k<min_curvature)) s= -num/denom;

                } else {
                    //std::cout << "warning!!!dlkajf" << std::endl;
                    if (num>0) {
                        s=-std::numeric_limits<double>::max();
                    } else {
                        s=std::numeric_limits<double>::max();
                    }
                }
            }

            return s;


            /*//double denomG=std::sqrt(denomM)*denomM;


            if (denom!=0) {
                if (mode==0) {   //mean curvature

                    double num=     0.5*PhiX*PhiX*(PhiYY+PhiZZ)-PhiY*PhiZ*PhiYZ+        //mean curvature
                                    0.5*PhiY*PhiY*(PhiXX+PhiZZ)-PhiX*PhiZ*PhiXZ+
                                    0.5*PhiZ*PhiZ*(PhiXX+PhiYY)-PhiX*PhiY*PhiXY;



                    double k=num/(std::sqrt(denom)*denom);

                    if ((k<0.5) && (k>-0.5)) return 0;

                    return -num/denom*((PhiX*PhiX+PhiZ*PhiZ)/std::sqrt(denom));



                } else if (mode==1) {

                    double num=     PhiX*PhiX*(PhiYY*PhiZZ-PhiYZ*PhiYZ)+            //Gaussian curvature
                                    PhiY*PhiY*(PhiXX*PhiZZ-PhiXZ*PhiXZ)+
                                    PhiZ*PhiZ*(PhiXX*PhiYY-PhiXY*PhiXY)+
                                    2*PhiX*PhiY*(PhiXZ*PhiYZ-PhiXY*PhiZZ)+
                                    2*PhiY*PhiZ*(PhiXY*PhiXZ-PhiYZ*PhiXX)+
                                    2*PhiX*PhiZ*(PhiXY*PhiYZ-PhiXZ*PhiYY);
                    return -num/(std::sqrt(denom)*denom);

                }

            } else {
                return 0;
            }*/

        }
    };

}


#endif /*INTEGRATION_SCHEMES_HPP_*/
