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

        std::vector<typename LevelSetType::const_iterator_runs_offset> it_neighbors;       //the neighbor iterators

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

            value_type grad=0.;
            value_type dissipation=0.;

            for (int i=0;i<D;i++) {

                const value_type pos    =   LS.grid().grid_position_of_local_index(i,it.start_indices(i));

                const value_type d_p    =   math::abs(LS.grid().grid_position_of_global_index(i,it.start_indices(i)+1)-pos);
                const value_type d_n    =   -math::abs(LS.grid().grid_position_of_global_index(i,it.start_indices(i)-1)-pos);

                const value_type phi_0=it.value();
                const value_type phi_p=it_neighbors[i*order].value();
                const value_type phi_n=it_neighbors[(i+D)*order].value();

                value_type f_p= (phi_p-phi_0)/d_p;
                value_type f_n= (phi_n-phi_0)/d_n;

                const value_type f__0 = (((d_n*phi_p    -d_p    *phi_n) /   (d_p    -d_n)   +phi_0))    /(d_p   *d_n);

                if (order==2) {         //if second order time integration scheme is used

                    const value_type d_pp   =   math::abs(LS.grid().grid_position_of_global_index(i,it.start_indices(i)+2)-pos);
                    const value_type d_nn   =   -math::abs(LS.grid().grid_position_of_global_index(i,it.start_indices(i)-2)-pos);

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
                    std::cout << "warning!!!dlkajf" << std::endl;
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
