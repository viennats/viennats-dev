#ifndef SURFACE2LEVELSET_HPP_
#define SURFACE2LEVELSET_HPP_

/* =========================================================================
   Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.

                            -----------------
                 ViennaTS - The Vienna Topography Simulator
                            -----------------

   Contact:         viennats@iue.tuwien.ac.at

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


#include <cassert>
#include <vector>
#include <algorithm>
#include <bitset>

#include "output.hpp"

#include "vector.hpp"
#include "kernel.hpp"
#include "math.hpp"

#include <iostream>  //TODO test

#include "misc.hpp"

namespace lvlset {

    ///this function checks if a line containing the point with coordinates "Point"
    ///and parallel to axis "dir" (x=0, y=1, z=2) intersects the surface element given by the nodes "c"
    ///additionally the normal-vector "SurfaceNormalVector" of this surface element has to be given
    ///if there is an intersection this function returns true, otherwise false
    ///the intersection coordinate is returned by "intersection"
    template <class T>
    int calculate_gridline_triangle_intersection(
                                        vec<T,3> Point,
                                        const vec<T,3> *c,
                                        int dir,
                                        T& intersection
                                    ) {



        bool inside_pos=true;
        bool inside_neg=true;

        T A[3];

        const int dirA=(dir+1)%3;
        const int dirB=(dir+2)%3;

        for (int k=0;k<3;k++) {
            /*A[k]=(c[(k+1)%3][dirA]-c[(k+2)%3][dirA])*((c[(k+2)%3][dirB]-Point[dirB])+(c[(k+1)%3][dirB]-Point[dirB]))+
                 (c[(k+2)%3][dirB]-c[(k+1)%3][dirB])*((c[(k+1)%3][dirA]-Point[dirA])+(c[(k+2)%3][dirA]-Point[dirA]));*/

          bool swapped=(c[(k+1)%3]<c[(k+2)%3]);          //necessary to guarantee anti commutativity
          const vec<T,3>& v1=(swapped)?c[(k+2)%3]:c[(k+1)%3];
          const vec<T,3>& v2=(swapped)?c[(k+1)%3]:c[(k+2)%3];

          A[k]=(v1[dirA]-Point[dirA])*(v2[dirB]-Point[dirB])-(v2[dirA]-Point[dirA])*(v1[dirB]-Point[dirB]);

          if (swapped) A[k]=-A[k];

            if (A[k]<T(0)) inside_pos=false;
            if (A[k]>T(0)) inside_neg=false;

        }

        if ((!inside_pos && inside_neg) || (inside_pos && !inside_neg)) {
            T sum=A[0]+A[1]+A[2];
            int k=0;
            for (int i=1;i<3;++i) {
              if (inside_pos) {
                if (A[i]>A[k]) k=i;
              } else {
                if (A[i]<A[k]) k=i;
              }
            }
            intersection=c[k][dir]+(c[(k+1)%3][dir]-c[k][dir])*(A[(k+1)%3]/sum)+(c[(k+2)%3][dir]-c[k][dir])*(A[(k+2)%3]/sum);
            //intersection=c[0][dir]+(A[1]/sum)*(c[1][dir]-c[0][dir])+(A[2]/sum)*(c[2][dir]-c[0][dir]);
            //intersection=(A[0]/sum)*c[0][dir]+(A[1]/sum)*c[1][dir]+(A[2]/sum)*c[2][dir];
            return (inside_pos)?1:-1;
        } else {
            return 0;
        }
    }


    ///this function checks if a line containing the point with coordinates "Point"
    ///and parallel to axis "dir" (x=0, y=1) intersects the surface element given by the nodes "c"
    ///additionally the normal-vector "SurfaceNormalVector" of this surface element has to be given
    ///if there is an intersection this function returns true, otherwise false
    ///the intersection coordinate is returned by "intersection"
    template <class T>
    int calculate_gridline_triangle_intersection(
                                        vec<T,2> Point,
                                        const vec<T,2> *c,
                                        int dir,
                                        T& intersection
                                    ) {

        bool inside_pos=true;
        bool inside_neg=true;

        T A[2];

        const int dirA=(dir+1)%2;

        for (int k=0;k<2;k++) {
            A[k]=c[(k+1)%2][dirA]-Point[dirA];
            if (k==dir) A[k]=-A[k];
            if (A[k]<T(0)) inside_pos=false;
            if (A[k]>T(0)) inside_neg=false;
        }

        if ((!inside_pos && inside_neg) || (inside_pos && !inside_neg)) {
          T sum=A[0]+A[1];
      int k=0;
      for (int i=1;i<2;++i) {
        if (inside_pos) {
          if (A[i]>A[k]) k=i;
        } else {
          if (A[i]<A[k]) k=i;
        }
      }
      intersection=c[k][dir]+(c[(k+1)%2][dir]-c[k][dir])*(A[(k+1)%2]/sum);
            //intersection=(A[0]/sum)*c[0][dir]+(A[1]/sum)*c[1][dir];
            return (inside_pos)?1:-1;
        } else {
            return 0;
        }
    }

    ///Initialise level set function
    template<class TriangulationType, class LevelSetType>
    void init(      LevelSetType& l,            //the level set function which should be initialized
                        const TriangulationType& srf,    //the triangulated surface
                        bool report_import_errors,        //to ignore errors in the imported geometry mesh
                        typename LevelSetType::value_type eps_sign=1e-6,    //"eps_sign" is used to determine
                                                    //the sign of the distance to the surface
                        typename LevelSetType::value_type eps_boundary=1e-5,  //"eps_boundary" defines the range at grid boundaries
                                                    //within which all nodes snap exactly to the boundary coordinates
                        typename LevelSetType::value_type eps_distance=1e-4,  //"eps_distance" ensures that all grid points are initialized which are
                                                    //are on a grid segment which is intersected by the surface element
                        typename LevelSetType::value_type eps_distance2=1e-7

        ) {

        typedef typename TriangulationType::element_index_type element_index_type;
//        typedef typename TriangulationType::node_index_type node_index_type;
        typedef typename LevelSetType::index_type index_type;
        typedef typename LevelSetType::value_type value_type;

        const int D=LevelSetType::dimensions;

        std::vector<std::pair<vec<index_type,D>,value_type> > points2;

        //############################################################################

        //determine boundar-epsilons for the real coordinates for all axis directions

        value_type boundary_eps_min[D];
        value_type boundary_eps_max[D];

        for (int i=0;i<D;++i) {

            if (l.grid().is_neg_boundary_infinite(i)) {
                boundary_eps_min[i]=static_cast<value_type>(0);
            } else {
                boundary_eps_min[i] = eps_boundary*l.grid().grid_delta();
            }

            if (l.grid().is_pos_boundary_infinite(i)) {
                boundary_eps_max[i]=static_cast<value_type>(0);
            } else {
                boundary_eps_max[i] = eps_boundary*l.grid().grid_delta();
            }

            if (l.grid().parity(i)) {
                std::swap(boundary_eps_min[i],boundary_eps_max[i]);
                boundary_eps_min[i]=-boundary_eps_min[i];
                boundary_eps_max[i]=-boundary_eps_max[i];
            }

            assert(boundary_eps_min[i]>=0);
            assert(boundary_eps_max[i]>=0);
        }

        //###############################################################################################

        //setup list of grid points with distances to surface elements

        {
            typedef typename std::vector< std::pair< vec<index_type,D>, std::pair<value_type,value_type> > > point_vector;
            point_vector points;

            //for each surface element do
            for (element_index_type e=0;e<srf.number_of_elements();e++) {
                vec<value_type,D> c[D];            //nodes of element
                vec<value_type,D> center(value_type(0));  //center point of triangle

                std::bitset<2*D> flags;
                flags.set();

                for (int dim=0;dim<D;dim++) {

                    for (int q=0;q<D;q++) {
                        c[q][dim]=srf.node_coordinate(srf.element_node_id(e, q), dim);
                        if (math::abs(c[q][dim]-l.grid().min_local_coordinate(dim))<boundary_eps_min[dim])
                          c[q][dim]= l.grid().min_local_coordinate(dim);
                        if (math::abs(c[q][dim]-l.grid().max_local_coordinate(dim))<boundary_eps_max[dim])
                          c[q][dim]= l.grid().max_local_coordinate(dim);

                        if (c[q][dim]>l.grid().min_local_coordinate(dim))
                          flags.reset(dim);      //TODO
                        if (c[q][dim]<l.grid().max_local_coordinate(dim))
                          flags.reset(dim+D);    //TODO

                        center[dim]+=c[q][dim];  //center point calculation
                    }
                }

                if (flags.any()){
                  continue;
                }   //triangle is outside of domain


                center/=static_cast<value_type>(D);  //center point calculation

                //vec<value_type,D> normal=NormalVector(c);  //normalvector calculation

                //determine min and max of nodes
                vec<value_type,D> min_c=c[0];
                vec<value_type,D> max_c=c[0];
                for (int i=1;i<D;++i) {
                    min_c=Min(min_c,c[i]);
                    max_c=Max(max_c,c[i]);
                }

                vec<value_type,D> min, max;
                vec<index_type,D> min_idx,max_idx;
                for (int q=0;q<D;q++) {
                    min[q]=l.grid().local_coordinate_2_local_index(q, min_c[q]);
                    max[q]=l.grid().local_coordinate_2_local_index(q, max_c[q]);
//                    if (l.grid().parity(q)) std::swap(min[q], max[q]);
                    assert(min[q]<=max[q]);
                    assert(min[q]>=l.grid().min_grid_index(q));
                    assert(max[q]<=l.grid().max_grid_index(q));
                    //min_idx[q]=std::max(static_cast<index_type>(std::ceil(min[q])),l.grid().min_grid_index(q));
                    //max_idx[q]=std::min(static_cast<index_type>(std::floor(max[q])),l.grid().max_grid_index(q));
                    min_idx[q]=static_cast<index_type>(std::ceil(min[q]));
                    max_idx[q]=static_cast<index_type>(std::floor(max[q]));
                    assert(min_idx[q]>=l.grid().min_grid_index(q));
                    assert(max_idx[q]<=l.grid().max_grid_index(q));
                }



                for (int z=0;z<D;z++) {         //for each axis direction do

                    vec<index_type,D-1> min_bb, max_bb;

                    for (int h=0;h<D-1;++h) {
                        min_bb[h]=min_idx[(z+h+1)%D];
                        max_bb[h]=max_idx[(z+h+1)%D];
                    }

                    box<index_type, D-1> bb(min_bb, max_bb);

                    if (bb.is_empty()) continue;

                    for (typename box<index_type,D-1>::iterator it_bb(bb);!it_bb.is_finished();it_bb++) {

                        vec<index_type,D> it_b;
                        for (int h=0;h<D-1;++h) it_b[(z+h+1)%D]=(*it_bb)[h];

                        vec<value_type,D> p;
                        for (int k=1;k<D;k++) p[(k+z)%D]=l.grid().grid_position_of_global_index((k+z)%D, it_b[(k+z)%D]);

                        value_type intersection;
                        int intersection_status=calculate_gridline_triangle_intersection(  p,
                                                                                            c,
                                                                                            z,
                                                                                            intersection
                                                                                        );

                        if (intersection_status!=0) {
                            //if there is an intersection

                          if (intersection<min_c[z]) assert(0);    //TODO
                          if (intersection>max_c[z]) assert(0);    //TODO
                            intersection=std::max(intersection, min_c[z]);
                            intersection=std::min(intersection, max_c[z]);

                            if (intersection>l.grid().max_local_coordinate(z)) continue;
                            if (intersection<l.grid().min_local_coordinate(z)) continue;

                            value_type intersection2=l.grid().local_coordinate_2_local_index(z, intersection);

                            index_type floor=static_cast<index_type>(std::floor(intersection2-eps_distance));
                            index_type ceil=static_cast<index_type>(std::ceil(intersection2+eps_distance));

                            floor=std::max(floor, min_idx[z]-1);
                            ceil=std::min(ceil, max_idx[z]+1);
                            floor=std::max(floor,l.grid().min_grid_index(z));
                            ceil=std::min(ceil, l.grid().max_grid_index(z));

                            vec<value_type, D> t=center;
                            t[z]-=intersection;
                            for (int k=1;k<D;k++) t[(z+k)%D]-=p[(z+k)%D];
                            t=Normalize(t);

                            for (it_b[z]=floor;it_b[z]<=ceil;++(it_b[z])) {

                                value_type RealDistance=it_b[z]-intersection2;

                                if (l.grid().parity(z)) RealDistance=-RealDistance;

                                value_type SignDistance=RealDistance;

                                SignDistance-=eps_sign*t[z];

                                //if (normal[z]<0) {
                                if (intersection_status<0) {
                                    RealDistance=-RealDistance;
                                    SignDistance=-SignDistance;
                                }

                                if (math::signbit(RealDistance)) SignDistance=-SignDistance;

                                RealDistance*=(1.-eps_distance2);
                                if (RealDistance>1.) RealDistance=1.;
                                if (RealDistance<-1.) RealDistance=-1.;

                                if (RealDistance==0.) RealDistance=0.;          //to avoid zeros with negative sign

                                points.push_back(   std::make_pair(
                                                        l.grid().global_indices_2_local_indices(it_b),
                                                        std::make_pair(SignDistance, RealDistance)
                                                    )
                                                );
                            }
                        }
                    }
                }
            }

            std::sort(points.begin(),points.end());         //sort points lexicographically

            //setup list of index/distance pairs for level set initialization
            typename point_vector::iterator it_points=points.begin();
            while(it_points!=points.end()) {
                vec<index_type,D> tmp=it_points->first;
                points2.push_back(std::make_pair(it_points->first, it_points->second.second));

                do {
                    ++it_points;
                } while ((it_points!=points.end()) && (it_points->first==tmp));
            }
        }

        l.insert_points(points2);    //initialize level set function

        if (report_import_errors) {
          std::string err= misc::test(l);     //check level set function
          if (err.size()) {                   //if inconsistent print out error message
            std::cout << "Initialization of level set function from triangulated surface failed!" << std::endl;
            #ifdef VERBOSE
            std::cout << err << std::endl;
            #endif // VERBOSE
            write_explicit_levelset(l, "importError.vtp");
            std::cout << "Level set points have been written to importError.vtp" << std::endl;
            abort();
          }
        }

        l.prune();       //remove active grid point which have no opposite signed neighbor grid point
        l.segment();    //distribute points evenly across threads
    }

}

#endif /*SURFACE2LEVELSET_HPP_*/
