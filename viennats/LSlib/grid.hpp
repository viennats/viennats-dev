#ifndef GRID_HPP_
#define GRID_HPP_

/* =========================================================================
   Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.

                            -----------------
                 ViennaTS - The Vienna Topography Simulator
                            -----------------

   Contact:         viennats@iue.tuwien.ac.at

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


#include <memory>
#include <vector>
#include "vector.hpp"
#include "math.hpp"
#include <cmath>


namespace lvlset {

    typedef unsigned int boundary_type;

    const boundary_type SYMMETRIC_BOUNDARY=0;           //reflective boundary conditions

    const boundary_type INFINITE_BOUNDARY=1;            //no boundary, the maximum allowed index is only
                                                        //limited by the constant INF_EXTENSION

    const boundary_type PERIODIC_BOUNDARY=2;            //periodic boundary conditions

    const boundary_type POS_INFINITE_BOUNDARY=3;        //one-sided reflective boundary conditions,
                                                        //no boundary in positive direction (for positive
                                                        //indices)

    const boundary_type NEG_INFINITE_BOUNDARY=4;        //one-sided reflective boundary conditions,
                                                        //no boundary in positive direction (for negative
                                                        //indices)

    template <class GridTraitsType> class grid_type {
        //the template parameter
        //has to be a class defined by the user
        //which must contain certain constants,
        //type and member definitions, specifiying the grid as shown below:
        //
        //  class GridTraitsType {
        //  public:
        //      typedef int index_type;         //the data type used for the grid indices
        //      typedef double coord_type;      //the data type used for the coordinates
        //
        //      const static int dimensions=2;  //number of dimensions (2 or 3 allowed)
        //
        //      index_type min_index(int dir) const;    //for each grid direction dir (x=0,y=1,z=2)
        //                                              //the minimum grid index has to be returned
        //      index_type max_index(int dir) const;    //for each grid direction dir (x=0,y=1,z=2)
        //                                              //the minimum grid index has to be returned
        //
        //       coord_type grid_position(int dir, index_type Index) const;     //this function returns the coordinates
        //                                                                      //of the grid planes (3D) or grid lines (2D)
        //                                                                      //with index "Index" and
        //                                                                      //which are normal to the axis given by "dir"
        //                                                                      //in this way arbritrary rectilinear grids can be defined
        //                                                                      //IMPORTANT: this grid_position-function has to be strictly monotonic increasing or decreasing
        //                                                                      //and has to be defined for all indices. In case of infinite boundaries the function has to be defined up
        //                                                                      //to +INF_EXTENSION and/or down to -INF_EXTENSION (INF_EXTENSION is defined below)
        //
        //      boundary_type boundary_condition(int dir) const;    //for each grid direction "dir" the boundary condition type
        //                                                          //has to be specified, returning one of the boundary_type-constants
        //                                                          //defined above
        //
        //IMPORTANT: the GridTraitsType has to be copyable, make sure that you define a appropriate copy-constructor if the default copy-constructor does not work for your class

    public:
        typedef typename GridTraitsType::index_type index_type;
        typedef typename GridTraitsType::coord_type coord_type;
        //typedef GridTraitsType grid_traits_type;

        static const int dimensions=GridTraitsType::dimensions;

    private:

        static const index_type INF_EXTENSION;
        GridTraitsType GridTraits;
        static const int D=GridTraitsType::dimensions;

        vec<index_type,D> Min_,Max_,Ext_;                           //Minimum, maximum, extension of whole grid for all grid directions

        vec<boundary_type,D> BoundaryConditions_;                   //here the boundary conditions for all grid directions are stored

        vec<index_type,D> MinGridPointCoord_, MaxGridPointCoord_;   //effective maximum and minimum grid point coordinates
                                                                    //due to periodic boundary conditions the grid points at opposite boundaries are the same
                                                                    //therefore these indices are less for periodic boundary conditions
                                                                    //for periodic boundary conditions in direction k, MaxGridPointCoord_[k]=Max_[k]-1 holds

        vec<bool, D> parities;                                      //here the parities of the grid are stored, the parity is false if the grid_position function is strictly monotonic increasing
                                                                    //and true if the grid_position function is strictly monotonic decreasing

    public:

        grid_type(const GridTraitsType& gt) : GridTraits(gt) {      //constructor
            for (int i=0;i<D;i++) {
                BoundaryConditions_[i]=GridTraits.boundary_condition(i);

                MinGridPointCoord_[i]=GridTraits.min_index(i);
                MaxGridPointCoord_[i]=GridTraits.max_index(i);
                Min_[i]=GridTraits.min_index(i);
                Max_[i]=GridTraits.max_index(i);

                if (    (GridTraits.boundary_condition(i)==INFINITE_BOUNDARY) ||
                        (GridTraits.boundary_condition(i)==NEG_INFINITE_BOUNDARY)) {
                    MinGridPointCoord_[i]=-INF_EXTENSION;
                    Min_[i]=-INF_EXTENSION;
                }
                if (    (GridTraits.boundary_condition(i)==INFINITE_BOUNDARY) ||
                        (GridTraits.boundary_condition(i)==POS_INFINITE_BOUNDARY)) {
                    MaxGridPointCoord_[i]=INF_EXTENSION;
                    Max_[i]=INF_EXTENSION;
                }
                if (GridTraits.boundary_condition(i)==PERIODIC_BOUNDARY) MaxGridPointCoord_[i]--;
            }
            Ext_=Max_-Min_;
            //initialize parities
            for (int i=0;i<D;i++) {
                parities[i]=(GridTraits.grid_position(i, Min_[i]) > GridTraits.grid_position(i,Max_[i]));
            }
        }

        //copy constructor
        grid_type(const grid_type& gt): GridTraits(gt.GridTraits),
                                        Min_(gt.Min_),
                                        Max_(gt.Max_),
                                        Ext_(gt.Ext_),
                                        BoundaryConditions_(gt.BoundaryConditions_),
                                        MinGridPointCoord_(gt.MinGridPointCoord_),
                                        MaxGridPointCoord_(gt.MaxGridPointCoord_),
                                        parities(gt.parities) {}
        //empty constructor
        grid_type():GridTraits(GridTraitsType()){}

        void print() const{
          std::cout << "Min_: " << Min_ << std::endl;
          std::cout << "Max_: " << Max_ << std::endl;
          std::cout << "Ext_: " << Ext_ << std::endl;
          std::cout << "MinGridPoint: " << MinGridPointCoord_ << std::endl;
          std::cout << "MaxGridPoint: " << MaxGridPointCoord_ << std::endl;
          std::cout << "BNC: " << BoundaryConditions_ << std::endl;
          GridTraits.print();
        }

        bool parity(int dim) const {
            //parity is false/true if the "grid_position" function in GridTraitsType is
            //monotonic increasing/decreasing respectively for the given grid direction
            return parities[dim];
        }

        bool parity() const {
            //returns the total parity of the grid
            bool b=parity(0);
            for (int i=1;i<D;i++) b ^=parity(i);
            return b;
        }

        index_type grid_extent(int dim) const{
          return Ext_[dim];
        }

        index_type min_grid_index(int dim) const {
            return Min_[dim];
        }

        index_type max_grid_index(int dim) const {
            return Max_[dim];
        }

        inline const vec<index_type,D>& min_grid_index() const {
            return Min_;
        }

        inline const vec<index_type,D>& max_grid_index() const {
            return Max_;
        }

        inline const vec<boundary_type,D> & boundary_conditions() const {
            return BoundaryConditions_;
        }

        inline boundary_type boundary_conditions(int dir) const {
            return BoundaryConditions_[dir];
        }

        bool is_boundary_periodic(int dim) const {
            return BoundaryConditions_[dim]==PERIODIC_BOUNDARY;
        }

        bool is_pos_boundary_infinite(int dim) const {
            //returns true if the boundary is infinite for positive indices in "dim"-direction
            return (    (BoundaryConditions_[dim]==INFINITE_BOUNDARY) ||
                        (BoundaryConditions_[dim]==POS_INFINITE_BOUNDARY));
        }

        bool is_neg_boundary_infinite(int dim) const {
            //returns true if the boundary is infinite for negative indices in "dim"-direction
            return (    (BoundaryConditions_[dim]==INFINITE_BOUNDARY) ||
                        (BoundaryConditions_[dim]==NEG_INFINITE_BOUNDARY));
        }

        template <class V>
        bool is_cell_member(const V& vec) const {
          for (int i=0;i<D;++i) {
            if ((vec[i]<min_grid_index(i)) || (vec[i]>=max_grid_index(i))) return false;
          }
          return true;
        }

        inline index_type max_point_index(int dim) const {
            return MaxGridPointCoord_[dim];
        }

        inline index_type min_point_index(int dim) const {
            return MinGridPointCoord_[dim];
        }

        inline const vec<index_type,D>& max_point_index() const {
            return MaxGridPointCoord_;
        }

        inline const vec <index_type,D>& min_point_index() const {
            return MinGridPointCoord_;
        }

        template <class V> inline bool is_at_infinity(const V& v) const {
            //this function returns if one of the indices of
            //the index-vector "v" is infinite
            for (int i=0;i<D;i++) {
                if (math::abs(v[i])==INF_EXTENSION) return true;
            }
            return false;
        }

        inline index_type local_index_2_global_index(int dim, index_type relative_coord, int cycles, index_type offset=0) const {
            //this function transforms a local index to the corresponding global index
            //the global index is the index of an infinite grid in all directions
            //if symmetric or periodic boundary conditions are used more than one global index is mapped
            //to the same local index
            //for the local index always (min_point_index<=local index<=max_point_index holds)

            if (cycles==0) {
                return relative_coord-offset;
            } else {
                if (is_boundary_periodic(dim)) {
                    return relative_coord-offset+cycles*(Max_[dim]-Min_[dim]);
                } else {
                    if ((cycles & 1)==0) { //if cycles is even
                        return cycles*Ext_[dim]+relative_coord-offset;
                    } else {                //if cycle is odd
                        return cycles*Ext_[dim]+Max_[dim]+Min_[dim]-relative_coord-offset;
                    }
                }
            }
        }

        inline index_type global_index_2_local_index(int dim, index_type absolute_coord, index_type offset, int &cycles) const {
            //this function transforms a global index to the corresponding local index
            //the global index is the index of an infinite grid in all directions
            //if symmetric or periodic boundary conditions are used more than one global index are mapped
            //to the same local index
            //for the local index always (min_point_index<=local index<=max_point_index holds)

            absolute_coord+=offset;
            cycles=0;
            while (absolute_coord<Min_[dim]) {
                cycles--;
                absolute_coord+=Ext_[dim];
            }
            while (absolute_coord>=Max_[dim]) {
                cycles++;
                absolute_coord-=Ext_[dim];
            }
            if (((cycles & 1)==0) || (is_boundary_periodic(dim))) {
                return absolute_coord;
            } else {
                return Min_[dim]+Max_[dim]-absolute_coord;
            }
        }

        inline index_type global_index_2_local_index(int dim, index_type absolute_coord, index_type offset=0) const {
            //this function transforms a global index to the corresponding local index
            //the global index is the index of an infinite grid in all directions
            //if symmetric or periodic boundary conditions are used more than one global index are mapped
            //to the same local index
            //for the local index always (min_point_index<=local index<=max_point_index holds)

            absolute_coord+=offset;
            bool b=true;
            while (absolute_coord<Min_[dim]) {
                b=!b;
                absolute_coord+=Ext_[dim];
            }
            while (absolute_coord>=Max_[dim]) {
                b=!b;
                absolute_coord-=Ext_[dim];
            }
            if (b || is_boundary_periodic(dim)) {
                return absolute_coord;
            } else {
                return Min_[dim]+Max_[dim]-absolute_coord;
            }
        }

        template<class V>
        inline vec<index_type,D> global_indices_2_local_indices(const V& v) const {
            //this function transforms a global index vector to the corresponding local index vector

            vec<index_type,D> tmp;
            for (int i=0;i<D;i++) {
                tmp[i]=global_index_2_local_index(i,v[i]);
            }
            return tmp;

        }

        coord_type grid_position_of_local_index(int dir, index_type Index) const {
            return GridTraits.grid_position(dir, Index);
        }

        coord_type local_index_2_local_coordinate(int dir, coord_type c) const {
            //this function transforms the coordinate c in respect to the rectilinear grid into the
            //real coordinates

            coord_type lc=std::floor(c);
            coord_type uc=std::ceil(c);

            if (lc!=uc) {
                return (    grid_position_of_local_index(dir, static_cast<index_type>(lc))*((uc-c))+    //TODO
                            grid_position_of_local_index(dir, static_cast<index_type>(uc))*((c-lc))        );
            } else {
                return grid_position_of_local_index(dir,  static_cast<index_type>(lc));
            }
        }

        //[Josef] This function certainly has a role to play
        coord_type global_index_2_global_coordinate(int dir, coord_type c) const {
            //this function transforms the coordinate c in respect to the rectilinear grid into the
            //real coordinates

            coord_type lc=std::floor(c);
            coord_type uc=std::ceil(c);
            if (lc!=uc) {
                return (    grid_position_of_global_index(dir, static_cast<index_type>(lc))*((uc-c)/(uc-lc))+        //TODO
                            grid_position_of_global_index(dir, static_cast<index_type>(uc))*((c-lc)/(uc-lc))        );
            } else {
                return grid_position_of_global_index(dir,  static_cast<index_type>(lc));
            }
        }

        template<class V>
        vec<coord_type, D> global_indices_2_global_coordinates(const V& v) const {
          vec<coord_type, D> tmp;
          for(unsigned i=0; i<D; ++i) tmp[i] = global_index_2_global_coordinate(i, v[i]);
          return tmp;
        }

        index_type global_coordinate_2_global_index(int dir, coord_type c) const {
            return index_type( round(c / GridTraits.grid_delta()) );
        }

        template<class V>
        vec<index_type, D> global_coordinates_2_global_indices(const V& v) const {
          vec<index_type, D> tmp;
          for(unsigned i=0; i<D; ++i) tmp[i] = global_coordinate_2_global_index(i, v[i]);
          return tmp;
        }

        coord_type local_coordinate_2_local_index(int dir, coord_type c, index_type a, index_type b) const {            //TODO global/local coordinates
            //this function transforms the real coordinate c into
            //the coordinate in respect to the indices of the rectilinear grid
            //a and b allow to restrict the search between two grid indices

            coord_type ac=grid_position_of_local_index(dir, a);
            coord_type bc=grid_position_of_local_index(dir, b);

            if (ac>bc) {
                std::swap(ac,bc);
                std::swap(a,b);
            }

            if (c>bc) return b;
            if (c<ac) return a;

            while (math::abs(a-b)>index_type(1)) {

                index_type mid=(a+b)/2;
                coord_type midc=grid_position_of_local_index(dir, mid);

                if (c<=midc) {
                    b=mid;
                    bc=midc;
                } else {
                    a=mid;
                    ac=midc;
                }
            }

            return a+(b-a)*((c-ac)/(bc-ac));

        }

        coord_type local_coordinate_2_local_index(int dir, coord_type c) const {
            //this function transforms the real coordinate c into
            //the coordinate in respect to the indices of the rectilinear grid

            index_type a=min_grid_index(dir);
            index_type b=max_grid_index(dir);
            return local_coordinate_2_local_index(dir,c,a,b);

        }


        /*coord_type global_coordinate_2_global_index(int dim, coord_type absolute_coord) const {
            int cycles=0;
            if (is_boundary_periodic(dim)) {
                while (absolute_coord<min_local_coordinate(dim)) {
                    cycles--;
                    absolute_coord+=(max_local_coordinate(dim)-min_local_coordinate(dim));
                }
                while (absolute_coord>=max_local_coordinate(dim)) {
                    cycles++;
                    absolute_coord-=(max_local_coordinate(dim)-min_local_coordinate(dim));
                }
                if (parity(dim)) cycles=-cycles;
                return local_coordinate_2_local_index(dim, absolute_coord)+cycles*Ext_[dim];
            } else {
                do {
                    if (absolute_coord<min_local_coordinate(dim)) {
                        cycles--;
                        absolute_coord=min_local_coordinate(dim)+(min_local_coordinate(dim)-absolute_coord);
                        continue;
                    }
                    if (absolute_coord>max_local_coordinate(dim)) {
                        cycles++;
                        absolute_coord=max_local_coordinate(dim)+(max_local_coordinate(dim)-absolute_coord);
                        continue;
                    }
                } while(false);
            }
            if (parity(dim)) cycles=-cycles;

            if ((cycles & 1)==0) {
                return local_coordinate_2_local_index(dim, absolute_coord)+cycles*Ext_[dim];
            } else {
                return cycles*Ext_[dim]+Max_[dim]+Min_[dim]-local_coordinate_2_local_index(dim, absolute_coord);
            }
        }*/

        const GridTraitsType& grid_traits() const {
            return GridTraits;
        }

        const typename GridTraitsType::coord_type grid_delta() const{
            return GridTraits.grid_delta();
        }


        coord_type grid_position_of_global_index(int dir, index_type offset) const {    //TODO check
            //returns the grid position for a global index
            //while taking account the boundary conditions

            int cycles=0;

            while (offset<Min_[dir]) {
                cycles--;
                offset+=Ext_[dir];
            }
            while (offset>=Max_[dir]) {
                cycles++;
                offset-=Ext_[dir];
            }

            if (((cycles & 1)==0) || is_boundary_periodic(dir)) {
                return grid_position_of_local_index(dir, offset)+cycles*(grid_position_of_local_index(dir, max_grid_index(dir))-grid_position_of_local_index(dir, min_grid_index(dir)));
            } else {
                return (1+cycles)*grid_position_of_local_index(dir, max_grid_index(dir))+(1-cycles)*grid_position_of_local_index(dir, min_grid_index(dir))-grid_position_of_local_index(dir, max_grid_index(dir)+min_grid_index(dir)-offset);
            }
        }

        coord_type min_local_coordinate(int dir) const {
            //returns the minimum coordinate of the domain of the grid
            if (parity(dir)) {
                return grid_position_of_local_index(dir, max_grid_index(dir));
            } else {
                return grid_position_of_local_index(dir, min_grid_index(dir));
            }
        }

        coord_type max_local_coordinate(int dir) const {
            //returns the maximum coordinate of the domain of the grid
            if (parity(dir)) {
                return grid_position_of_local_index(dir, min_grid_index(dir));
            } else {
                return grid_position_of_local_index(dir, max_grid_index(dir));
            }
        }

        template<class V>
        V increment_indices(V v) const {
            //this function transforms a global index vector to the corresponding local index vector

            int dim=0;
            for (;dim<D-1;++dim) {
                if (v[dim]<max_point_index(dim)) break;
                v[dim]=min_point_index(dim);
            }
            ++v[dim];
            return v;
        }

        template<class V>
        V decrement_indices(V v) const {
            //this function transforms a global index vector to the corresponding local index vector

            int dim=0;
            for (;dim<D-1;++dim) {
                if (v[dim]>min_point_index(dim)) break;
                v[dim]=max_point_index(dim);
            }
            --v[dim];
            return v;
        }

        // determine whether index is on border of simulation domain
        template<class V>
        bool is_border_point(V v) const {
          for(unsigned i=0; i<D; ++i){
            if(v[i] <= Min_[i] || v[i] >= Max_[i]) return true;
          }
          return false;
        }

        // determine whether point is outside of domain in direction other than infinite boundary
        template<class V>
        bool is_outside_domain(V v) const{
          for(unsigned i=0; i<D; ++i){
            if(BoundaryConditions_[i]==INFINITE_BOUNDARY) continue;
            if(BoundaryConditions_[i]==PERIODIC_BOUNDARY && v[i]==Max_[i]) return true; 
            if(v[i]<Min_[i] || v[i]>Max_[i]) return true;
          }
          return false;
        }
    };

    template <class GridTraitsType> const typename grid_type<GridTraitsType>::index_type
    grid_type<GridTraitsType>::INF_EXTENSION=std::numeric_limits<typename grid_type<GridTraitsType>::index_type>::max()/3;  //in this way the addition of 2 indices
                                                                                                                            //is still defined (INF_EXTENSION+INF_EXTENSION)
}


#endif /*GRID_HPP_*/
