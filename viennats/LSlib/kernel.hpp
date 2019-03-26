#ifndef KERNEL_HPP_
#define KERNEL_HPP_

/* =========================================================================
   Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.

                            -----------------
                 ViennaTS - The Vienna Topography Simulator
                            -----------------

   Contact:         viennats@iue.tuwien.ac.at

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


#include <vector>
#include <limits>
#include <cassert>
#include <algorithm>
#include <string>
#include <sstream>
#include <bitset>

#include <iostream>     //TODO
#include "misc.hpp"     //TODO test
#include <iomanip>

#include "vector.hpp"
#include "grid.hpp"
#include "math.hpp"

#ifdef _OPENMP
    #include <omp.h>
#endif


static const unsigned int lvlset_omp_max_num_threads=1000;


//#include "misc.hpp"

//in this file the basic functions and data structures of the level set library are defined
//as data structure hierarchical run-length-encoding (HRLE) as described in
//      B.~Houston, M.B.~Nielsen, C.~Batty, O.~Nilsson, K.~Museth,
//      "Hierarchical RLE-Level Set - A Compact and Versatile Deformable Surface Representation",
//      ACM Trans. Graph. 25/1, pp. 151-175, 2006.
//the sparse field level set method, which is used for time integration is described in
//      R.T.~Whitaker, A Level-Set Approach to 3D Reconstruction from Range Data
//      J. Comp. Vision 29/3, pp. 203-231, 1998.
//      in this file functions are defined which rebuild the level set function after the active grid points are time integrated
//      the functions for the time integration itself are defined in "timeintegration.hpp"

namespace lvlset {

    double allocation_factor=1.2;

    // types and constants

    enum sign_type {POS_SIGN=0, NEG_SIGN=1};        //type used for the sign

    typedef unsigned int direction_type;            //direction type
    const direction_type X_DIRECTION=0;             //constants for the directions
    const direction_type Y_DIRECTION=1;
    const direction_type Z_DIRECTION=2;


    class DefaultLevelSetTraitsType {               //the default level set traits type
    public:
        typedef unsigned int size_type;             //type used for indexing all arrays
        typedef double value_type;                  //data type used for the level set values
    };

    // the level set class
    template <class GridTraitsType, class LevelSetTraitsType=DefaultLevelSetTraitsType>
    class levelset {
        // the class GridTraitsType has to have the properties as described in "grid.hpp"

        //the LevelSetTraitsType defines the used data types within the level
        //set data structure, and has to have the certain type and static constant definitions
        //as shown in the "DefaultLevelSetTraitsType" above

        //type and constant redefinitions

        static const int D=GridTraitsType::dimensions;

        void operator=(const levelset& l) const;

    public:
         //type redefinitions
        typedef typename LevelSetTraitsType::value_type value_type;
        typedef typename LevelSetTraitsType::size_type size_type;
        typedef typename GridTraitsType::index_type index_type;

        typedef lvlset::grid_type<GridTraitsType> grid_type2;

        typedef vec<index_type,D> point_type;
        typedef std::vector<point_type> points_type;

    private:

        class allocation_type {
        public:
            vec<size_type,D> num_values;
            vec<size_type,D> num_runs;

            template<class X>
            allocation_type& operator*=(const X& x) {
                for (int i=0;i<D;++i) num_values[i]=static_cast<size_type>(num_values[i]*x);
                for (int i=0;i<D;++i) num_runs[i]=static_cast<size_type>(num_runs[i]*x);
                return *this;
            }

            allocation_type& operator+=(const allocation_type& x) {
                num_values+=x.num_values;
                num_runs+=x.num_runs;
                return *this;
            }

            template<class X>
            allocation_type&  operator/=(const X& x) {
                for (int i=0;i<D;++i) num_values[i]=static_cast<size_type>(num_values[i]/x)+1;
                for (int i=0;i<D;++i) num_runs[i]=static_cast<size_type>(num_runs[i]/x)+1;
                return *this;
            }

            template<class X>
            allocation_type operator*(const X& x) {
                allocation_type tmp(*this);
                tmp*=x;
                return tmp;
            }

            template<class X>
            allocation_type operator/(const X& x) {
                allocation_type tmp(*this);
                tmp/=x;
                return tmp;
            }

            template<class X>
            allocation_type operator+(const X& x) {
                allocation_type tmp(*this);
                tmp+=x;
                return tmp;
            }


            allocation_type():num_values(size_type(0)), num_runs(size_type(0)) {}
        };

        //NOTE: This vector is not used anywhere.
        typedef std::vector<allocation_type> allocations_type;


        // [Josef] The level set class, containing the Grid class.
        // the level set class
        class sub_levelset_type {
        public:

             //type redefinitions
            typedef typename LevelSetTraitsType::value_type value_type;
            typedef typename LevelSetTraitsType::size_type size_type;
            typedef typename GridTraitsType::index_type index_type;

        public:

            //the following vectors are used to store the run-lenght-encoded data structure
            //for more details on that data structure see B. Houston, M.B. Nielsen, C. Batty, O. Nilsson, K. Museth,
            //"Hierarchical RLE-Level Set - A Compact and Versatile Deformable Surface Representation", ACM Trans. Graph. 25/1, pp. 151-175, 2006.
            std::vector<size_type> start_indices[D];
            std::vector<size_type> runtypes[D];
            std::vector<index_type> runbreaks[D];
            std::vector<value_type> distances;          //this vector keeps the level set values of all defined grid points
                                                        //its size is therefore equal to the number of defined grid points

            std::vector<size_type> active;  //an additional vector, which has the same size as the distances-vector
                                            //in case the defined grid point is active in terms of the sparse field level set method
                                            //(the level set value is in the range [-0.5,0.5])
                                            //this vector keeps the active point id
                                            //the active point ids are given to the active grid points in lexicographical order starting with index 0
                                            //therefore the largest active point ID equals the number of active grid points-1
                                            //if the grid point is not active, its value is set to the constant INACTIVE


            size_type num_active_points;    //num_active_points stores the number of active points



            const grid_type2& Grid; //Grid stores the information about the grid, on which the level set function is defined

        public:

            allocation_type get_allocation() const {
                 // allocation_type allocates the requried sizes to num_values and num_runs
                 // num_values[0] is to contain level set values, num_values[i] contains the start indices at the i-th dimension
                 // num_runs[i] is to contain the run types at the i-th dimension
                 allocation_type a;
                 a.num_values[0]=distances.size();  // To contain level set values
                 a.num_runs[0]=runtypes[0].size();
                 for (int i=1;i<D;++i) {
                     a.num_values[i]=start_indices[i-1].size();
                     a.num_runs[i]=runtypes[i].size();
                 }

                 // Check that the vector sizes make sense for the H-RLE structure
                 assert(runbreaks[D-1].size()==a.num_runs[D-1]-1);
                 for (int i=1;i<D;++i) {
                     assert(runbreaks[i-1].size()==a.num_runs[i-1]-a.num_values[i]);
                 }
                 return a;
             }

            sub_levelset_type(const grid_type2& g, const allocation_type& a) :num_active_points(0),  Grid(g) {
                 // Constructor for the sub_levelset_type class
                 // Takes a pointer to a grid_type and an allocation_type to reserve the memory required for the
                 // start indeces, run types and run breaks arays, as well as the LS values

                 distances.reserve(a.num_values[0]);
                 runtypes[0].reserve(a.num_runs[0]);

                 for (int i=1;i<D;++i) {
                     start_indices[i-1].reserve(a.num_values[i]);
                     runbreaks[i-1].reserve(a.num_runs[i-1]-a.num_values[i]);
                     runtypes[i].reserve(a.num_runs[i]);
                 }

                 start_indices[D-1].push_back(0);
                 runbreaks[D-1].reserve(a.num_runs[D-1]-1);
             };

            const sub_levelset_type & operator=(const sub_levelset_type& s) {
                 for (int i=0;i<D;++i) start_indices[i]=s.start_indices[i];
                 for (int i=0;i<D;++i) runtypes[i]=s.runtypes[i];
                 for (int i=0;i<D;++i) runbreaks[i]=s.runbreaks[i];
                 distances=s.distances;
                 active=s.active;
                 num_active_points=s.num_active_points;
                 assert((&Grid)==(&s.Grid));
                 return *this;
             }

            index_type GetRunStartCoord(int dim, size_type start_indices_pos, size_type run_type_pos) const {
                //returns the start index of the run given by start_indices_pos and run_type_pos
                //see the definition of the HRLE-data structure for more details
                if (run_type_pos==start_indices[dim][start_indices_pos]) {
                    return Grid.min_point_index(dim);
                } else {
                    return runbreaks[dim][run_type_pos-start_indices_pos-1];
                }
            }

            index_type GetRunEndCoord(int dim, size_type start_indices_pos, size_type run_type_pos) const {
                //returns the end index of the run given by start_indices_pos and run_type_pos
                //NOTE: the end index is not included by the run
                //see the definition of the HRLE-data structure for more details
                if (run_type_pos+1<GetStartIndex(dim,start_indices_pos+1)) {
                    return runbreaks[dim][run_type_pos-start_indices_pos]-1;
                } else {
                    return Grid.max_point_index(dim);
                }
            }

            size_type GetStartIndex(int dim, size_type start_indices_pos) const {
                //returns the starting index of the runtypes array
                //if start_indices_pos is equal to the size of start_indices array
                //the size of the run_types array is returned
                if (start_indices_pos==start_indices[dim].size()) {
                    return runtypes[dim].size();
                } else {
                    return start_indices[dim][start_indices_pos];
                }
            }

            index_type getRunBreak(int dim, int runbreak = std::numeric_limits<int>::max()) const{
                if(runbreak==std::numeric_limits<int>::max()) return runbreaks[dim].back();
                return runbreaks[dim][runbreak];
            }

            index_type getMaxRunBreak(int dim) const{
                return *std::max_element(runbreaks[dim].begin(), runbreaks[dim].end());
            }

            index_type getMinRunBreak(int dim) const{
                return *std::min_element(runbreaks[dim].begin(), runbreaks[dim].end());
            }

            unsigned long int used_memory() const {
                //this function gives an estimation of the used memory of the level
                //set function. however, the allocated memory is much higher because the
                //STL-vector allocates more memory than needed, to accelerate push_back-operations
                //the allocated memory can be obtained by the "allocated_memory"-function defined below
                //NOTE: since this function is based on  the sizeof-operator it cannot take the memory needed for the
                //      GridTraitsType (see grid.hpp) into account if dynamic data structures are used for example
                //      Furthermore due to data structure alignment the result can be somehow inaccurate

                unsigned long int x=sizeof(sub_levelset_type);

                for (int i=0;i<D;i++) {
                    x+=sizeof(size_type)*start_indices[i].size();
                    x+=sizeof(size_type)*runtypes[i].size();
                    x+=sizeof(index_type)*runbreaks[i].size();
                }
                x+=sizeof(value_type)*distances.size();
                x+=sizeof(size_type)*active.size();

                return x;

            }

            unsigned long int allocated_memory() const {
                //this function gives an estimation of the allocated memory for the level
                //set function.
                //NOTE: since this function is based on  the sizeof-operator it cannot take the memory needed for the
                //      GridTraitsType (see grid.hpp) into account if dynamic data structures are used for example
                //      Furthermore due to data structure alignment the result can be somehow inaccurate

                unsigned long int x=sizeof(sub_levelset_type);

                for (int i=0;i<D;i++) {
                    x+=sizeof(size_type)*start_indices[i].capacity();
                    x+=sizeof(size_type)*runtypes[i].capacity();
                    x+=sizeof(index_type)*runbreaks[i].capacity();
                }
                x+=sizeof(value_type)*distances.capacity();
                x+=sizeof(size_type)*active.capacity();

                return x;
            }

            size_type num_pts() const {
                //this function returns the number of defined grid points
                return distances.size();
            };

            size_type num_active_pts() const {
                //this function returns the number of active grid points
                 return num_active_points;
            }

            size_type number_of_runs(int level) const {
                if (level>=0) {
                    return runtypes[level].size();
                } else {
                    return distances.size();
                }
            }


            template <class V>
            void push_back_undefined(V start_point, const V& end_point, size_type rt) {


                if (start_point>end_point) return;  // in this case, do not add the point

                for (int dim=0;dim<D-1;++dim) {

                    if (start_point[dim]!=Grid.min_point_index(dim)) {

                        if (start_point[dim]<=Grid.max_point_index(dim)){
                          push_back_undefined(start_point, rt);
                        }
                        start_point[dim]=Grid.min_point_index(dim);
                        ++start_point[dim+1];
                    }

                    if (start_point>end_point) return;
                }

                if (start_point[D-1]<=Grid.max_point_index(D-1)) push_back_undefined(start_point, rt);
            }

            template <class V>
            void push_back_undefined(const V& point, size_type rt) {

              int level;
              for (level=0;level<D;++level) {
                if (point[level]!=Grid.min_point_index(level)) break;
              }

              size_type old_sign=0;
              int dim;

              for(dim=D-1;dim>level;--dim) {

                if (runtypes[dim].size()==start_indices[dim].back()) {    //if there is no run
                  if(point[dim]!=Grid.min_point_index(dim)) {
                    runtypes[dim].push_back(old_sign);
                    runbreaks[dim].push_back(point[dim]);
                  }
                  runtypes[dim].push_back(start_indices[dim-1].size());
                  start_indices[dim-1].push_back(runtypes[dim-1].size());
                } else if (!levelset::is_defined(runtypes[dim].back())) {          //if there is an defined run
                  old_sign=runtypes[dim].back();
                  if (old_sign==rt) return;
                  if (runtypes[dim].size()==start_indices[dim].back()+1) {   //if there is a single run
                    if (point[dim]==Grid.min_point_index(dim)) {
                      runtypes[dim].back()=start_indices[dim-1].size();
                    } else {
                      runbreaks[dim].push_back(point[dim]);
                      runtypes[dim].push_back(start_indices[dim-1].size());
                    }
                  } else {                                                    //if there are more than one runs
                      if (point[dim]==runbreaks[dim].back()) {
                        runtypes[dim].pop_back();
                        if (!levelset::is_defined(runtypes[dim].back())){
                          runtypes[dim].push_back(start_indices[dim-1].size());
                        } else runbreaks[dim].pop_back();
                    } else {
                      runbreaks[dim].push_back(point[dim]);
                      runtypes[dim].push_back(start_indices[dim-1].size());
                    }
                  }
                  start_indices[dim-1].push_back(runtypes[dim-1].size());
                }
              }

              if (runtypes[dim].size()==start_indices[dim].back()) {    //if there is no run
                if(point[dim]!=Grid.min_point_index(dim)) {
                  runtypes[dim].push_back(old_sign);
                  runbreaks[dim].push_back(point[dim]);
                }
                runtypes[dim].push_back(rt);
              } else if (!levelset::is_defined(runtypes[dim].back())) {          //if there is an defined run
                old_sign=runtypes[dim].back();
                if (old_sign==rt) return;
                if (runtypes[dim].size()==start_indices[dim].back()+1) {   //if there is a single run
                  if (point[dim]==Grid.min_point_index(dim)) {
                    runtypes[dim].back()=rt;
                  } else {
                    runbreaks[dim].push_back(point[dim]);
                    runtypes[dim].push_back(rt);
                  }
                } else {                                                    //if there are more than one runs
                  if (point[dim]==runbreaks[dim].back()) {
                    runtypes[dim].back()=rt;
                  } else {
                    runbreaks[dim].push_back(point[dim]);
                    runtypes[dim].push_back(rt);
                  }
                }
              } else {
                runbreaks[dim].push_back(point[dim]);
                runtypes[dim].push_back(rt);
              }
            }

            template <class V>
            void push_back(const V& point, value_type distance) {

                int level;
                for (level=0;level<D;++level) {
                    if (point[level]!=Grid.min_point_index(level)) break;
                }

                size_type old_sign=0;

                for(int dim=D-1;dim>0;--dim) {
                    if (runtypes[dim].size()==start_indices[dim].back()) {    //if there is no run
                       if(point[dim]!=Grid.min_point_index(dim)) {
                           runtypes[dim].push_back(old_sign);
                           runbreaks[dim].push_back(point[dim]);
                       }
                       runtypes[dim].push_back(start_indices[dim-1].size());
                       start_indices[dim-1].push_back(runtypes[dim-1].size());
                    } else if (!is_defined(runtypes[dim].back())) {          //if there is an defined run
                        old_sign=runtypes[dim].back();
                        if (runtypes[dim].size()==start_indices[dim].back()+1) {   //if there is a single run
                            if (point[dim]==Grid.min_point_index(dim)) {
                                runtypes[dim].back()=start_indices[dim-1].size();
                            } else {
                                runbreaks[dim].push_back(point[dim]);
                                runtypes[dim].push_back(start_indices[dim-1].size());
                            }
                        } else {                                                    //if there are more than one runs
                            if (point[dim]==runbreaks[dim].back()) {
                                runtypes[dim].pop_back();
                                if (!is_defined(runtypes[dim].back())) {
                                    runtypes[dim].push_back(start_indices[dim-1].size());
                                } else {
                                    runbreaks[dim].pop_back();
                                }
                            } else {
                                runbreaks[dim].push_back(point[dim]);
                                runtypes[dim].push_back(start_indices[dim-1].size());
                            }
                        }
                        start_indices[dim-1].push_back(runtypes[dim-1].size());
                    } else {
                        if (dim<=level) start_indices[dim-1].push_back(runtypes[dim-1].size());
                    }
                }

                if (runtypes[0].size()==start_indices[0].back()) {    //if there is no run
                   if(point[0]!=Grid.min_point_index(0)) {
                       runtypes[0].push_back(old_sign);
                       runbreaks[0].push_back(point[0]);
                   }
                   runtypes[0].push_back(distances.size());
                } else if (!is_defined(runtypes[0].back())) {          //if there is an defined run
                    old_sign=runtypes[0].back();
                    if (runtypes[0].size()==start_indices[0].back()+1) {   //if there is a single run
                        if (point[0]==Grid.min_point_index(0)) {
                            runtypes[0].back()=distances.size();
                        } else {
                            runbreaks[0].push_back(point[0]);
                            runtypes[0].push_back(distances.size());
                        }
                    } else {                                                    //if there are more than one runs
                        if (point[0]==runbreaks[0].back()) {
                            runtypes[0].pop_back();
                            if (!is_defined(runtypes[0].back())) {
                                runtypes[0].push_back(distances.size());
                            } else {
                                runbreaks[0].pop_back();
                            }
                        } else {
                            runbreaks[0].push_back(point[0]);
                            runtypes[0].push_back(distances.size());
                        }
                    }
                }

                distances.push_back(distance);

                if (math::abs(distance)<=value_type(0.5)) {
                    active.push_back(num_active_points);
                    ++num_active_points;

                } else {
                    active.push_back(INACTIVE);
                }

            }

            void invert() {
                //this function inverts a level set function
                //by switching the sign of all level set values
                //and by interchanging for all undefined runs
                //POS_PT and NEG_PT

                for (size_type a=0;a<distances.size();a++) {
                    distances[a]=-distances[a];
                }

                for (int b=0;b<D;b++) {
                    for (size_type c=0;c<runtypes[b].size();c++) {
                        if (runtypes[b][c]==POS_PT) {
                            runtypes[b][c]=NEG_PT;
                        } else if (runtypes[b][c]==NEG_PT) {
                            runtypes[b][c]=POS_PT;
                        }
                    }
                }
            }

            void print(std::ostream& out = std::cout) const {
                std::ostringstream oss;

                out << std::endl;
                out << "levelset data structure" << std::endl << std::endl;

                for (int dim=D-1;dim>=0;--dim) {
                    out <<  dim <<  " start_indices: " << start_indices[dim].size();
                    int c=0;
                    for (typename std::vector<size_type>::const_iterator it=start_indices[dim].begin();it!=start_indices[dim].end();++it) {
                        if (c%10==0) out << std::endl;
                        ++c;
                        out << std::setw(8) << *it;
                    }
                    out << std::endl;

                    out << dim << " run_types: " << runtypes[dim].size();
                    c=0;
                    for (typename std::vector<size_type>::const_iterator it=runtypes[dim].begin();it!=runtypes[dim].end();++it) {
                        if (c%10==0) out << std::endl;
                        ++c;
                        if ((*it)==POS_PT) {
                            out << std::setw(8) << "+oo";
                        } else if ((*it)==NEG_PT) {
                            out << std::setw(8) << "-oo";
                        } else if ((*it)==UNDEF_PT) {
                            out << std::setw(8) << "UND";
                        } else if (!is_defined(*it)) {
                            oss.str("");
                            oss << "S" << ((*it) - SEGMENT_PT);
                            out << std::setw(8) << oss.str();
                        } else {
                            out << std::setw(8) << (*it);
                        }
                    }
                    out << std::endl;

                    out << dim << " run_breaks: " << runbreaks[dim].size();
                    c=0;
                    for (typename std::vector<index_type>::const_iterator it=runbreaks[dim].begin();it!=runbreaks[dim].end();++it) {
                        if (c%10==0) out << std::endl;
                        ++c;
                        out << std::setw(8) << *it;
                    }
                    out << std::endl;
                }

                out << "distances: " << distances.size();
                int c=0;
                for (typename std::vector<value_type>::const_iterator it=distances.begin();it!=distances.end();++it) {
                    if (c%10==0) out << std::endl;
                    ++c;
                    out << std::setw(8) << std::fixed << std::setprecision(3) << *it;
                }
                out << std::endl;
                out << std::endl;
            }

        };

        class sub_levelsets_type {
        // contains the vector of sub_levelset_type items that are used for parallelization

            typedef sub_levelset_type* sub_levelset_ptr_type;  // pointer to memory containing sub_levelset_type object

            typedef typename std::vector< sub_levelset_ptr_type> sub_levelsets_intern_type;
            sub_levelsets_intern_type subs;

            const sub_levelsets_type& operator=(const sub_levelsets_type& s) {}

        public:

            typedef typename sub_levelsets_intern_type::size_type size_type;

            const sub_levelset_type & operator[](size_type i) const {
                return *(subs[i]);
            }

            sub_levelset_type & operator[](size_type i) {
                return *(subs[i]);
            }

            const sub_levelset_type & back() const {
                return *(subs.back());
            }

            sub_levelset_type & back() {
                return *(subs.back());
            }

            void swap(sub_levelsets_type& s) {
                subs.swap(s.subs);
            }

            void initialize(size_type i, const grid_type2& g,  allocation_type a=allocation_type()) {
                // initialize the parallelization by implementing i number of sub_levelset_type objects
                // through a sub_levelset_ptr_type object, that each correspond to a parallel core


                for (typename sub_levelsets_intern_type::iterator it=subs.begin();it!=subs.end();++it) if (*it!=0) delete *it;

                subs.clear();
                subs.resize(i);//,sub_levelset_ptr_type(0));

                a*=allocation_factor;
                a/=i;

                #pragma omp parallel for schedule(static,1)  // parallelization - Iterations divided into chunks of size 1. Each chunk is assigned to a thread
                for (int k=0;k<static_cast<int>(subs.size());++k) {
                    subs[k]= sub_levelset_ptr_type(new sub_levelset_type(g,a));
                }
            }

            size_type size() const {
                // return the size of the vector for parallelization
                return subs.size();
            }

            sub_levelsets_type(const sub_levelsets_type& s) {
                // builds a new sub_levelsets_type that copies the one passed to the function as s
                subs.resize(s.size());

                #pragma omp parallel for schedule(static,1)// parallelization - Iterations divided into chunks of size 1. Each chunk is assigned to a thread
                for (int k=0;k<static_cast<int>(s.size());++k) {
                    subs[k]= sub_levelset_ptr_type(new sub_levelset_type(s[k]));
                }
            }

            sub_levelsets_type() {}

            ~sub_levelsets_type() {
                // sub_levelsets_type destructor used to delete all pointers for clean-up
                for (typename sub_levelsets_intern_type::iterator it=subs.begin();it!=subs.end();++it) if (*it!=0) delete *it;
            }
        };

        /*class sub_levelsets_type {

            typedef typename std::vector<sub_levelset_type> sub_levelsets_intern_type;
            sub_levelsets_intern_type subs;

        public:

            typedef typename sub_levelsets_intern_type::size_type size_type;

            const sub_levelset_type & operator[](size_type i) const {
                return (subs[i]);
            }

            sub_levelset_type & operator[](size_type i) {
                return (subs[i]);
            }

            const sub_levelset_type & back() const {
                return (subs.back());
            }

            sub_levelset_type & back() {
                return (subs.back());
            }

            void swap(sub_levelsets_type& s) {
                subs.swap(s.subs);
            }

            void initialize(size_type i, const grid_type2& g, const allocation_type& a=allocation_type()) {

                subs.clear();
                subs.resize(i,sub_levelset_type(g));
            }

            size_type size() const {
                return subs.size();
            }

        };*/





    public:



        //typedef std::vector<sub_levelset_type> sub_levelsets_type;

        //typedef typename GridTraitsType::coord_type coord_type;

        //###############################################
        // Declaration of contants
        //###############################################

        static const size_type INACTIVE;    //this constant is used for grid points which are not active

        static const value_type POS_VALUE;  //this constant is used for positive undefined grid points as level set value
        static const value_type NEG_VALUE;  //this constant is used for negative undefined grid points as level set value

        static const size_type POS_PT;      //this constant is used as runtype for positive runs
        static const size_type NEG_PT;      //this constant is used as runtype for negative runs
        static const size_type UNDEF_PT;    //this constant is used as runtype for uninitialized runs
        static const size_type SEGMENT_PT;    //this constant is used as runtype for uninitialized runs

        //Public dimension, D is private
        static const int dimensions=GridTraitsType::dimensions;

        //########################################################################
        //filter definitions
        //-------------------
        //these filters which are basically functors returning an boolean
        //are used for certain iterators, to specify at which grid points an iterator should stop

        class filter_all_defined {      //stop at all defined grid points
        public:
            template <class C> bool operator()(const C& c) const {return c.is_defined();}
        };

        class filter_all {              //stop at all possible stop points (beginning of each run)
        public:                         //in case of undefined runs of the data structure, which can consist of more than 1 grid point
                                        //the iterator will also stop at the beginning of such runs
            template <class C> bool operator()(const C& c) const {return true;}
        };

        class filter_none {             //never stop
        public:
            template <class C> bool operator()(const C& c) const {return false;}
        };

        class filter_active {           //stop at active grid points
        public:
            template <class C> bool operator()(const C& c) const {return c.is_active();}
        };

        class filter_value {            //stop at (defined) grid points which have an absolute level set value
            value_type limit;           //less than a given limit
        public:
            filter_value(value_type x) : limit(x) {}
            template <class C> bool operator()(const C& c) const {return (math::abs(c.value())<=limit);}
        };

        //######################################################################


        /*class local_pt_id_type {
        public:
            int sub;
            size_type id;
            local_pt_id_type(int s, size_type i):sub(s), id(i) {}
        };

        class global_pt_id_type {
        public:
            size_type id;
            global_pt_id_type(size_type i):id(i) {}
        };*/



    private:

        //the following vectors are used to store the run-lenght-encoded data structure
        //for more details on that data structure see B. Houston, M.B. Nielsen, C. Batty, O. Nilsson, K. Museth,
        //"Hierarchical RLE-Level Set - A Compact and Versatile Deformable Surface Representation", ACM Trans. Graph. 25/1, pp. 151-175, 2006.
        //std::vector<size_type> start_indices[D];
        //std::vector<size_type> runtypes[D];
        //std::vector<index_type> runbreaks[D];
        //std::vector<value_type> distances;          //this vector keeps the level set values of all defined grid points
                                                    //its size is therefore equal to the number of defined grid points

        //std::vector<size_type> active;  //an additional vector, which has the same size as the distances-vector
                                        //in case the defined grid point is active in terms of the sparse field level set method
                                        //(the level set value is in the range [-0.5,0.5])
                                        //this vector keeps the active point id
                                        //the active point ids are given to the active grid points in lexicographical order starting with index 0
                                        //therefore the largest active point ID eqauls the number of active grid points-1
                                        //if the grid point is not active, its value is set to the constant INACTIVE

        points_type segmentation;

        sub_levelsets_type sub_levelsets;


        int num_layers; //num_layers keeps the number of layers of which the level set function
                        //currently consists of at least
                        //the minimum valid value is 1, that means that at least all grid points are explicitly stored (defined grid points)
                        //which have a value in the range [-0.5,0.5]
                        //if num_layers=2 then the level set data structures stores at least all grid points, which are
                        //in the range [-1,1]
                        //if num_layers=n at least all grid points are explicitly stored which are in the range [-n/2,n/2]

        const grid_type<GridTraitsType>& Grid; //Grid stores the information about the grid, on which the level set function is defined

        std::vector<size_type> active_pt_id_offsets;
        std::vector<size_type> pt_id_offsets;

        unsigned levelset_id;

        static unsigned top_levelset_id;




        //size_type num_active_points;    //num_active_points stores the number of active points



        /*template <class V> V GetNextCoordinate(V coords, int& level) const {
            //this functions returns the next indices of the grid in lexicographical order
            //first coords[0] is increased, if it is out of range of the grid, it is set to the minimum index
            //and coords[1] is incresed and so on
            //"level" returns the maximum direction, for which the index is altered
            //level==dimensions is returned if the end of the grid is reached
            level=0;
            while(true) {
                ++coords[level];
                if (coords[level]<=(Grid.max_point_index(level))) break;
                coords[level]=Grid.min_point_index(level);
                ++level;
                if (level==D) {
                    coords[D-1]=Grid.max_point_index(D-1)+1;
                    break;
                }
            }
            return coords;
        }*/


        /*index_type GetRunStartCoord(int dim, size_type start_indices_pos, size_type run_type_pos) const {
            //returns the start index of the run given by start_indices_pos and run_type_pos
            //see the definition of the HRLE-data structure for more details
            if (run_type_pos==start_indices[dim][start_indices_pos]) {
                return Grid.min_point_index(dim);
            } else {
                return runbreaks[dim][run_type_pos-start_indices_pos-1];
            }
        }

        index_type GetRunEndCoord(int dim, size_type start_indices_pos, size_type run_type_pos) const {
            //returns the end index of the run given by start_indices_pos and run_type_pos
            //NOTE: the end index is not included by the run
            //see the definition of the HRLE-data structure for more details
            if (run_type_pos+1<GetStartIndex(dim,start_indices_pos+1)) {
                return runbreaks[dim][run_type_pos-start_indices_pos]-1;
            } else {
                return Grid.max_point_index(dim);
            }
        }


        size_type GetStartIndex(int dim, size_type start_indices_pos) const {
            //returns the starting index of the runtypes array
            //if start_indices_pos is equal to the size of start_indices array
            //the size of the run_types array is returned
            if (start_indices_pos==start_indices[dim].size()) {
                return runtypes[dim].size();
            } else {
                return start_indices[dim][start_indices_pos];
            }
        }*/

        /*size_type FindInferiorRunTypePos(int dim, size_type start_indices_pos, index_type relative_coord) const {
            //this functions finds the right run by looking up the given index "relative_coord" in the
            //runbreaks-array, and returns the index to that run
            const size_type start=GetStartIndex(dim, start_indices_pos);
            const size_type end=GetStartIndex(dim, start_indices_pos+1);

            const size_type start_break=start-start_indices_pos;
            const size_type end_break=end-start_indices_pos-1;

            return start+(std::upper_bound(runbreaks[dim].begin()+start_break,runbreaks[dim].begin()+end_break,relative_coord)-(runbreaks[dim].begin()+start_break));
        }*/



    public:

        void print_segment_sizes(){
            for (typename sub_levelsets_type::size_type i=0;i!=sub_levelsets.size();++i) {
                std::cout << sub_levelsets[i].distances.size() << std::endl;
            }
        }

        void print(std::ostream& out = std::cout) const {
            for (typename sub_levelsets_type::size_type i=0;i!=sub_levelsets.size();++i) {
                sub_levelsets[i].print(out);
                out << std::endl;
            }
        }

        void print_without_segmentation(std::ostream& out = std::cout){
          serialize();
          print(out);
          segment();
        }

        //returns a non-const reference to the start indices of a segment
        std::vector<size_type>& start_indices(unsigned int dim, unsigned int segment = 0){
          return sub_levelsets[segment].start_indices[dim];
        }
        //returns a non-const reference to the runtypes of a segment
        std::vector<size_type>& runtypes(unsigned int dim, unsigned int segment = 0){
          return sub_levelsets[segment].runtypes[dim];
        }
        //returns a non-const reference to the runbreaks of a segment
        std::vector<index_type>& runbreaks(unsigned int dim, unsigned int segment = 0){
          return sub_levelsets[segment].runbreaks[dim];
        }
        //returns a non-const reference to the distances of a segment
        std::vector<value_type>& distances(unsigned int segment = 0){
          return sub_levelsets[segment].distances;
        }

        //writes the levelset to a file
        levelset& export_levelset(const std::string& path, int bits_per_distance = 8){
          export_levelset_to_file(*this, path, bits_per_distance);
          return *this;
        }

        //reads the levelset from a levelset file; NOTE: grid has to be read first!!
        levelset& import_levelset(const std::string& path){
          import_levelset_from_file(*this, path);
          return *this;
        }


        void swap(levelset<GridTraitsType, LevelSetTraitsType>& l) {
            //this function swaps the level set function with another level set function "l"
            segmentation.swap(l.segmentation);
            sub_levelsets.swap(l.sub_levelsets);
            pt_id_offsets.swap(l.pt_id_offsets);
            active_pt_id_offsets.swap(l.active_pt_id_offsets);

            //sub_levelsets.update_levelset_ref(*this);
            //l.sub_levelsets.update_levelset_ref(l);

            /*for (typename sub_levelsets_type::iterator it=sub_levelsets.begin();it!=sub_levelsets.end();++it) {
                it->l=this;
            }

            for (typename sub_levelsets_type::iterator it=l.sub_levelsets.begin();it!=l.sub_levelsets.end();++it) {
                it->l=&l;
            }*/


            std::swap(num_layers, l.num_layers);

            assert(&Grid==&l.Grid);
            //std::swap(Grid, l.Grid);
        }

        unsigned long int used_memory() const {
            //this function gives an estimation of the used memory of the level
            //set function. however, the allocated memory is much higher because the
            //STL-vector allocates more memory than needed, to accelerate push_back-operations
            //the allocated memory can be obtained by the "allocated_memory"-function defined below
            //NOTE: since this function is based on  the sizeof-operator it cannot take the memory needed for the
            //      GridTraitsType (see grid.hpp) into account if dynamic data structures are used for example
            //      Furthermore due to data structure alignment the result can be somehow inaccurate

            //TODO

            unsigned long int x=sizeof(levelset<GridTraitsType, LevelSetTraitsType>);

            for (unsigned int i=0;i<sub_levelsets.size();i++) {
                x+=sub_levelsets[i].used_memory();
            }

            return x;

        }

        unsigned long int allocated_memory() const {
            //this function gives an estimation of the allocated memory for the level
            //set function.
            //NOTE: since this function is based on  the sizeof-operator it cannot take the memory needed for the
            //      GridTraitsType (see grid.hpp) into account if dynamic data structures are used for example
            //      Furthermore due to data structure alignment the result can be somehow inaccurate

            unsigned long int x=sizeof(levelset<GridTraitsType, LevelSetTraitsType>);

            for (unsigned int i=0;i<sub_levelsets.size();i++) {
                x+=sub_levelsets[i].allocated_memory();
            }
            return x;

        }



        size_type number_of_runs(int level, int sub) const {
            //this function returns the number of runs for the given sub_levelset and given level

            return sub_levelsets[sub].number_of_runs(level);
        }

        typename sub_levelsets_type::size_type number_of_segments() const {
            return sub_levelsets.size();
        }

        index_type get_runbreak(int dim, int runbreak = std::numeric_limits<int>::max()) const{
            return sub_levelsets[0].getRunBreak(dim, runbreak);
        }

        index_type get_max_runbreak(int dim) const{
            index_type maxBreak = sub_levelsets[0].getMaxRunBreak(dim);
            for(unsigned int i=1; i<number_of_segments(); ++i) maxBreak = std::max(maxBreak, sub_levelsets[i].getMaxRunBreak(dim));
            return maxBreak;
        }

        index_type get_min_runbreak(int dim) const{
            index_type minBreak = sub_levelsets[0].getMinRunBreak(dim);
            for(unsigned int i=1; i<number_of_segments(); ++i) minBreak = std::min(minBreak, sub_levelsets[i].getMinRunBreak(dim));
            return minBreak;
        }

        const grid_type<GridTraitsType>& grid() const {
            //this function returns a constant reference to the grid,
            //on which the level set function is defined
            return Grid;
        }

        const unsigned get_levelset_id() const{
          return levelset_id;
        }

        void set_levelset_id(){
          levelset_id = top_levelset_id;
          ++top_levelset_id;
        }


        //###################################################
        // functions to get the point-ID for a given indices

        template <class CoordType>
        size_type pt_id2(const CoordType& coords) const {
            //this function returns the ID of a defined grid point given by its local indices "coords"
            //or in case of undefined grid points it returns POS_PT or NEG_PT
            //here, local indices means that the indices must be contained by the grid
            //see also the function pt_id below

            //size_type sub=std::upper_bound(segmentation.begin(), segmentation.end(), coords)-segmentation.begin();

            //size_type distances_pos(0);     //start_indices_pos

          int sub=std::upper_bound(segmentation.begin(), segmentation.end(), coords)-segmentation.begin();

            //shfdhsfhdskjhgf assert(sub>=0);
            //shfdhsfhdskjhgf assert(sub<sub_levelsets.size());

            size_type id=0;

            for (int level=D-1;level>=0;--level) {

                size_type run_type_pos=sub_levelsets[sub].start_indices[level][id];

                typename std::vector<index_type>::const_iterator start_breaks = sub_levelsets[sub].runbreaks[level].begin()+(run_type_pos-id);
                typename std::vector<index_type>::const_iterator end_breaks= sub_levelsets[sub].runbreaks[level].begin()+(sub_levelsets[sub].GetStartIndex(level, id+1)-(id+1));
                typename std::vector<index_type>::const_iterator pos_breaks=std::upper_bound(
                                start_breaks,
                                end_breaks,
                                coords[level]);

                run_type_pos+=(pos_breaks-start_breaks);

                id=sub_levelsets[sub].runtypes[level][run_type_pos];

                if (!is_defined(id)) {
                    //shfdhsfhdskjhgf assert(id==POS_PT || id==NEG_PT);
                    return id;
                }

                if (pos_breaks==start_breaks) {
                    id+=(coords[level]-Grid.min_point_index(level));
                } else {
                    id+=(coords[level]- (*(pos_breaks-1)));
                }

            }

            //shfdhsfhdskjhgf assert(is_defined(id));

            return id;
        }

        template <class CoordType>
        size_type pt_id(const CoordType& coords) const {
            //this function returns the ID of a grid point given by its global indices "coords"
            //or in case of undefined grid points it returns POS_PT or NEG_PT
            //here, global indices means that also indices can be passed to this function, which are not included
            //by the grid, or in other words the boundary conditions are taken into account to get the right grid point
            return pt_id2(Grid.global_indices_2_local_indices(coords));
        }

        //#############################################
        // functions to access level set values
        //##############################################

        const value_type& value(size_type pt_id) const {       //this function returns the level set value
            if (pt_id==POS_PT) return POS_VALUE;               //of a grid point given by its pt_id
            if (pt_id==NEG_PT) return NEG_VALUE;
            return value2(pt_id);
        }

        const value_type& value(int sub, size_type pt_id) const {       //this function returns the level set value
            if (pt_id==POS_PT) return POS_VALUE;                        //of a grid point given by its pt_id
            if (pt_id==NEG_PT) return NEG_VALUE;
            return value2(sub,pt_id);
        }

        const value_type& value2(size_type pt_id) const {
            typename std::vector<size_type>::const_iterator it= std::upper_bound(pt_id_offsets.begin(),pt_id_offsets.end(),pt_id)-1;
            return sub_levelsets[it-pt_id_offsets.begin()].distances[pt_id-(*it)];
        }

        const value_type& value2(int sub,size_type pt_id) const {
            return sub_levelsets[sub].distances[pt_id];
        }

        void set_value(size_type pt_id, value_type v) {     //this function sets the level set value of a grid point
            typename std::vector<size_type>::const_iterator it= std::upper_bound(pt_id_offsets.begin(),pt_id_offsets.end(),pt_id)-1;
            sub_levelsets[it-pt_id_offsets.begin()].distances[pt_id-(*it)]=v;
        }

        void set_value(int sub, size_type pt_id, value_type v) {     //this function sets the level set value of a grid point
            sub_levelsets[sub].distances[pt_id]=v;
        }

        template <class CoordType>
        const value_type& value(const CoordType& c) const { //this function returns the level set value
            return value(pt_id(c));                         //for the grid point given by its indices "c"
        }                                                   //in case that this grid point it is undefined
                                                            //this function returns POS_VALUE or NEG_VALUE
                                                            //depending on its sign

        //#############################################
        // functions to access signs of level set values
        //##############################################

        static sign_type sign(value_type value) {                           //this function returns the sign
            return (sign_type(math::signbit(value))?NEG_SIGN:POS_SIGN);     //of a given level set value "value"
        }

        sign_type sign(size_type pt_id) const {     //this function returns the sign of a grid point given by
            if (pt_id==POS_PT) return POS_SIGN;     //its point-ID, and also checks if the point-ID is a positive or
            if (pt_id==NEG_PT) return NEG_SIGN;     //negative run
            //shfdhsfhdskjhgf assert(is_defined(pt_id));
            return sign_type(sign(value2(pt_id)));
        }

        sign_type sign2(size_type pt_id) const {    //this function returns the sign of a grid point given by
            return sign(value2(pt_id));             //its-point-ID
        }


        template <class CoordType>                  //returns the sign for a grid point given by its
        sign_type sign(const CoordType& c) const {  //indices "c"
            return sign(pt_id(c));
        }

        //################################################
        // functions to get information of a certain point
        // given by the point-ID number
        //################################################

        static bool is_defined(size_type r) {           //returns if the grid point given by the "pt_id"
            return (r<UNDEF_PT);                          //is a defined grid point
        }


        size_type active_pt_id(size_type pt_id) const {     //returns the active_pt_id for a grid point given by the
            if (is_defined(pt_id)) {                        //"pt_id". If it is an undefined (pt_id=POS_PT or NEG_PT)
                return active_pt_id2(pt_id);                       //or inactive point INACTIVE is returned
            } else {
                return INACTIVE;
            }
        }

        bool is_active(size_type pt_id) {           //returns if the grid point given by the "pt_id"
          return (active_pt_id(pt_id)!=INACTIVE);    //is an active grid point
        }


        size_type active_pt_id2(size_type pt_id) const {    //this function assumes that pt_id is a valid point-ID of

            //shfdhsfhdskjhgf assert(pt_id_offsets.size()!=0);
            //shfdhsfhdskjhgf assert(active_pt_id_offsets.size()!=0);

            typename std::vector<size_type>::const_iterator it= std::upper_bound(pt_id_offsets.begin(),pt_id_offsets.end(),pt_id)-1;

            //shfdhsfhdskjhgf assert(*it<=pt_id);

            int sub=it-pt_id_offsets.begin();

            size_type a=sub_levelsets[sub].active[pt_id-(*it)];

            if (a!=INACTIVE) a+=active_pt_id_offsets[sub];

            return a;
        }

        /*bool is_active(local_pt_id_type pt_id) const {             //returns if the grid point given by its ID "pt_id"
            return (active_pt_id(pt_id)!=INACTIVE);         //is an active grid point, using the function active_pt_id
        }

        bool is_active2(local_pt_id_type pt_id) const {            //returns if the grid point given by its ID "pt_id"
            return (active_pt_id2(pt_id)!=INACTIVE);        //is an active grid point, using the function active_pt_id
        }*/

        //#########################################################


        //void clear(const points_type & =points_type());    //this function deletes all defined grid points and
                                                //assigns to the whole grid the sign given by "sign"
                                                //(an undefined run of the specified sign is created containing
                                                //all grid points of the grid)
                                                //by default, all grid points are initialized to positive

        void initialize(const points_type & =points_type(), const allocation_type & = allocation_type());

        void finalize(int);

        allocation_type get_allocation() const {
            // allocation_type allocates the requried sizes to num_values and num_runs for all the sub_levelsets members
            // num_values[0] is to contain level set values, num_values[i] contains the start indices at the i-th dimension
            // num_runs[i] is to contain the run types at the i-th dimension

             allocation_type a;
             a.num_values=vec<index_type,D>(0);
             a.num_runs=vec<index_type,D>(0);
             for (typename sub_levelsets_type::size_type i=0; i!=sub_levelsets.size();++i ) {
                 allocation_type b=sub_levelsets[i].get_allocation();
                 a.num_values=Max(a.num_values, b.num_values);
                 a.num_runs=Max(a.num_runs, b.num_runs);
            }

             return (a*sub_levelsets.size());
             //return a;
        }

        template <class PointsType> levelset(const grid_type2&, const PointsType& point_defs);
        //this is a constructor, creating a new level set from a sorted list of index/level set value-pairs

        levelset(const levelset& l):    segmentation(l.segmentation),
                                        sub_levelsets(l.sub_levelsets),
                                        num_layers(l.num_layers),
                                        Grid(l.Grid),
                                        active_pt_id_offsets(l.active_pt_id_offsets),
                                        pt_id_offsets(l.pt_id_offsets)
                                        {

             //for (typename sub_levelsets_type::iterator it = sub_levelsets.begin(); it!=sub_levelsets.end(); ++it) {
             //    it->l=this;
             //}
            //sub_levelsets.update_levelset_ref(*this);
        }

        levelset(const grid_type2&);
        //initializes a levelset function without any defined points
        //if sign = POS_SIGN(default) one  undefined run of specified sign
        //is created containing all grid points

        levelset(const grid_type2&,  value_type, int, bool );
        //initializes a plain levelset with its normal pointing into the given axis direction

        /*const levelset & operator=(const levelset& l) {

             segmentation=l.segmentation;
             sub_levelsets=l.sub_levelsets;
             num_layers=l.num_layers;
             Grid=l.Grid;
             for (typename sub_levelsets_type::iterator it = sub_levelsets.begin(); it!=sub_levelsets.end(); ++it) {
                 it->l=this;
             }
             //shfdhsfhdskjhgf assert(sub_levelsets.size()==segmentation.size()+1);

             return *this;
        }*/



        template <class V> void push_back(int sub, const V& point, value_type distance) {
            //this function adds a new defined grid point to the level set function
            //the indices of the new grid point "indices" must be greater (lexiographically greater) than
            //the indices of the current "last" grid point
            sub_levelsets[sub].push_back(point,distance);
        }

        template <class V> void push_back_undefined(int sub, const V& point, size_type rt) {
            //this function sets the sign of an undefined run starting at indices "point"
            //the indices of the new grid point "indices" must be greater (lexiographically greater) than
            //the indices of the current "last" grid point
            sub_levelsets[sub].push_back_undefined(point,rt);
        }

        template <class V> int find_sub_levelset(const V& point) {
            return std::upper_bound(segmentation.begin(), segmentation.end(), point)-segmentation.begin();
        }

        vec<index_type, D> find_coordinate_of_pt(size_type pt) const {

            size_type rem_pt=pt;

            vec<index_type, D> tmp;

            //if (pt==num_pts()) return grid().max_point_index();

            //shfdhsfhdskjhgf assert(pt<num_pts());

             const int sub= (std::upper_bound(pt_id_offsets.begin()+1,pt_id_offsets.end(),pt)-(pt_id_offsets.begin()+1));

             pt-=pt_id_offsets[sub];

             //shfdhsfhdskjhgf assert(sub>=0);
             //shfdhsfhdskjhgf assert(sub<=segmentation.size());
             //shfdhsfhdskjhgf assert(segmentation.size()+1==pt_id_offsets.size());

             const sub_levelset_type & s=sub_levelsets[sub];

             for (int g=0;g<D;++g) {

                 size_type min=0;
                 size_type max=s.runtypes[g].size()-1;

                 while (!is_defined(s.runtypes[g][min])) ++min;
                 while (!is_defined(s.runtypes[g][max])) --max;

                 while (min!=max) {
                     size_type mid=(max+min+1)/2;
                     while (!is_defined(s.runtypes[g][mid])) ++mid;

                     if (s.runtypes[g][mid]<=pt) {
                         min=mid;
                     } else {
                         max=(max+min-1)/2;
                         while (!is_defined(s.runtypes[g][max])) --max;
                     }
                 }

                 //shfdhsfhdskjhgf assert(s.runtypes[g][min]<=pt);
                 if (min+1!=s.runtypes[g].size()) {
                     //shfdhsfhdskjhgf assert(s.runtypes[g][min+1]>pt);
                 }


                 //shfdhsfhdskjhgf assert(is_defined(s.runtypes[g][min]));

                 //shfdhsfhdskjhgf assert(pt>=s.runtypes[g][min]);

                 tmp[g]=pt-s.runtypes[g][min];

                 pt=std::upper_bound(s.start_indices[g].begin()+1,s.start_indices[g].end(),min)-(s.start_indices[g].begin()+1);

                 //shfdhsfhdskjhgf assert(s.start_indices[g][pt]<=min);
                 if (pt+1!=s.start_indices[g].size()) {
                     //shfdhsfhdskjhgf assert(s.start_indices[g][pt+1]>min);
                 }



                 //shfdhsfhdskjhgf assert(pt>=0);

                 tmp[g]+=s.GetRunStartCoord(g , pt, min);

                 //shfdhsfhdskjhgf assert(tmp[g]>=grid().min_point_index(g));
                 //shfdhsfhdskjhgf assert(tmp[g]<=grid().max_point_index(g));


             }

             const_iterator_runs it(*this, tmp);

             if (it.pt_id()!=rem_pt) {
                 std::cout << num_pts() << std::endl;
                 std::cout << rem_pt << std::endl;
                 std::cout << it.pt_id() << std::endl;
                 //shfdhsfhdskjhgf assert(it.pt_id()==rem_pt);
             }

             return tmp;
        }

        points_type get_new_segmentation() const {

            points_type tmp_segmentation;

            int n=1;    //TODO
            #ifdef _OPENMP
                n=omp_get_max_threads();
            #endif

            size_type n_pts=num_pts();  // number of defined grid points
            size_type sum=0;

            for (unsigned int g=0;g<static_cast<unsigned int>(n-1);++g) {
                sum+=n_pts/n+((n_pts%n)>g);
                //shfdhsfhdskjhgf assert(sum==(n_pts/n)*(g+1)+std::min(g+1,n_pts%n));
                if (sum!=n_pts) tmp_segmentation.push_back(find_coordinate_of_pt(sum));
            }

            if (n_pts==0) assert(tmp_segmentation.size()==0);

            //shfdhsfhdskjhgf assert(sum+n_pts/n==n_pts);
            ////shfdhsfhdskjhgf assert(tmp_segmentation.size()==n-1);

            return tmp_segmentation;
        }


        template <class V> void push_back(const V& point, value_type distance) {
            //this function adds a new defined grid point to the level set function
            //the indices of the new grid point "indices" must be greater (lexiographically greater) than the
            //the indices of the current "last" grid point
            push_back(find_sub_levelset(point), point,distance);
        }

        template <class V> void push_back_undefined(const V& point, size_type rt) {
            //this function sets the sign of an undefined run starting at indices "point"
            //the indices of the new grid point "indices" must be greater (lexiographically greater) than the
            //the indices of the current "last" grid point
            push_back_undefined(find_sub_levelset(point),point,rt);
        }

        size_type num_pts() const {
            //this function returns the number of defined grid points
            return pt_id_offsets.back()+sub_levelsets.back().num_pts();
        };

        size_type num_pts(int i) const {
            //this function returns the number of defined grid points
            return sub_levelsets[i].num_pts();
        };


        size_type num_active_pts() const {
            //this function returns the number of active grid points
             return active_pt_id_offsets.back()+sub_levelsets.back().num_active_pts();
        }

        size_type num_active_pts(int i) const {
            return sub_levelsets[i].num_active_pts();
        }

        void invert();      //this function inverts the level set function, by
                            //changing the signs of the level set values, and
                            //by interchanging positive and negative undefined runs

        void serialize();   //this function is used to serialize a segmented levelset

        void rebuild();     //this function is used to rebuild the level set function
                            //after time integration, after rebuilding the level set function
                            //will consist of two layers of defined grid points (number_of_layers()==2)

        template <class DataType>                          //this function does the same as the Rebuild()-function
        void rebuild(const DataType& Data);    //however, it allows the user to transfer user-data asigned to
                                                                                //active grid points, to the new active grid points
                                                                                //the size of the array Data must be equal to the number of active grid points
                                                                                //after the call of this functions the size of Data is again equal to
                                                                                //the number of active grid points after rebuilding
                                                                                //DataSize defines how many consecutive entries of the array Data
                                                                                //are assigned to each active grid point


        int number_of_layers() const { return num_layers;}    //this function returns the number of layers
                                                              //a level set function exists of at least

        /*void set_number_of_layers(int x) {                  //this function sets the number of layers
            //shfdhsfhdskjhgf assert(x<=num_layers);
            //shfdhsfhdskjhgf assert(x>0);
            num_layers=x;
        }*/

        void expand(int end_width, int start_width);        //this function expands the level set function
                                                            //"start_width" is the number of layers of the level set function
                                                            //which should be used, then the level set function is expanded
                                                            //where "end_width" defines the final number of layers
                                                            //HINT: start_width<=number_of_layers() must be fulfilled

        void expand(int end_width) {                        //this function expands the level set function
            expand(end_width, number_of_layers());          //to the given number of layers "end_width"
        }

    private:

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

        template<class BinaryType>
        void minmax(const levelset&, const BinaryType& min_or_max);

    public:
        void min(const levelset& l) {
            minmax(l, min_wrapper());
        }

        void max(const levelset& l) {
            minmax(l, max_wrapper());
        }



        void reduce(int width);         //this function reduces the level set function
                                        //to the given number of layers (>=1)

        void segment();                 //segments the levelset in the available sub_levelsets to balance correctly

        void prune();                 //removes all grid points, which do not have at least one opposite signed neighbour and distributes points across sub_levelsets to approximately balance the number of points across threads

        void add_voxel_corners();       //this function expands the level set functions to 5 layers, in addition it adds defined grid points, so that
                                        //each voxel (grid cell) which has at least one active grid point as corner, has all its corner grid points defined

        void add_voxel_corners2();

        template <class PointsType> void insert_points(const PointsType& point_defs); //this functions clears the level set function
                                                                                      //and setups the function from new
                                                                                      //from a sorted list "point_defs" of index/level set value-pairs

        //iterator declarations
        //iterators are used to process the data sequentially and to obtain linear scaling algorithms

        class const_iterator_base;

        class const_iterator_runs;  //this iterator is the main iterator, on which all other iterators base on
                                    //this iterator iterates over all runs and stops at every start index of each run
                                    //this run can be defined or undefined, if it is a defined run it corresponds to a single defined grid point


        class const_iterator_runs_offset;   //this iterator is derived from the const_iterator_runs class and
                                            //allows to iterate over the level set function with a given offset-index
                                            //with this iterator it is possible to create neighbor-iterators, which consist of a
                                            //const_iterator_runs and a const_iterator_runs_offset-iterator for each neighbor grid point
                                            //a movement of the neighbor-iterator corresponds to a synchronous movement of
                                            //all individual iterators. In this way linear scaling algorithms can be obtained


        template <class F, int NumN> class const_iterator_neighbor_filtered;    //this neighbor-iterator allows to access up to "NumN"-neighbors in all positive and negative axis directions
                                                                                //the filter "F" defines where the iterator for the center should stop

        template <class F> class const_iterator_cells_filtered;   //this iterator consist of 2^Dimensions const_iterator_runs_offset-iterators
                                                                  //one for each corner of a grid cell
                                                                  //all these 2^Dimensions iterators are synchronously moved
                                                                  //whenever the filter is fulfilled for one of
                                                                  //these iterators, the const_iterator_cells_filtered iterator stops

        class const_iterator_neighbor;      //this neighbor iterator consists of 2*Dimenions const_iterator_runs_offset for the neighbors and a
                                            //const_iterator_runs-iterator for the center
                                            //whenever one of these (2*Dimensions+1) iterators reach a defined grid point
                                            //the const_iterator_neighbor iterator stops
                                            //this iterator is internally used for rebuilding and expanding the level set function

    private:
         template <class I>
         void next(I &) const;

         template <class I>
         void previous(I &) const;

         template <class V, class I>
         void go_to_indices(const V&, I&) const;

         template <class V, class I>
         void go_to_indices(int, const V&, I&) const;

         template <class V, class I>
         void go_to_indices_sequential(const V&, I&) const;

    public:



    };


    //######################################
    // Definitions of constants
    //######################################

    template <class GridTraitsType, class LevelSetTraitsType> const typename levelset<GridTraitsType, LevelSetTraitsType>::size_type levelset<GridTraitsType, LevelSetTraitsType>::SEGMENT_PT =  std::numeric_limits<size_type>::max()-(lvlset_omp_max_num_threads-1);
    template <class GridTraitsType, class LevelSetTraitsType> const typename levelset<GridTraitsType, LevelSetTraitsType>::size_type levelset<GridTraitsType, LevelSetTraitsType>::POS_PT     =  std::numeric_limits<size_type>::max()-(lvlset_omp_max_num_threads+0);
    template <class GridTraitsType, class LevelSetTraitsType> const typename levelset<GridTraitsType, LevelSetTraitsType>::size_type levelset<GridTraitsType, LevelSetTraitsType>::NEG_PT     =  std::numeric_limits<size_type>::max()-(lvlset_omp_max_num_threads+1);
    template <class GridTraitsType, class LevelSetTraitsType> const typename levelset<GridTraitsType, LevelSetTraitsType>::size_type levelset<GridTraitsType, LevelSetTraitsType>::UNDEF_PT   =  std::numeric_limits<size_type>::max()-(lvlset_omp_max_num_threads+2);

    template <class GridTraitsType, class LevelSetTraitsType> const typename levelset<GridTraitsType, LevelSetTraitsType>::value_type levelset<GridTraitsType, LevelSetTraitsType>::POS_VALUE=std::numeric_limits<value_type>::max();
    template <class GridTraitsType, class LevelSetTraitsType> const typename levelset<GridTraitsType, LevelSetTraitsType>::value_type levelset<GridTraitsType, LevelSetTraitsType>::NEG_VALUE=-std::numeric_limits<value_type>::max();

    template <class GridTraitsType, class LevelSetTraitsType> const typename levelset<GridTraitsType, LevelSetTraitsType>::size_type levelset<GridTraitsType, LevelSetTraitsType>::INACTIVE=std::numeric_limits<size_type>::max();

    template<class GridTraitsType, class LevelSetTraitsType> unsigned levelset<GridTraitsType, LevelSetTraitsType>::top_levelset_id=0;



    template <class GridTraitsType, class LevelSetTraitsType> template <class PointsType> void levelset<GridTraitsType, LevelSetTraitsType>::insert_points(
                            const PointsType& point_defs)
                            {
        //this functions clears the level set function
        //and setups the function from new
        //from a sorted list "point_defs" of index/level set value-pairs
        //NOTE: a valid list of points must include all grid points, which are connected by edges of the grid, which are intersected
        //by the surface

        //TODO: is not parallelized yet
        assert(!point_defs.empty());

        //std::vector<vec<index_type, D > > tmp_segmentation;
        points_type tmp_segmentation;

        initialize(tmp_segmentation);

        if (point_defs.front().first!=Grid.min_point_index()) {
            push_back_undefined(0, Grid.min_point_index(),(sign_type(sign(point_defs.front().second))==POS_SIGN)?POS_PT:NEG_PT);     //TODO
        }

        typename PointsType::const_iterator begin_pt_defs=point_defs.begin();//+starts[t_num];
        typename PointsType::const_iterator end_pt_defs=point_defs.end();//point_defs.begin()+starts[t_num+1];

        vec<index_type,D> min_index=begin_pt_defs->first;
        vec<sign_type,D> signs(sign_type(sign(begin_pt_defs->second)));


        typename PointsType::const_iterator it_pt_defs=begin_pt_defs;
        while (it_pt_defs!=end_pt_defs) {

            if (math::abs(it_pt_defs->second)>=std::numeric_limits<typename PointsType::value_type::second_type>::max()) {
                push_back_undefined(0,it_pt_defs->first,(sign_type(sign(it_pt_defs->second))==POS_SIGN)?POS_PT:NEG_PT);
            } else {
                push_back(0,it_pt_defs->first, it_pt_defs->second);
            }

            //determine signs for next undefined runs
            {
                bool tmb_b=false;
                for (int l=D-1;l>=0;--l) {
                    tmb_b=tmb_b || (it_pt_defs->first[l]>min_index[l]);
                    if (tmb_b) {
                        signs[l]=sign_type(sign(it_pt_defs->second));
                        min_index[l]=it_pt_defs->first[l];
                    }
                }
            }

            vec<index_type,D> index=it_pt_defs->first;
            vec<index_type,D> next_index;

            it_pt_defs++;

            if (it_pt_defs==end_pt_defs) {
                next_index=Grid.max_point_index();
                next_index[D-1]++;
            } else {
                next_index=it_pt_defs->first;
            }

            for (int q=0;q<D;q++) {
                vec<index_type,D> tmp=index;
                tmp[q]++;
                if (tmp[q]>Grid.max_point_index(q)) continue;
                for (int r=0;r<q;++r) tmp[r]=Grid.min_point_index(r);

                if (tmp>=next_index) break;

                push_back_undefined(0,tmp,(signs[q]==POS_SIGN)?POS_PT:NEG_PT);

            }
        }

        finalize(2);
    }


    template <class GridTraitsType, class LevelSetTraitsType> template <class PointsType> levelset<GridTraitsType, LevelSetTraitsType>::levelset(
                        const grid_type2& g,
                        const PointsType& point_defs
                        ) : Grid(g) {
       //this  constructor takes the grid definition "g" and a sorted list of
       //index/level set value-pairs "point_defs" to
       //initialize the level set function

        insert_points(point_defs);
    }


    template <class GridTraitsType, class LevelSetTraitsType> void levelset<GridTraitsType, LevelSetTraitsType>::initialize(const points_type & p, const allocation_type & a) {

        //assert(p.size()==a.size());

        segmentation=p;

        //sub_levelsets.clear();
        //sub_levelsets.resize(segmentation.size()+1, sub_levelset_type(this));
        sub_levelsets.initialize(segmentation.size()+1, Grid, a);

        for (typename sub_levelsets_type::size_type i=1;i<sub_levelsets.size();++i) {

            sub_levelset_type & s=sub_levelsets[i];

            s.push_back_undefined(Grid.min_point_index(), Grid.decrement_indices(segmentation[0]),SEGMENT_PT);

            for (typename sub_levelsets_type::size_type j=1;j<i;++j) {
                s.push_back_undefined(segmentation[j-1],Grid.decrement_indices(segmentation[j]),SEGMENT_PT+j);
            }
        }
    }

    template <class GridTraitsType, class LevelSetTraitsType> void levelset<GridTraitsType, LevelSetTraitsType>::finalize(int num_lay) {


        //finalize sub-levelsets
        for (typename sub_levelsets_type::size_type i=0;i<sub_levelsets.size()-1;++i) {

            sub_levelset_type & s=sub_levelsets[i];

            for (typename sub_levelsets_type::size_type j=i+1;j<segmentation.size();++j) {
                s.push_back_undefined(segmentation[j-1],Grid.decrement_indices(segmentation[j]),SEGMENT_PT+j);
            }
            s.push_back_undefined(segmentation.back(),Grid.max_point_index(),SEGMENT_PT+segmentation.size());

        }

        //set number of layers
        num_layers=num_lay;

        //calculate id-offsets
        pt_id_offsets.clear();
        active_pt_id_offsets.clear();
        pt_id_offsets.push_back(0);
        active_pt_id_offsets.push_back(0);
        for (typename sub_levelsets_type::size_type i=0;i<sub_levelsets.size()-1;++i) {
            pt_id_offsets.push_back(pt_id_offsets.back()+sub_levelsets[i].num_pts());
            active_pt_id_offsets.push_back(active_pt_id_offsets.back()+sub_levelsets[i].num_active_pts());
        }

        //shfdhsfhdskjhgf assert(active_pt_id_offsets.size()==sub_levelsets.size());
        //shfdhsfhdskjhgf assert(pt_id_offsets.size()==sub_levelsets.size());
    }




    template <class GridTraitsType, class LevelSetTraitsType> levelset<GridTraitsType, LevelSetTraitsType>::levelset(const grid_type2& g)  :  Grid(g)   {

        //initialize empty level set
        /*sub_levelsets.resize(segmentation.size()+1, sub_levelset_type(this));
        for (typename sub_levelsets_type::iterator it=sub_levelsets.begin();it!=sub_levelsets.end();++it) {
            //it->start_indices[D-1].push_back(0);
            //it->num_active_points=0;
        }*/
    }

    template <class GridTraitsType, class LevelSetTraitsType> levelset<GridTraitsType, LevelSetTraitsType>::levelset(const grid_type2& g, value_type c, int direction, bool is_direction_negative)  :  segmentation(points_type()), Grid(g)  {

        //TODO

        initialize();

        for (int v=0;v<D;++v) {
            if (v!=direction) {
                //shfdhsfhdskjhgf assert(grid().boundary_conditions(v)!=POS_INFINITE_BOUNDARY);
                //shfdhsfhdskjhgf assert(grid().boundary_conditions(v)!=NEG_INFINITE_BOUNDARY);
                //shfdhsfhdskjhgf assert(grid().boundary_conditions(v)!=INFINITE_BOUNDARY);
            } else {
                //shfdhsfhdskjhgf assert(grid().boundary_conditions(v)==INFINITE_BOUNDARY);
            }
        }

        //sub_levelsets.initialize(1,*this);

        //sub_levelsets[0].num_active_points=0;
        //sub_levelsets[0].start_indices[D-1].push_back(0);
        //runtypes[D-1].push_back((is_direction_negative)?POS_PT:NEG_PT);

        vec<index_type,D> min_=grid().min_point_index();
        vec<index_type,D> max_=grid().max_point_index();

        min_[direction]=static_cast<index_type>(std::floor(c));
        max_[direction]=static_cast<index_type>(std::ceil(c));

        if (min_[direction]==max_[direction]) --(min_[direction]);

        value_type min_value=min_[direction]-c;
        value_type max_value=max_[direction]-c;

        if (is_direction_negative) {
            min_value=-min_value;
            max_value=-max_value;
        }

        //shfdhsfhdskjhgf assert(std::signbit(min_value)!=std::signbit(max_value));

        --(min_[direction]);
        ++(max_[direction]);

        //shfdhsfhdskjhgf assert(min_[direction]+3==max_[direction]);

        vec<index_type,D> it=min_;
        while(true) {

            if (it[direction]==min_[direction]) {
                vec<index_type,D> t=it;
                t[direction]=grid().min_point_index(direction);
                sub_levelsets[0].push_back_undefined(t, (is_direction_negative)?POS_PT:NEG_PT);
            } else if (it[direction]==min_[direction]+1) {
                sub_levelsets[0].push_back(it, min_value);
            } else if (it[direction]==max_[direction]-1) {
                vec<index_type,D> t=it;
                ++t[direction];
                sub_levelsets[0].push_back(it, max_value);
            } else {
                //shfdhsfhdskjhgf assert(it[direction]==max_[direction]);
                sub_levelsets[0].push_back_undefined(it, (is_direction_negative)?NEG_PT:POS_PT);
            }

            if (it==max_) break;

            //increment

            for (int j=0;j<D;++j) {
                if (it[j]!=max_[j]) {
                    it[j]++;
                    break;
                } else {
                    it[j]=min_[j];
                }
            }
        }

        finalize(2);

        std::string err= misc::test(*this);     //check level set function
        if (err.size()) {                   //if inconsistent print out error message
            std::cout << "plain init  failed!" << std::endl;
            std::cout << err << std::endl;
            assert(0);
        }
    }

    template <class GridTraitsType, class LevelSetTraitsType> void levelset<GridTraitsType, LevelSetTraitsType>::segment(){
        //segments the levelset in the available sub_levelsets to balance correctly
        //create new Level-Set
        levelset<GridTraitsType, LevelSetTraitsType> old_lvlset(grid());
        swap(old_lvlset);

        initialize(old_lvlset.get_new_segmentation(), old_lvlset.get_allocation());

        #pragma omp parallel num_threads(sub_levelsets.size()) //use num_threads(sub_levelsets.size()) threads
        {
            int p=0;
            #ifdef _OPENMP
            p=omp_get_thread_num();
            #endif

            sub_levelset_type & s=sub_levelsets[p];

            const_iterator_runs it(old_lvlset, (p==0)?grid().min_point_index():segmentation[p-1]);
            vec<index_type, D> end_v=(p!=static_cast<int>(segmentation.size()))?segmentation[p]:grid().max_point_index();

            for (;it.start_indices()<end_v;it.next()) {
                if (it.is_defined()) {
                    s.push_back(it.start_indices(),it.value());
                } else {
                    s.push_back_undefined(it.start_indices(),(it.sign()==POS_SIGN)?POS_PT:NEG_PT);
                }
            }
        }

        finalize(std::min(old_lvlset.num_layers,static_cast<int>(2)));
    }

    template <class GridTraitsType, class LevelSetTraitsType> void levelset<GridTraitsType, LevelSetTraitsType>::prune() {
        //removes all grid points, which do not have at least one opposite signed neighbour and distributes points across sub_levelsets to approximately balance the number of points
        levelset<GridTraitsType, LevelSetTraitsType> old_lvlset(grid()); //create new Level-Set
        swap(old_lvlset);

        initialize(old_lvlset.get_new_segmentation(), old_lvlset.get_allocation());

        #pragma omp parallel num_threads(sub_levelsets.size()) //use num_threads(sub_levelsets.size()) threads
        {
            int p=0;
            #ifdef _OPENMP
            p=omp_get_thread_num();
            #endif

            sub_levelset_type & s=sub_levelsets[p];

            const_iterator_neighbor_filtered<filter_all,1> it(old_lvlset, filter_all(), (p==0)?grid().min_point_index():segmentation[p-1]);
            vec<index_type, D> end_v=(p!=static_cast<int>(segmentation.size()))?segmentation[p]:grid().increment_indices(grid().max_point_index());

            for (;it.indices()<end_v;it.next()) {

                if (it.center().is_defined()) {
                    int i=0;
                    sign_type sgn=it.center().sign2();
                    for(;i<2*D;i++) {
                        if (it.neighbor(i).sign()!=sgn) break;
                    }
                    if (i!=2*D) {
                        s.push_back(it.indices(),it.center().value());
                    } else {
                        s.push_back_undefined(it.indices(),(sgn==POS_SIGN)?POS_PT:NEG_PT);
                    }
                } else {
                    s.push_back_undefined(it.indices(),(it.center().sign()==POS_SIGN)?POS_PT:NEG_PT);
                }
            }
        }

        finalize(std::min(old_lvlset.num_layers,static_cast<int>(2)));
    }

    template <class GridTraitsType, class LevelSetTraitsType>
    template <class DataType>
    void levelset<GridTraitsType, LevelSetTraitsType>::rebuild(const DataType& data) {

        //create new Level-Set
        levelset<GridTraitsType, LevelSetTraitsType> old_lvlset(grid());
        swap(old_lvlset);

        initialize(old_lvlset.get_new_segmentation(), old_lvlset.get_allocation()*(2./old_lvlset.num_layers));

        #pragma omp parallel num_threads(sub_levelsets.size()) //use num_threads(sub_levelsets.size()) threads
        {
            std::vector<size_type> data_transfer;

            {
                const size_type tmp=old_lvlset.num_active_pts()*(allocation_factor/sub_levelsets.size());     //TODO
                data_transfer.reserve(tmp);
            }

            int p=0;
            #ifdef _OPENMP
            p=omp_get_thread_num();
            #endif


            sub_levelset_type & s=sub_levelsets[p];

            const_iterator_neighbor it(old_lvlset, (p==0)?grid().min_point_index():segmentation[p-1]);
            vec<index_type, D> end_v=(p!=static_cast<int>(segmentation.size()))?segmentation[p]:grid().increment_indices(grid().max_point_index());

            for (;it.start_indices()!=end_v;it.next()) {

                if (it.center().is_active()) {  //if the center is an active grid point

                    int k=0;
                    for(;k<2*D;k++) if (it.neighbor(k).sign()!=it.center().sign()) break;

                    if (k!=2*D) {       //if there is at least one neighbor of opposite sign

                        if (it.center().value2()>0.5) {
                            int j=0;
                            for (;j<2*D;j++) {
                                if (it.neighbor(j).is_active()) if (it.neighbor(j).value2()<-0.5) break;
                            }

                            if (j==2*D) {
                                s.push_back(it.start_indices(),it.center().value2());
                            } else {   //if there is at least one active grid point, which is < -0.5
                                s.push_back(it.start_indices(),0.5);
                                data_transfer.push_back(it.center().active_pt_id2());
                                //std::copy(old_data.begin()+it.center().active_pt_id2()*data_size,old_data.begin()+it.center().active_pt_id2()*data_size+data_size,std::back_inserter(data));
                            }
                        } else if (it.center().value2()<-0.5) {
                            int j=0;
                            for (;j<2*D;j++) {
                                if (it.neighbor(j).is_active()) if (it.neighbor(j).value2()>0.5) break;
                            }

                            if (j==2*D) {
                                s.push_back(it.start_indices(),it.center().value2());
                            } else {    //if there is at least one active grid point, which is > 0.5
                                s.push_back(it.start_indices(),-0.5);
                                data_transfer.push_back(it.center().active_pt_id2());
                                //std::copy(old_data.begin()+it.center().active_pt_id2()*data_size,old_data.begin()+it.center().active_pt_id2()*data_size+data_size,std::back_inserter(data));
                            }
                        } else {
                            s.push_back(it.start_indices(),it.center().value2());
                            data_transfer.push_back(it.center().active_pt_id2());
                            //std::copy(old_data.begin()+it.center().active_pt_id2()*data_size,old_data.begin()+it.center().active_pt_id2()*data_size+data_size,std::back_inserter(data));
                        }
                    } else {
                        s.push_back_undefined(it.start_indices(),(it.center().sign()==POS_SIGN)?POS_PT:NEG_PT);
                    }

                } else {    //if the center is not an inactive grid point

                    int n_pt=-1;
                    if (it.center().sign()==POS_SIGN) {
                        value_type distance=POS_VALUE;
                        for(int i=0;i<2*D;i++) {
                            if (it.neighbor(i).is_active() && (it.neighbor(i).sign()==NEG_SIGN)) {
                                if (distance>it.neighbor(i).value()+1.) {
                                    distance=copysign(it.neighbor(i).value()+1.,1.);
                                    if (distance<=0.5) n_pt=i;
                                }
                            }
                        }
                        if (n_pt>=0) {
                            //size_type tmp=it.neighbor(n_pt).active_pt_id2()*data_size;
                            //std::copy(old_data.begin()+tmp,old_data.begin()+(tmp+data_size),std::back_inserter(data));
                            data_transfer.push_back(it.neighbor(n_pt).active_pt_id2());
                            s.push_back(it.start_indices(),distance);
                        } else if (distance<=1.) {
                            s.push_back(it.start_indices(),distance);
                        } else {
                            s.push_back_undefined(it.start_indices(),POS_PT);
                        }

                    } else {
                        value_type distance=NEG_VALUE;
                        for(int i=0;i<2*D;i++) {
                            if (it.neighbor(i).is_active() && (it.neighbor(i).sign()==POS_SIGN)) {
                                if (distance<it.neighbor(i).value()-1.) {
                                    distance=copysign(it.neighbor(i).value()-1.,-1.);
                                    if (distance>=-0.5) n_pt=i;
                                }
                            }
                        }
                        if (n_pt>=0) {
                            //size_type tmp=it.neighbor(n_pt).active_pt_id2()*data_size;
                            data_transfer.push_back(it.neighbor(n_pt).active_pt_id2());
                            //std::copy(old_data.begin()+tmp,old_data.begin()+(tmp+data_size),std::back_inserter(data));
                            s.push_back(it.start_indices(),distance);
                        } else if (distance>=-1.) {
                            s.push_back(it.start_indices(),distance);
                        } else {
                            s.push_back_undefined(it.start_indices(),NEG_PT);
                        }
                    }
                }
            }

      #pragma omp barrier //wait until all other threads in section reach the same point.
      #pragma omp single //section of code that must be run by a single available thread.
            {
                finalize(2);
                data.set_size(num_active_pts());
            }

            for (size_type i=0;i<s.num_active_pts();++i) data.copy_data(i+active_pt_id_offsets[p],data_transfer[i]);
        }

    }

    namespace {
        class tmp_data {
        public:
            template <class I> void set_size(I i) const {}
            template <class I> void copy_data(I i, I j) const {}
        };
    }


    template <class GridTraitsType, class LevelSetTraitsType> void levelset<GridTraitsType, LevelSetTraitsType>::rebuild() {
        //std::vector<int> dummy;
        rebuild(tmp_data());
    }


    template <class GridTraitsType, class LevelSetTraitsType>
    template <class BinaryType>
    void levelset<GridTraitsType, LevelSetTraitsType>::minmax(const levelset<GridTraitsType, LevelSetTraitsType>& lB, const BinaryType& min_or_max) {
        //this function is called by the min and the max function defined below
        //this function should not be called directly by the user,
        //      and is therefore within an anonymous namespace
        //depending on which binary functor is passed to the library (min or max)
        //      this function calculates the minimum or maximum of two level set functions
        //      with an optimal complexity O(lA.num_pts()+lB.num_pts())


        typedef levelset<GridTraitsType, LevelSetTraitsType> LevelSetType;

        LevelSetType lA(grid());

        swap(lA);

        //shfdhsfhdskjhgf assert(lA.number_of_layers()>0);
        //shfdhsfhdskjhgf assert(lB.number_of_layers()>0);

        initialize(lA.get_new_segmentation(), lA.get_allocation()+lB.get_allocation()); //TODO

        #pragma omp parallel num_threads(sub_levelsets.size()) //use num_threads(sub_levelsets.size()) threads
        {
            //#pragma omp for schedule(static, 1)
            //for (int p=0;p<static_cast<int>(sub_levelsets.size());++p) {
            {
                int p=0;
                #ifdef _OPENMP
                p=omp_get_thread_num();
                #endif


                sub_levelset_type & s=sub_levelsets[p];

                vec<index_type, D> pos=(p==0)?grid().min_point_index():segmentation[p-1];
                vec<index_type, D> end_v=(p!=static_cast<int>(segmentation.size()))?segmentation[p]:grid().increment_indices(grid().max_point_index());

                const_iterator_runs itA(lA,pos);
                const_iterator_runs itB(lB,pos);

                while (pos<end_v) {

                    typename LevelSetType::value_type d=min_or_max(itA.value(),itB.value());
//                  typename LevelSetType::value_type d;

//                    if (std::abs(itA.value()-itB.value())<0.001) d=POS_VALUE;

//                  if ((itA.value()-itB.value())<-0.001) {
//                    d=POS_VALUE;
//                  } else if (std::abs(itA.value())<std::abs(itB.value())) {
//                    d=itA.value();
//                  } else {
//                    d=itB.value();
//                  }

//                  min_or_max(itA.value(),itB.value());

//                    if ((itA.value()<0)&&(itB.value()<0)) {
//                      std::cout << "itA.value(): " << itA.value() << ", itB.value(): " << itB.value() << ", d: " << d << std::endl;
//                    }

                    if (math::abs(d)<std::numeric_limits<typename LevelSetType::value_type>::max()) {
                        s.push_back(pos, d);
                    } else {
                        s.push_back_undefined(pos, (sign(d)==POS_SIGN)?POS_PT:NEG_PT);
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
                    pos=std::max(itA.start_indices(),itB.start_indices());
                }
            }
        }

        finalize(std::min(lA.number_of_layers(), lB.number_of_layers()));
    }




    template <class GridTraitsType, class LevelSetTraitsType> void levelset<GridTraitsType, LevelSetTraitsType>::reduce(int width) {
        //this function reduces the level set function
        //to the given number of layers (>=1)

        levelset<GridTraitsType, LevelSetTraitsType> old_lvlset(grid());

        swap(old_lvlset);

        initialize(old_lvlset.get_new_segmentation(), old_lvlset.get_allocation()*std::min(1.,(static_cast<double>(width)/old_lvlset.number_of_layers())));
        //initialize();



        const value_type limit=width*0.5;

        //const int  factor=std::min(width, old_lvlset.number_of_layers())

        #pragma omp parallel num_threads(sub_levelsets.size()) //use num_threads(sub_levelsets.size()) threads
        {
            //#pragma omp for schedule(static, 1)
            //for (int p=0;p<static_cast<int>(sub_levelsets.size());++p) {
            {
                int p=0;
                #ifdef _OPENMP
                p=omp_get_thread_num();
                #endif


                sub_levelset_type & s=sub_levelsets[p];


                vec<index_type, D> begin_v=(p==0)?grid().min_point_index():segmentation[p-1];
                vec<index_type, D> end_v=(p!=static_cast<int>(segmentation.size()))?segmentation[p]:grid().increment_indices(grid().max_point_index());

                for (const_iterator_runs it(old_lvlset,begin_v);it.start_indices()<end_v;it.next()) {
                //for (const_iterator_runs it(old_lvlset);!it.is_finished();it.next()) {
                    if (it.is_defined()) {
                        value_type d=it.value2();
                        if (math::abs(d)<=limit) {
                            s.push_back(it.start_indices(),d);
                            continue;
                        }
                    }
                    s.push_back_undefined(it.start_indices(),(it.sign()==POS_SIGN)?POS_PT:NEG_PT);
                }
            }
        }

        finalize(std::min(width, old_lvlset.number_of_layers()));
        segment();
        /*std::string err= misc::test(*this);     //check level set function      //TODO
        if (err.size()) {                   //if inconsistent print out error message
            std::cout << "reduce failed!" << std::endl;
            std::cout << err << std::endl;
            //shfdhsfhdskjhgf assert(0);
        }*/

    }


    template <class GridTraitsType, class LevelSetTraitsType> void levelset<GridTraitsType, LevelSetTraitsType>::serialize() {

        if (sub_levelsets.size()>1) {

            levelset<GridTraitsType, LevelSetTraitsType> old_lvlset(grid());        //TODO

            swap(old_lvlset);

            initialize();

            sub_levelset_type & s=sub_levelsets[0];

            for (const_iterator_runs it(old_lvlset);!it.is_finished();it.next()) {
                if (it.is_defined()) {
                    s.push_back(it.start_indices(),it.value2());
                } else {
                    s.push_back_undefined(it.start_indices(),(it.sign()==POS_SIGN)?POS_PT:NEG_PT);
                }
            }

            finalize(old_lvlset.number_of_layers());
        }
    }


    template <class GridTraitsType, class LevelSetTraitsType> void levelset<GridTraitsType, LevelSetTraitsType>::expand(int end_width, int start_width) {

        //shfdhsfhdskjhgf assert(start_width<=number_of_layers());

        const value_type total_limit=end_width*0.5;

        levelset<GridTraitsType, LevelSetTraitsType> old_lvlset(grid());

        int num_required_cycles=(1+end_width-start_width)/2;

        for (int a=0; a< num_required_cycles; ++a) {

            const int num_layers=end_width+2*(1+a-num_required_cycles);
            const value_type limit = num_layers*value_type(0.5);

            swap(old_lvlset);

            //double tmp;
            //tmp=-my::time::GetTime();
            initialize(old_lvlset.get_new_segmentation(), old_lvlset.get_allocation()*(static_cast<double>(num_layers)/old_lvlset.num_layers));
            //tmp+=my::time::GetTime();
            //std::cout << "time1 = " << tmp;

            #pragma omp parallel num_threads(sub_levelsets.size()) //use num_threads(sub_levelsets.size()) threads
            {
                //#pragma omp for schedule(static, 1)
                //for (int p=0;p<static_cast<int>(sub_levelsets.size());++p) {
                {
                    int p=0;
                    #ifdef _OPENMP
                    p=omp_get_thread_num();
                    #endif


                    sub_levelset_type & s=sub_levelsets[p];

                    vec<index_type, D> begin_v=(p==0)?grid().min_point_index():segmentation[p-1];
                    vec<index_type, D> end_v=(p!=static_cast<int>(segmentation.size()))?segmentation[p]:grid().increment_indices(grid().max_point_index());

                    for (const_iterator_neighbor it(old_lvlset, begin_v);it.start_indices()!=end_v;it.next()) {

                        if (math::abs(it.center().value())<=total_limit) {
                            s.push_back(it.start_indices(),it.center().value());
                        } else {
                            if (it.center().sign()==POS_SIGN) {
                                value_type distance=POS_VALUE;
                                for (int i=0;i<2*D;i++) {
                                    //if (it.neighbor(i).sign()>0) //shfdhsfhdskjhgf assert(it.neighbor(i).value()!=0);
                                    distance=std::min(distance,it.neighbor(i).value()+value_type(1));
                                }
                                if (distance<=limit) {
                                    s.push_back(it.start_indices(),distance);
                                } else {
                                    s.push_back_undefined(it.start_indices(),POS_PT);
                                }
                            } else {
                                value_type distance=NEG_VALUE;
                                for (int i=0;i<2*D;i++) {
                                    //if (it.neighbor(i).sign()<0) //shfdhsfhdskjhgf assert(it.neighbor(i).value()!=0);
                                    distance=std::max(distance,it.neighbor(i).value()-value_type(1));
                                }
                                if (distance>=-limit) {
                                    s.push_back(it.start_indices(),distance);
                                } else {
                                    s.push_back_undefined(it.start_indices(),NEG_PT);
                                }
                            }
                        }
                    }
                }
            }

            //tmp=-my::time::GetTime();
            finalize(std::max(start_width, end_width+2*(1+a-num_required_cycles)));
            //tmp+=my::time::GetTime();
            //std::cout << "     time2 = " << tmp << std::endl;

        }
        //setting up the segmentation correctly
        segment();
    }


    /*
    //this function expands the level set functions to 5 layers, in addition it adds defined grid points, so that
    //each voxel (grid cell) which has at least one active grid point as corner, has all its corner grid points defined
    */

    template <class GridTraitsType, class LevelSetTraitsType> void levelset<GridTraitsType, LevelSetTraitsType>::add_voxel_corners() {



        expand(5);

        //shfdhsfhdskjhgf assert(sub_levelsets.size()==active_pt_id_offsets.size());

        //create new Level-Set
        levelset<GridTraitsType, LevelSetTraitsType> old_lvlset(grid());
        swap(old_lvlset);

        //std::cout << "alt: " << old_lvlset.num_pts() << std::endl;

        //shfdhsfhdskjhgf assert(old_lvlset.sub_levelsets.size()==old_lvlset.active_pt_id_offsets.size());

        initialize();

        sub_levelset_type & s=sub_levelsets[0];       //TODO

        std::vector<vec<index_type,D> > new_pts;

        const_iterator_runs it_mid(old_lvlset);

        std::vector<const_iterator_runs_offset> it_corners;
        for (int i=0;i<(1<<D);i++) {
            vec<index_type,D> offset;
            for (int j=0;j<D;j++) offset[j]=((i & (1<<j)) == 0)?-1:1;
            it_corners.push_back(const_iterator_runs_offset(old_lvlset, offset));
        }

        std::vector<const_iterator_runs_offset> ix;
        for (int m=0;m<2*D;m++) {
            vec<index_type,D> tv(index_type(0));
            if (m<D) tv[m]++; else tv[m-D]--;
            ix.push_back(const_iterator_runs_offset(old_lvlset, tv));
        }

        vec<index_type,D> coords=Grid.min_point_index();

        while (!it_mid.is_finished()) {

            if (it_mid.is_defined()) {
                s.push_back(coords, it_mid.value2());
            } else {
                //check corners
                int k=0;
                for (;k<(1<<D);++k) {
                    if (it_corners[k].is_active()) break;
                }
                if (k!=(1<<D)) {
                    //insert new point
                    value_type distance;
                    for (int m=0;m<2*D;++m) ix[m].go_to_indices_sequential(coords);
                    if (it_mid.sign()==POS_SIGN) {
                        distance=POS_VALUE;
                        for (int m=0;m<2*D;++m) distance=std::min(distance, ix[m].value()+value_type(1));
                    } else {
                        distance=NEG_VALUE;
                        for (int m=0;m<2*D;++m) distance=std::max(distance, ix[m].value()-value_type(1));
                    }
                    s.push_back(coords, distance);
                } else {
                    s.push_back_undefined(coords,  (it_mid.sign()==POS_SIGN)?POS_PT:NEG_PT);
                }
            }

            std::bitset<(1<<D)+1> increment;
            increment.set(1<<D);

            coords=it_mid.end_indices();
            for (int i=0;i<(1<<D);i++) {
                switch (compare(coords, it_corners[i].end_indices())) {
                    case 1:
                        coords=it_corners[i].end_indices();
                        increment.reset();
                    case 0:
                        increment.set(i);
                }
            }

            for (int i=0;i<(1<<D);i++) if (increment.test(i)) it_corners[i].next();
            if (increment.test(1<<D)) it_mid.next();

            coords=old_lvlset.grid().increment_indices(coords);


           /* vec<index_type,D> tmp=it_mid.end_indices();
            for (int i=0;i<(1<<D);i++) {
                if (it_corners[i].end_indices()<tmp) tmp=it_corners[i].end_indices();
            }

            //go to next position
            if (tmp==it_mid.end_indices()) {
                it_mid.next();
                coords=it_mid.start_indices();
            }
            for (int i=0;i<(1<<D);i++) {
                if (tmp==it_corners[i].end_indices()) {
                    it_corners[i].next();
                    coords=it_corners[i].start_indices();
                }
            }*/

        }

        finalize(old_lvlset.num_layers);

        //std::cout << "neu: " << num_pts() << std::endl;

    }


    template <class GridTraitsType, class LevelSetTraitsType> void levelset<GridTraitsType, LevelSetTraitsType>::add_voxel_corners2() {

        //shfdhsfhdskjhgf assert(start_width<=number_of_layers());

        expand(5);

        levelset<GridTraitsType, LevelSetTraitsType> old_lvlset(Grid.grid_traits());

        swap(old_lvlset);
        initialize(old_lvlset.get_new_segmentation());

        //std::cout << "alt: " << old_lvlset.num_pts() << std::endl;

        #pragma omp parallel num_threads(sub_levelsets.size()) //use num_threads(sub_levelsets.size()) threads
        {
            #pragma omp for schedule(static, 1) // parallelization - Iterations divided into chunks of size 1. Each chunk is assigned to a thread
            for (int p=0;p<static_cast<int>(sub_levelsets.size());++p) {

                sub_levelset_type & s=sub_levelsets[p];
                const_iterator_neighbor it(old_lvlset, (p==0)?grid().min_point_index():segmentation[p-1]);
                vec<index_type, D> end_v=(p!=static_cast<int>(segmentation.size()))?segmentation[p]:grid().increment_indices(grid().max_point_index());

                for (;it.start_indices()!=end_v;it.next()) {

                    if (it.center().is_defined()) {
                        s.push_back(it.start_indices(),it.center().value());
                    } else {
                        if (it.center().sign()==POS_SIGN) {
                            value_type distance=POS_VALUE;
                            std::bitset<D> x;
                            for (int i=0;i<2*D;i++) {
                                if (it.neighbor(i).is_defined()) {
                                    if (it.neighbor(i).value2()<=2.5) x.set(i%D);
                                //if (it.neighbor(i).sign()>0) //shfdhsfhdskjhgf assert(it.neighbor(i).value()!=0);
                                    distance=std::min(distance,it.neighbor(i).value2()+value_type(1));
                                }
                            }
                            if ((distance<=3.5) && (x.count()==D)) {
                                s.push_back(it.start_indices(),distance);
                            } else {
                                s.push_back_undefined(it.start_indices(),POS_PT);
                            }
                        } else {
                            value_type distance=NEG_VALUE;
                            std::bitset<D> x;
                            for (int i=0;i<2*D;i++) {
                                if (it.neighbor(i).is_defined()) {
                                    if (it.neighbor(i).value2()>=-2.5) x.set(i%D);
                                //if (it.neighbor(i).sign()<0) //shfdhsfhdskjhgf assert(it.neighbor(i).value()!=0);
                                    distance=std::max(distance,it.neighbor(i).value2()-value_type(1));
                                }
                            }
                             if ((distance>= -3.5) && (x.count()==D)) {
                                s.push_back(it.start_indices(),distance);
                            } else {
                                s.push_back_undefined(it.start_indices(),NEG_PT);
                            }
                        }
                    }
                }
            }
        }

        finalize(old_lvlset.number_of_layers());

        //std::cout << "neu: " << num_pts() << std::endl;
    }

    template <class GridTraitsType, class LevelSetTraitsType> template<class F, int NumN> class levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_neighbor_filtered {
        //this neighbor-iterator allows to access up to "NumN"-neighbors in all positive and negative axis directions
        //the filter "F" defines where the iterator for the center should stop

        const_iterator_runs it_mid;                                 //the iterator for the center
        std::vector<const_iterator_runs_offset> it_neighbors;       //the neighbor iterators
        const F filter;                                             //the filter for the center grid point
        const levelset<GridTraitsType, LevelSetTraitsType>& l;      //reference to the level set function

    public:

        typedef levelset<GridTraitsType, LevelSetTraitsType> levelset_type;

        const levelset_type & get_levelset() const {
            return l;
        }

        const_iterator_neighbor_filtered(const levelset<GridTraitsType, LevelSetTraitsType>& lx, const F& f=F()) : it_mid(lx), filter(f), l(lx) {
            //constructor, initializes the iterator and sets to the first position fulfilling the "filter"-criterion

            while ((!filter(it_mid)) && (!it_mid.is_finished())) {
                it_mid.next();
            }

            if(it_mid.is_finished()) return;

            it_neighbors.reserve(2*D*NumN);

            for (int i=0;i<2*D;i++) {
                vec<index_type,D> tv(index_type(0));
                for (int j=0;j<NumN;j++) {
                    if (i<D) tv[i]++; else tv[i-D]--;
                    it_neighbors.push_back(const_iterator_runs_offset(lx, tv,it_mid.start_indices()));
                }
            }
        }

        template <class V>
        const_iterator_neighbor_filtered(const levelset<GridTraitsType, LevelSetTraitsType>& lx, const F& f, const V& v) : it_mid(lx,v), filter(f), l(lx) {
            //constructor, initializes the iterator and sets to the first position fulfilling the "filter"-criterion

            while ((!filter(it_mid)) && (!it_mid.is_finished())) {
                it_mid.next();
            }

            if(it_mid.is_finished()) return;

            it_neighbors.reserve(2*D*NumN);

            for (int i=0;i<2*D;i++) {
                vec<index_type,D> tv(index_type(0));
                for (int j=0;j<NumN;j++) {
                    if (i<D) tv[i]++; else tv[i-D]--;
                    it_neighbors.push_back(const_iterator_runs_offset(lx, tv,it_mid.start_indices()));
                }
            }
        }


        void next() {   //move to the next run fulfilling "filter"
            do {
                it_mid.next();
            } while ((!filter(it_mid)) && (!it_mid.is_finished()));
            if (!it_mid.is_finished()) {
                for (int i=0;i<2*D*NumN;i++) it_neighbors[i].go_to_indices_sequential(it_mid.start_indices());
            }
        }


        const const_iterator_runs_offset& neighbor(int direction, int item=0) const {
            //access to a certain neighbor grid point
            return it_neighbors[NumN*direction+item];
        }

        const const_iterator_runs& center() const {     //access to the center grid point
            return it_mid;
        }

        bool is_finished() const {                  //returns true if iterator is at the end
            return center().is_finished();
        }

        vec<index_type,D> indices() const { //returns the indices of the center
            return it_mid.start_indices();
        }

        index_type indices(int dir) const {         //returns the index of the center for the given axis direction
            return it_mid.start_indices(dir);
        }

        value_type gradient(int dir) const {
            //returns the derivation in respect to the real coordinates for the given axis direction
            //this function requires that the level set function consists at least of 3 layers of grid points,
            //which means that all neighbor grid points of active grid points are defined

            const value_type pos  =  l.Grid.grid_position_of_local_index(dir,indices(dir));
            const value_type d_p  =  l.Grid.grid_position_of_global_index(dir,indices(dir)+1)-pos;
            const value_type d_n  =  l.Grid.grid_position_of_global_index(dir,indices(dir)-1)-pos;

            const value_type phi_0=center().value();
            const value_type phi_p=neighbor(dir,0).value();
            const value_type phi_n=neighbor(dir+D,0).value();

            return ((d_p/d_n)*(phi_p-phi_0)-(d_n/d_p)*(phi_n-phi_0))/(d_n-d_p);
        }

        vec<value_type,D>  gradient() const {
            vec<value_type, D> tmp;
            for (int i=0;i<D;++i) tmp[i]=gradient(i);
            return tmp;
        }

        value_type gradient2(int dir) const {       //returns the derivation of the level set function in respect to the index for the given axis direction
            const value_type phi_p=neighbor(dir,0).value();
            const value_type phi_n=neighbor(dir+D,0).value();

            return (phi_p-phi_n)/2.;
        }

        vec<value_type,D>  gradient2() const {
            vec<value_type, D> tmp;
            for (int i=0;i<D;++i) tmp[i]=gradient2(i);
            return tmp;
        }

    };

    template <class GridTraitsType, class LevelSetTraitsType> class levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_base {
    protected:
        const levelset<GridTraitsType, LevelSetTraitsType> &l;
        vec<size_type,D+1> start_indices_pos;
        vec<size_type,D> run_type_pos;
        vec<index_type,D> start_run_abs_coords;
        vec<index_type,D> end_run_abs_coords;
        vec<index_type,D> abs_coords;
        vec<index_type,D> end_abs_coords;
        int r_level;
        int s_level;
        int sub;

        void go_up_BA() {
            ++r_level;
            //shfdhsfhdskjhgf assert(r_level==s_level);
        }

        void go_up_AB() {
            //shfdhsfhdskjhgf assert((r_level==s_level) || (s_level==r_level+1));
            s_level=r_level+1;
            abs_coords[r_level]=l.grid().min_point_index(r_level);
            end_abs_coords[r_level]=l.grid().max_point_index(r_level);
            //shfdhsfhdskjhgf assert(s_level==r_level+1);
        }


    public:

        void print() const {            //TODO
            std::cout << "start_indices_pos: " << start_indices_pos << std::endl;
            std::cout << "run_type_pos: " << run_type_pos << std::endl;
            std::cout << "start_run_abs_coords: " << start_run_abs_coords << std::endl;
            std::cout << "end_run_abs_coords: " << end_run_abs_coords << std::endl;
            std::cout << "abs_coords: " << abs_coords << std::endl;
            std::cout << "end_abs_coords: " << end_abs_coords << std::endl;
            std::cout << "r_level: " << r_level << std::endl;
            std::cout << "s_level: " << s_level << std::endl;
            std::cout << "sub: " << sub << std::endl;
            std::cout << "value: " << value() << std::endl;
            std::cout << "sign: " << sign() << std::endl;
            std::cout << "active: " << is_active() << std::endl;
            std::cout << "current_segment: " << current_segment() << std::endl;

        }



        typedef levelset<GridTraitsType, LevelSetTraitsType> levelset_type;

        const_iterator_base(const levelset& lx):l(lx),abs_coords(l.grid().min_point_index()), end_abs_coords(l.grid().max_point_index()), r_level(D), s_level(D), sub(0) {
            start_indices_pos[D]=0;
        }

        const levelset_type & get_levelset() const {
            return l;
        }

        const const_iterator_base& operator=(const const_iterator_base& it) {
            //copy assignment
            start_indices_pos=it.start_indices_pos;
            run_type_pos=it.run_type_pos;
            start_run_abs_coords=it.start_run_abs_coords;
            end_run_abs_coords=it.end_run_abs_coords;
            abs_coords=it.abs_coords;
            end_abs_coords=it.end_abs_coords;
            s_level=it.s_level;
            r_level=it.r_level;
            sub=it.sub;
            return *this;
        }

        bool is_finished() const {
            //returns true if iterator reached the end
            return (r_level==D);
        }

        index_type  start_indices(int dir) const {
            //returns the start index of a run for a certain axis direction
            return abs_coords[dir];
        }

        const vec<index_type,D> & start_indices() const {
            return abs_coords;
        }

        index_type  end_indices(int dir) const {
            //returns the end index of a run for a certain axis direction
            return end_abs_coords[dir];
        }


        const vec<index_type,D> & end_indices() const {
            return end_abs_coords;
        }

        size_type run_type_position() const {
            if (s_level==0) {
                //shfdhsfhdskjhgf assert(s_level==r_level);
                return start_indices_pos[0];
            } else {
                //shfdhsfhdskjhgf assert(s_level>r_level);
                return run_type_pos[r_level];
            }
        }

        int get_level() const {
            return s_level;
        }

        /*int get_sub() const {
            return sub;
        }*/

        size_type current_run_code() const {
            return start_indices_pos[r_level];
        }

        size_type current_segment() const {
            const size_type r=current_run_code();
            //shfdhsfhdskjhgf assert(POS_PT<SEGMENT_PT);
            //shfdhsfhdskjhgf assert(NEG_PT<SEGMENT_PT);
            if (r>=SEGMENT_PT) {
                return r-SEGMENT_PT;
            } else {
                return sub;
            }
        }

        size_type get_segment_num() const {
            return sub;
        }

        size_type pt_id() const {
            //returns the pt_id if it is a defined run
            //returns POS_PT or NEG_PT for undefined runs
            size_type tmp=current_run_code();
            if (is_defined()) tmp+=l.pt_id_offsets[sub];
            return tmp;
        }

        bool is_defined() const {
            //returns if a run is defined or not
            //NOTE: if a run is defined, it always has length 1 and therefore it always
            //      represents a single grid point

            return (s_level==0);
        }

        size_type active_pt_id() const {
            //returns the active_pt_id if the current run is an active grid point
            //returns INACTIVE if not
            if (is_defined()) {
                return active_pt_id2();
            } else {
                return INACTIVE;
            }
        }

        size_type active_pt_id2() const {               //TODO
            //the same as "active_pt_id", however this function assumes
            //that the current run is defined
            //shfdhsfhdskjhgf assert(l.sub_levelsets.size()>sub);
            //shfdhsfhdskjhgf assert(sub>=0);
            //shfdhsfhdskjhgf assert(l.sub_levelsets.size()==l.active_pt_id_offsets.size());
            //shfdhsfhdskjhgf assert(current_run_code()>=0);
            //shfdhsfhdskjhgf assert(current_run_code()<l.sub_levelsets[sub].distances.size());

            //shfdhsfhdskjhgf assert(l.sub_levelsets[sub].active.size()>current_run_code());

            size_type tmp=l.sub_levelsets[sub].active[current_run_code()];
            if (tmp!=INACTIVE) {
                tmp+=l.active_pt_id_offsets[sub];
                //shfdhsfhdskjhgf assert(tmp>=0);
                //shfdhsfhdskjhgf assert(tmp<l.num_active_pts());
            }
            return tmp;
        }

        value_type value() const {
            //returns the level set value for the current run
            //if the run is undefined either POS_VALUE or NEG_VALUE is returned
            if (is_defined()) {
                return value2();
            } else if (current_run_code()==POS_PT) {
                return POS_VALUE;
            } else {
                //assert(current_run_code()==NEG_PT);                                                  //TODO
                return NEG_VALUE;
            }
        }

        value_type value2() const {
            //the same as "value", however this function assumes
            //that the current run is defined, and therefore no check is required
            //if the run is undefined or not
            //shfdhsfhdskjhgf assert(current_run_code()>=0);
            //shfdhsfhdskjhgf assert(current_run_code()<l.sub_levelsets[sub].distances.size());
            return l.sub_levelsets[sub].distances[current_run_code()];
        }

        value_type value_p() const{
            //returns the level set value for the current run
            //if the run is undefined positive, POS_VALUE is returned, else NEG_VALUE
            if (is_defined()) {
                return value2();
            } else if (current_run_code()==POS_PT) {
                return POS_VALUE;
            } else {
                //assert(current_run_code()==NEG_PT);
                if(any(abs_coords, l.Grid.max_point_index())){
                    return -0.5;             //TODO: add border checking
                }
                else if(any(abs_coords, l.Grid.min_point_index())){
                    return -0.5;             //TODO: add border checking
                }
                return NEG_VALUE;
            }
        }

        bool is_active() const {
            //this function returns if a grid point is active or not
            return (active_pt_id()!=INACTIVE);
        }

        bool is_active2() const {
            //the same as "is_active", however this function assumes that
            //the current run is a defined run/grid point
            return (active_pt_id2()!=INACTIVE);
        }


        sign_type sign() const {
            //returns the sign of the current run
            return l.sign(value());
        }

        sign_type sign2() const {
            //the same as "sign", however it is assumed that the current run
            //is defined
            return l.sign(value2());
        }

        sign_type sign_p() const{
            return l.sign(value_p());
        }


    };



    template <class GridTraitsType, class LevelSetTraitsType>
    class levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_runs:public levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_base {
        //this iterator is the main iterator, on which all other iterators base on
        //this iterator iterates over all runs and stops at every start index of each run
        //this run can be defined or undefined, if it is a defined run it corresponds to a single defined grid point

        using levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_base::l;
        using levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_base::r_level;
        using levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_base::s_level;
        using levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_base::sub;
        using levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_base::abs_coords;
        using levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_base::end_abs_coords;
        using levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_base::start_run_abs_coords;
        using levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_base::end_run_abs_coords;
        using levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_base::run_type_pos;
        using levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_base::start_indices_pos;
        using levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_base::go_up_AB;
        using levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_base::go_up_BA;
    public:
        using levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_base::is_finished;
    private:

        template <class I>
        friend void levelset<GridTraitsType, LevelSetTraitsType>::next(I &) const;

        template <class I>
        friend void levelset<GridTraitsType, LevelSetTraitsType>::previous(I &) const;

        template <class V, class I>
        friend void levelset<GridTraitsType, LevelSetTraitsType>::go_to_indices(const V&, I&) const;

        template <class V, class I>
        friend void levelset<GridTraitsType, LevelSetTraitsType>::go_to_indices(int, const V&, I&) const;


        template <class V, class I>
        friend void levelset<GridTraitsType, LevelSetTraitsType>::go_to_indices_sequential(const V&, I&) const;



    protected:

        void go_down_AB(index_type c)  {          //find right run_type_pos

            //shfdhsfhdskjhgf assert(s_level==r_level);

            const sub_levelset_type & sl= l.sub_levelsets[sub];


            const size_type & s=start_indices_pos[s_level];

            --r_level;

            size_type r=sl.start_indices[r_level][s];

            typename std::vector<index_type>::const_iterator start_breaks = sl.runbreaks[r_level].begin()+(r-s);
            typename std::vector<index_type>::const_iterator end_breaks= sl.runbreaks[r_level].begin()+(sl.GetStartIndex(r_level, s+1)-(s+1));

            //shfdhsfhdskjhgf assert(start_breaks<=end_breaks);

            typename std::vector<index_type>::const_iterator pos_breaks=std::upper_bound(
                            start_breaks,
                            end_breaks,
                            c);

            r+=(pos_breaks-start_breaks);

            if (pos_breaks==start_breaks) {
                start_run_abs_coords[r_level]=l.grid().min_point_index(r_level);
            } else {
                //shfdhsfhdskjhgf assert(pos_breaks>start_breaks);
                start_run_abs_coords[r_level]=*(pos_breaks-1);
            }

            if (pos_breaks==end_breaks) {
                end_run_abs_coords[r_level]=l.grid().max_point_index(r_level);
            } else {
                //shfdhsfhdskjhgf assert(pos_breaks<end_breaks);
                end_run_abs_coords[r_level]=(*pos_breaks)-1;
            }

            run_type_pos[r_level]=r;

            //shfdhsfhdskjhgf assert(s_level>=1);
            //shfdhsfhdskjhgf assert(s_level==r_level+1);
        }

        void go_down_AB_first()  {          //find right run_type_pos

            //shfdhsfhdskjhgf assert(s_level==r_level);

            const sub_levelset_type & sl= l.sub_levelsets[sub];

            const size_type & s=start_indices_pos[s_level];

            --r_level;

            const size_type r=sl.start_indices[r_level][s];

            start_run_abs_coords[r_level]=l.grid().min_point_index(r_level);

            end_run_abs_coords[r_level]=sl.GetRunEndCoord(r_level, s ,r);

            run_type_pos[r_level]=r;

            //shfdhsfhdskjhgf assert(s_level>=1);
            //shfdhsfhdskjhgf assert(s_level==r_level+1);
        }

        void go_down_AB_last()  {          //find right run_type_pos

            //shfdhsfhdskjhgf assert(s_level==r_level);

            const sub_levelset_type & sl= l.sub_levelsets[sub];

            const size_type & s=start_indices_pos[s_level];

            --r_level;

            const size_type r=sl.GetStartIndex(r_level, s+1)-1;

            start_run_abs_coords[r_level]=sl.GetRunStartCoord(r_level, s ,r);

            end_run_abs_coords[r_level]=l.grid().max_point_index(r_level);

            run_type_pos[r_level]=r;

            //shfdhsfhdskjhgf assert(s_level>=1);
            //shfdhsfhdskjhgf assert(s_level==r_level+1);
        }

        bool go_down_BA(index_type c) {
            //shfdhsfhdskjhgf assert(s_level==r_level+1);
            //shfdhsfhdskjhgf assert(s_level>=1);

            const sub_levelset_type & sl= l.sub_levelsets[sub];

            start_indices_pos[s_level-1]=sl.runtypes[r_level][run_type_pos[r_level]];
            if (l.is_defined(start_indices_pos[s_level-1])) {
                --s_level;
                start_indices_pos[s_level]+=(c-start_run_abs_coords[r_level]);
                abs_coords[s_level]=c;
                end_abs_coords[s_level]=c;
                return true;
            } else {
                abs_coords[r_level]=start_run_abs_coords[r_level];
                end_abs_coords[r_level]=end_run_abs_coords[r_level];
                return false;
            }
        }

        bool go_next_A() {
            if (s_level==r_level) {
                //if (start_indices_pos[s_level]<l.runtypes[r_level][run_type_pos[r_level]]+(end_run_coords[r_level]-start_run_coords[r_level])) {
                if (abs_coords[s_level]<end_run_abs_coords[r_level]) {
                    ++start_indices_pos[s_level];
                    ++abs_coords[s_level];
                    ++end_abs_coords[s_level];
                    return true;
                }
            }
            return false;
        }


        bool go_next_B() {
            //shfdhsfhdskjhgf assert(s_level==r_level+1);
            //shfdhsfhdskjhgf assert(end_run_abs_coords[r_level]<=l.grid().max_point_index(r_level));

            const sub_levelset_type & sl= l.sub_levelsets[sub];

            if (end_run_abs_coords[r_level]!=l.grid().max_point_index(r_level)) {
                ++run_type_pos[r_level];
                start_run_abs_coords[r_level]=end_run_abs_coords[r_level]+1;
                end_run_abs_coords[r_level]=sl.GetRunEndCoord(r_level, start_indices_pos[s_level],run_type_pos[r_level]);
                return true;
            } else {
                return false;
            }
        }

        bool go_previous_A() {
            if (s_level==r_level) {
                if (abs_coords[s_level]>start_run_abs_coords[r_level]) {
                    --start_indices_pos[s_level];
                    --abs_coords[s_level];
                    --end_abs_coords[s_level];
                    return true;
                }
            }
            return false;
        }


        bool go_previous_B() {
            //shfdhsfhdskjhgf assert(s_level==r_level+1);
            //shfdhsfhdskjhgf assert(start_run_abs_coords[r_level]>=l.grid().min_point_index(r_level));

            const sub_levelset_type & sl= l.sub_levelsets[sub];

            if (start_run_abs_coords[r_level]!=l.grid().min_point_index(r_level)) {
                --run_type_pos[r_level];
                end_run_abs_coords[r_level]=start_run_abs_coords[r_level]-1;
                start_run_abs_coords[r_level]=sl.GetRunStartCoord(r_level, start_indices_pos[s_level],run_type_pos[r_level]);
                return true;
            } else {
                return false;
            }
        }


    public:

        void print() const {
            levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_base::print();
        }

        const_iterator_runs(const levelset& lx, bool reverse=false):const_iterator_base(lx) {
            if (reverse) {
                go_to_indices(l.grid().max_point_index());
            } else {
                go_to_indices(l.grid().min_point_index());
            }
        }


        template <class V>
        const_iterator_runs(const levelset& lx, const V& v):const_iterator_base(lx) {
            go_to_indices(v);
        }

        // TODO: it would make more sense if the call to next returned a bool whether it was possible to move or not. currently calling next on a finished iterator results in a segfault because r_level==D in end_run_abs_coords[r_level] in go_next_B, it would make more sense to check is_finished first, do nothing if it is and return false
        void next() {
            l.next(*this);
        }

        void previous() {
            l.previous(*this);
        }


        template <class V>
        void go_to_indices(const V& v) {
            l.go_to_indices(v, *this);
        }

        template <class V>
        void go_to_indices_sequential(const V& v) {
             l.go_to_indices_sequential(v, *this);
        }

    };

    template <class GridTraitsType, class LevelSetTraitsType>
    template <class I>
    void levelset<GridTraitsType, LevelSetTraitsType>::next(I & it) const {
         while(true) {
            if (it.go_next_A()) break;
            it.go_up_AB();
            if (it.go_next_B()) {
                it.go_down_BA(it.start_run_abs_coords[it.r_level]);
                break;
            }
            it.go_up_BA();
            //shfdhsfhdskjhgf assert(it.r_level==it.s_level);
            if (it.r_level==D) {
                it.abs_coords=it.l.grid().increment_indices(it.l.grid().max_point_index());
                it.end_abs_coords=it.abs_coords;        //TODO
                //it.start_indices_pos[D]=it.l.num_pts();
                return;
            }
        }

        while (it.r_level == it.s_level  && it.s_level>0) {
            const index_type c=it.l.grid().min_point_index(it.r_level-1);
            //go_down_AB(c);
            it.go_down_AB_first();
            it.go_down_BA(c);
        }

        int s=it.current_segment();
        if (s!=it.sub) {
            go_to_indices(s,it.start_indices(),it);
            /*if (s>0) {
                //shfdhsfhdskjhgf assert(it.l.segmentation[s-1]==it.start_indices());
            } else {
                //shfdhsfhdskjhgf assert(it.l.grid().min_point_index()==it.start_indices());
            }*/
        }
    }

    template <class GridTraitsType, class LevelSetTraitsType>
    template <class I>
    void levelset<GridTraitsType, LevelSetTraitsType>::previous(I & it) const {
         while(true) {
            if (it.go_previous_A()) break;
            it.go_up_AB();
            if (it.go_previous_B()) {
                it.go_down_BA(it.end_run_abs_coords[it.r_level]);
                break;
            }
            it.go_up_BA();
            //shfdhsfhdskjhgf assert(it.r_level==it.s_level);
            if (it.r_level==D) {
                it.abs_coords=it.l.grid().decrement_indices(it.l.grid().min_point_index());
                it.end_abs_coords=it.abs_coords;        //TODO
                return;
            }
        }

       while (it.r_level == it.s_level  && it.s_level>0) {
            const index_type c=it.l.grid().max_point_index(it.r_level-1);
            //go_down_AB(c);
            it.go_down_AB_last();
            it.go_down_BA(c);
       }

       int s=it.current_segment();
       if (s!=it.sub) {
           go_to_indices(s,it.end_indices(),it);
           /*if (s<it.l.segmentation.size()) {
               //shfdhsfhdskjhgf assert(it.l.grid().decrement_indices(it.l.segmentation[s])==it.end_indices());
           } else {
               //shfdhsfhdskjhgf assert(it.l.grid().max_point_index()==it.end_indices());
           }*/
       }
    }




    template <class GridTraitsType, class LevelSetTraitsType>
    template <class V, class I>
    void levelset<GridTraitsType, LevelSetTraitsType>::go_to_indices(int sub, const V& v, I& it) const {
        it.r_level=D;
        it.s_level=D;
        it.sub=sub;
        it.start_indices_pos[D]=0;
        do {
            //shfdhsfhdskjhgf assert(it.r_level==it.s_level);
            const index_type c=v[it.r_level-1];

            //shfdhsfhdskjhgf assert(is_defined(it.start_indices_pos[it.s_level]));

            it.go_down_AB(c);

            //shfdhsfhdskjhgf assert(c>=it.start_run_abs_coords[it.r_level]);
            //shfdhsfhdskjhgf assert(it.end_run_abs_coords[it.r_level]>=c);

            it.go_down_BA(c);
        } while (it.r_level == it.s_level  && it.s_level>0);
        //shfdhsfhdskjhgf assert(!it.is_finished());

        for (int h=0;h<it.r_level;++h) {
            it.abs_coords[h]=it.l.grid().min_point_index(h);
            it.end_abs_coords[h]=it.l.grid().max_point_index(h);
        }

        //TODO initialize abs_coords
    }

    template <class GridTraitsType, class LevelSetTraitsType>
    template <class V, class I>
    void levelset<GridTraitsType, LevelSetTraitsType>::go_to_indices(const V& v, I& it) const {
        go_to_indices(0, v, it);        //TODO
        int s=it.current_segment();
        if (s!=0) go_to_indices(s,v,it);
    }




    template <class GridTraitsType, class LevelSetTraitsType>
    template <class V, class I>
    void levelset<GridTraitsType, LevelSetTraitsType>::go_to_indices_sequential(const V& v, I& it) const {

        //shfdhsfhdskjhgf assert(it.r_level<D);
        //shfdhsfhdskjhgf assert(it.r_level<=it.s_level);

        int lvl=D-1;

        while(true) {

            //shfdhsfhdskjhgf assert(lvl>=it.r_level);

            if (it.end_run_abs_coords[lvl]<v[lvl]) {
                it.r_level=lvl;
                it.s_level=lvl+1;
                do {
                    if (!it.go_next_B()) {
                        //shfdhsfhdskjhgf assert(0);
                    }
                } while (it.end_run_abs_coords[it.r_level]<v[it.r_level]);
                //shfdhsfhdskjhgf assert(it.s_level==it.r_level+1);
                it.go_down_BA(v[lvl]);

            } else if (it.start_run_abs_coords[lvl]>v[lvl]) {
                it.r_level=lvl;
                it.s_level=lvl+1;
                do {
                    if (!it.go_previous_B()) {
                        //shfdhsfhdskjhgf assert(0);
                    }
                } while (it.start_run_abs_coords[it.r_level]>v[it.r_level]);
                //shfdhsfhdskjhgf assert(it.s_level==it.r_level+1);
                it.go_down_BA(v[lvl]);

            } else if (lvl>=it.s_level) {
                if (it.abs_coords[lvl]!=v[lvl]) {
                    it.r_level=lvl;
                    it.s_level=lvl+1;
                    //shfdhsfhdskjhgf assert(it.s_level==it.r_level+1);
                    it.go_down_BA(v[lvl]);

                } else if (lvl>0) {
                    --lvl;
                    continue;
                }
            }
            break;
        }

        while (it.r_level == it.s_level  && it.s_level>0) {
            //shfdhsfhdskjhgf assert(it.r_level==it.s_level);
            const index_type c=v[it.r_level-1];

            //shfdhsfhdskjhgf assert(is_defined(it.start_indices_pos[it.s_level]));

            it.go_down_AB(c);

            //shfdhsfhdskjhgf assert(c>=it.start_run_abs_coords[it.r_level]);
            //shfdhsfhdskjhgf assert(it.end_run_abs_coords[it.r_level]>=c);

            it.go_down_BA(c);
        }
        //shfdhsfhdskjhgf assert(!it.is_finished());

        for (int h=0;h<it.r_level;++h) {
            it.abs_coords[h]=it.l.grid().min_point_index(h);
            it.end_abs_coords[h]=it.l.grid().max_point_index(h);
        }

        int s=it.current_segment();
        if (s!=it.sub) go_to_indices(s,v,it);
    }


    template <class GridTraitsType, class LevelSetTraitsType>
    class levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_runs_offset:public levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_base {
        //this iterator is the main iterator, on which all other iterators base on
        //this iterator iterates over all runs and stops at every start index of each run
        //this run can be defined or undefined, if it is a defined run it corresponds to a single defined grid point

        using levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_base::l;
        using levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_base::r_level;
        using levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_base::s_level;
        using levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_base::sub;
        using levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_base::abs_coords;
        using levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_base::end_abs_coords;
        using levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_base::start_run_abs_coords;
        using levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_base::end_run_abs_coords;
        using levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_base::run_type_pos;
        using levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_base::start_indices_pos;
        using levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_base::go_up_AB;
        using levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_base::go_up_BA;
    public:
        using levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_base::is_finished;
    private:

        template <class I>
        friend void levelset<GridTraitsType, LevelSetTraitsType>::next(I &) const;

        template <class I>
        friend void levelset<GridTraitsType, LevelSetTraitsType>::previous(I &) const;

        template <class V, class I>
        friend void levelset<GridTraitsType, LevelSetTraitsType>::go_to_indices(const V&, I&) const;

        template <class V, class I>
        friend void levelset<GridTraitsType, LevelSetTraitsType>::go_to_indices(int, const V&, I&) const;

        template <class V, class I>
        friend void levelset<GridTraitsType, LevelSetTraitsType>::go_to_indices_sequential(const V&, I&) const;

        vec<index_type,D> start_run_rel_coords;
        vec<index_type,D> end_run_rel_coords;
        vec<index_type,D> rel_coords;
        vec<index_type,D> offset;
        std::bitset<D> move_inverse;

    protected:

        void go_down_AB(index_type abs_c)  { //X         //find right run_type_pos

            //shfdhsfhdskjhgf assert(s_level==r_level);

            //shfdhsfhdskjhgf assert(l.is_defined(start_indices_pos[s_level]));

            const sub_levelset_type & sl = l.sub_levelsets[sub];

            const size_type & s=start_indices_pos[s_level];

            --r_level;

            int cycles=0;
            index_type rel_c=l.grid().global_index_2_local_index(r_level,abs_c,offset[r_level],cycles);

            //shfdhsfhdskjhgf assert(s<sl.start_indices[r_level].size());
            //shfdhsfhdskjhgf assert(r_level>=0);

            size_type r=sl.start_indices[r_level][s];

            typename std::vector<index_type>::const_iterator start_breaks = sl.runbreaks[r_level].begin()+(r-s);
            typename std::vector<index_type>::const_iterator end_breaks= sl.runbreaks[r_level].begin()+(sl.GetStartIndex(r_level, s+1)-(s+1));

            //shfdhsfhdskjhgf assert(start_breaks<=end_breaks);

            typename std::vector<index_type>::const_iterator pos_breaks=std::upper_bound(
                            start_breaks,
                            end_breaks,
                            rel_c);

            r+=(pos_breaks-start_breaks);

            run_type_pos[r_level]=r;

            if (pos_breaks==start_breaks) {
                start_run_rel_coords[r_level]=l.grid().min_point_index(r_level);
            } else {
                //shfdhsfhdskjhgf assert(pos_breaks>start_breaks);
                start_run_rel_coords[r_level]=*(pos_breaks-1);
            }

            if (pos_breaks==end_breaks) {
                end_run_rel_coords[r_level]=l.grid().max_point_index(r_level);
            } else {
                //shfdhsfhdskjhgf assert(pos_breaks<end_breaks);
                end_run_rel_coords[r_level]=(*pos_breaks)-1;
            }

            //shfdhsfhdskjhgf assert(start_run_rel_coords[r_level]>=l.grid().min_point_index(r_level));
            //shfdhsfhdskjhgf assert(start_run_rel_coords[r_level]<=rel_c);
            //shfdhsfhdskjhgf assert(rel_c<=end_run_rel_coords[r_level]);
            //shfdhsfhdskjhgf assert(l.grid().max_point_index(r_level)>=end_run_rel_coords[r_level]);

            //calculate run_abs_coords


            if (l.grid().is_boundary_periodic(r_level)) {

                const index_type & rel_s = start_run_rel_coords[r_level];
                const index_type & rel_e =   end_run_rel_coords[r_level];


                move_inverse.reset(r_level);

                start_run_abs_coords[r_level]=std::max(
                                                        l.grid().local_index_2_global_index(r_level,rel_s,cycles, offset[r_level]),
                                                        l.grid().min_point_index(r_level)
                                                    );
                end_run_abs_coords[r_level]=std::min(
                                                            l.grid().local_index_2_global_index(r_level,rel_e,cycles, offset[r_level]),
                                                            l.grid().max_point_index(r_level)
                                                        );

            } else {


                move_inverse.set(r_level, cycles & 1);

                if (start_breaks==end_breaks) {

                    start_run_abs_coords[r_level]=l.grid().min_point_index(r_level);
                    end_run_abs_coords[r_level]=l.grid().max_point_index(r_level);

                } else {

                    index_type rel_s = start_run_rel_coords[r_level];
                    index_type rel_e =   end_run_rel_coords[r_level];


                    if (pos_breaks==start_breaks)
                        rel_s+=(start_run_rel_coords[r_level]-end_run_rel_coords[r_level]);
                    else if (pos_breaks==end_breaks)
                        rel_e+=(end_run_rel_coords[r_level]-start_run_rel_coords[r_level]);

                    if ((cycles & 1) !=0) std::swap(rel_s,rel_e);

                    start_run_abs_coords[r_level]=std::max(
                                                            l.grid().local_index_2_global_index(r_level,rel_s,cycles, offset[r_level]),
                                                            l.grid().min_point_index(r_level)
                                                        );

                    end_run_abs_coords[r_level]=std::min(
                                                                l.grid().local_index_2_global_index(r_level,rel_e,cycles, offset[r_level]),
                                                                l.grid().max_point_index(r_level)
                                                            );
                }
            }



            //shfdhsfhdskjhgf assert(start_run_abs_coords[r_level]>=l.grid().min_point_index(r_level));
            //shfdhsfhdskjhgf assert(start_run_abs_coords[r_level]<=abs_c);
            //shfdhsfhdskjhgf assert(abs_c<=end_run_abs_coords[r_level]);
            //shfdhsfhdskjhgf assert(l.grid().max_point_index(r_level)>=end_run_abs_coords[r_level]);

            //shfdhsfhdskjhgf assert(s_level>=1);
            //shfdhsfhdskjhgf assert(s_level==r_level+1);

            if (end_run_abs_coords[r_level]-start_run_abs_coords[r_level]> end_run_rel_coords[r_level]-start_run_rel_coords[r_level]) {
                //shfdhsfhdskjhgf assert(!(l.grid().is_boundary_periodic(r_level)));
                //shfdhsfhdskjhgf assert(end_run_rel_coords[r_level]==l.grid().max_point_index(r_level) || start_run_rel_coords[r_level]==l.grid().min_point_index(r_level));
            }

             //shfdhsfhdskjhgf assert(l.grid().global_index_2_local_index(r_level,start_run_abs_coords[r_level], offset[r_level])>=start_run_rel_coords[r_level]);
             //shfdhsfhdskjhgf assert(l.grid().global_index_2_local_index(r_level,start_run_abs_coords[r_level], offset[r_level])<=end_run_rel_coords[r_level]);

             //shfdhsfhdskjhgf assert(l.grid().global_index_2_local_index(r_level,end_run_abs_coords[r_level], offset[r_level])>=start_run_rel_coords[r_level]);
             //shfdhsfhdskjhgf assert(l.grid().global_index_2_local_index(r_level,end_run_abs_coords[r_level], offset[r_level])<=end_run_rel_coords[r_level]);



        }


        void go_down_AB_first()  {          //find right run_type_pos

            //shfdhsfhdskjhgf assert(s_level==r_level);

            const index_type c=l.grid().min_point_index(r_level-1);
            //shfdhsfhdskjhgf assert(l.is_defined(start_indices_pos[s_level]));
            go_down_AB(c);

            //shfdhsfhdskjhgf assert(s_level>=1);
            //shfdhsfhdskjhgf assert(s_level==r_level+1);
        }

        void go_down_AB_last()  {          //find right run_type_pos

            //shfdhsfhdskjhgf assert(s_level==r_level);

            const index_type c=l.grid().max_point_index(r_level-1);
            //shfdhsfhdskjhgf assert(l.is_defined(start_indices_pos[s_level]));
            go_down_AB(c);

            //shfdhsfhdskjhgf assert(s_level>=1);
            //shfdhsfhdskjhgf assert(s_level==r_level+1);
        }


        bool go_down_BA(index_type abs_c) {           //X


            //shfdhsfhdskjhgf assert(s_level==r_level+1);
            //shfdhsfhdskjhgf assert(s_level>=1);

            const sub_levelset_type & sl = l.sub_levelsets[sub];

            start_indices_pos[s_level-1]=sl.runtypes[r_level][run_type_pos[r_level]];
            if (l.is_defined(start_indices_pos[s_level-1])) {
                --s_level;
                //shfdhsfhdskjhgf assert(l.is_defined(start_indices_pos[s_level]));

                index_type rel_c=l.grid().global_index_2_local_index(r_level,abs_c,offset[r_level]);

                //shfdhsfhdskjhgf assert(rel_c>=start_run_rel_coords[r_level]);
                //shfdhsfhdskjhgf assert(rel_c<=end_run_rel_coords[r_level]);

                start_indices_pos[s_level]+=(rel_c-start_run_rel_coords[r_level]);
                rel_coords[s_level]=rel_c;
                abs_coords[s_level]=abs_c;
                end_abs_coords[s_level]=abs_c;


                //shfdhsfhdskjhgf assert(l.is_defined(start_indices_pos[s_level]));

                return true;
            } else {
                abs_coords[r_level]=start_run_abs_coords[r_level];
                end_abs_coords[r_level]=end_run_abs_coords[r_level];

                //shfdhsfhdskjhgf assert(s_level==r_level+1);
                return false;
            }
        }

        bool go_next_A() {                  //X

            if (s_level==r_level) {

                //shfdhsfhdskjhgf assert(rel_coords[s_level]<=end_run_rel_coords[r_level]);
                //shfdhsfhdskjhgf assert(rel_coords[s_level]>=start_run_rel_coords[r_level]);
                //shfdhsfhdskjhgf assert(abs_coords[s_level]<=end_run_abs_coords[r_level]);
                //shfdhsfhdskjhgf assert(abs_coords[s_level]>=start_run_abs_coords[r_level]);

                //if (start_indices_pos[s_level]<l.runtypes[r_level][run_type_pos[r_level]]+(end_run_coords[r_level]-start_run_coords[r_level])) {

                if (abs_coords[s_level]<end_run_abs_coords[r_level]) {
                    ++abs_coords[s_level];
                    ++end_abs_coords[s_level];

                    if (l.grid().is_boundary_periodic(r_level)) {
                        //shfdhsfhdskjhgf assert(rel_coords[s_level]!=l.grid().max_point_index(s_level));
                        ++rel_coords[s_level];
                        ++start_indices_pos[s_level];
                    } else {

                      if (rel_coords[s_level]==l.grid().max_point_index(s_level)) {
              move_inverse.set(s_level);
            } else if (rel_coords[s_level]==l.grid().min_point_index(s_level)) {
              move_inverse.reset(s_level);
            }

                      if (move_inverse.test(s_level)) {
                            //shfdhsfhdskjhgf assert(start_indices_pos[s_level]>0);
                            --start_indices_pos[s_level];
                            --rel_coords[s_level];
                        } else {
                            ++start_indices_pos[s_level];
                            ++rel_coords[s_level];
                        }

                    }

                    //shfdhsfhdskjhgf assert(rel_coords[s_level]<=end_run_rel_coords[r_level]);
                    //shfdhsfhdskjhgf assert(rel_coords[s_level]>=start_run_rel_coords[r_level]);
                    //shfdhsfhdskjhgf assert(abs_coords[s_level]<=end_run_abs_coords[r_level]);
                    //shfdhsfhdskjhgf assert(abs_coords[s_level]>=start_run_abs_coords[r_level]);

                    return true;
                }
            }
            return false;
        }

        bool go_previous_A() {

            if (s_level==r_level) {

                //shfdhsfhdskjhgf assert(rel_coords[s_level]<=end_run_rel_coords[r_level]);
                //shfdhsfhdskjhgf assert(rel_coords[s_level]>=start_run_rel_coords[r_level]);
                //shfdhsfhdskjhgf assert(abs_coords[s_level]<=end_run_abs_coords[r_level]);
                //shfdhsfhdskjhgf assert(abs_coords[s_level]>=start_run_abs_coords[r_level]);

                if (abs_coords[s_level]>start_run_abs_coords[r_level]) {
                    --abs_coords[s_level];
                    --end_abs_coords[s_level];

                    if (l.grid().is_boundary_periodic(r_level)) {
                        //shfdhsfhdskjhgf assert(rel_coords[s_level]!=l.grid().min_point_index(s_level));
                        --rel_coords[s_level];
                        --start_indices_pos[s_level];
                    } else {
                      if (rel_coords[s_level]==l.grid().max_point_index(s_level)) {
                        move_inverse.reset(s_level);
                      } else if (rel_coords[s_level]==l.grid().min_point_index(s_level)) {
                        move_inverse.set(s_level);
                      }
                        if (move_inverse.test(s_level)) {
              ++start_indices_pos[s_level];
              ++rel_coords[s_level];
                        } else {
              --start_indices_pos[s_level];
              --rel_coords[s_level];
                        }
                    }

                    //shfdhsfhdskjhgf assert(rel_coords[s_level]<=end_run_rel_coords[r_level]);
                    //shfdhsfhdskjhgf assert(rel_coords[s_level]>=start_run_rel_coords[r_level]);
                    //shfdhsfhdskjhgf assert(abs_coords[s_level]<=end_run_abs_coords[r_level]);
                    //shfdhsfhdskjhgf assert(abs_coords[s_level]>=start_run_abs_coords[r_level]);

                    return true;
                }
            }
            return false;
        }

        bool go_next_B() {

            //shfdhsfhdskjhgf assert(s_level==r_level+1);
            //shfdhsfhdskjhgf assert(end_run_abs_coords[r_level]<=l.grid().max_point_index(r_level));
            //shfdhsfhdskjhgf assert(start_run_abs_coords[r_level]>=l.grid().min_point_index(r_level));

            const sub_levelset_type & sl = l.sub_levelsets[sub];

            //shfdhsfhdskjhgf assert(run_type_pos[r_level]>=sl.GetStartIndex(r_level, start_indices_pos[s_level]));
            //shfdhsfhdskjhgf assert(run_type_pos[r_level]<sl.GetStartIndex(r_level, start_indices_pos[s_level]+1));

            if (end_run_abs_coords[r_level]!=l.grid().max_point_index(r_level)) {

                if (l.grid().is_boundary_periodic(r_level)) {       //periodic boundary conditions

                    if (end_run_rel_coords[r_level]==l.grid().max_point_index(r_level)) {
                        run_type_pos[r_level]=sl.GetStartIndex(r_level, start_indices_pos[s_level]);
                        start_run_rel_coords[r_level]=l.grid().min_point_index(r_level);
                    } else {
                         ++run_type_pos[r_level];
                        start_run_rel_coords[r_level]=end_run_rel_coords[r_level]+1;
                    }

                    end_run_rel_coords[r_level]=sl.GetRunEndCoord(r_level, start_indices_pos[s_level],run_type_pos[r_level]);

                    ++end_run_abs_coords[r_level];
                    start_run_abs_coords[r_level]=end_run_abs_coords[r_level];


                } else {

                    if (start_run_rel_coords[r_level]==l.grid().min_point_index(r_level)) {
                        move_inverse.reset(r_level);
                    } else if (end_run_rel_coords[r_level]==l.grid().max_point_index(r_level)) {
                        move_inverse.set(r_level);
                    }

                    if (move_inverse.test(r_level)) {

                         //shfdhsfhdskjhgf assert(run_type_pos[r_level]>sl.GetStartIndex(r_level, start_indices_pos[s_level]));

                        --run_type_pos[r_level];

                        end_run_rel_coords[r_level]=start_run_rel_coords[r_level]-1;
                        start_run_rel_coords[r_level]=sl.GetRunStartCoord(r_level, start_indices_pos[s_level],run_type_pos[r_level]);

                        ++end_run_abs_coords[r_level];
                        start_run_abs_coords[r_level]=end_run_abs_coords[r_level];

                        if (start_run_rel_coords[r_level]==l.grid().min_point_index(r_level)) {
                            end_run_abs_coords[r_level]+=(end_run_rel_coords[r_level]-start_run_rel_coords[r_level]);
                        }

                    } else {

                        //shfdhsfhdskjhgf assert(run_type_pos[r_level]<sl.GetStartIndex(r_level, start_indices_pos[s_level]+1)-1);

                        ++run_type_pos[r_level];

                        start_run_rel_coords[r_level]=end_run_rel_coords[r_level]+1;
                        end_run_rel_coords[r_level]=sl.GetRunEndCoord(r_level, start_indices_pos[s_level],run_type_pos[r_level]);

                        ++end_run_abs_coords[r_level];
                        start_run_abs_coords[r_level]=end_run_abs_coords[r_level];

                        if (end_run_rel_coords[r_level]==l.grid().max_point_index(r_level)) {
                            end_run_abs_coords[r_level]+=(end_run_rel_coords[r_level]-start_run_rel_coords[r_level]);
                        }
                    }

                }

                end_run_abs_coords[r_level]+=(end_run_rel_coords[r_level]-start_run_rel_coords[r_level]);
                if (end_run_abs_coords[r_level]>l.grid().max_point_index(r_level)) end_run_abs_coords[r_level]=l.grid().max_point_index(r_level);


                //shfdhsfhdskjhgf assert(run_type_pos[r_level]>=sl.GetStartIndex(r_level, start_indices_pos[s_level]));
                //shfdhsfhdskjhgf assert(run_type_pos[r_level]<sl.GetStartIndex(r_level, start_indices_pos[s_level]+1));


                //shfdhsfhdskjhgf assert(start_run_rel_coords[r_level]>=l.grid().min_point_index(r_level));
                //shfdhsfhdskjhgf assert(start_run_rel_coords[r_level]<=end_run_rel_coords[r_level]);
                //shfdhsfhdskjhgf assert(l.grid().max_point_index(r_level)>=end_run_rel_coords[r_level]);

                //shfdhsfhdskjhgf assert(start_run_abs_coords[r_level]>=l.grid().min_point_index(r_level));
                //shfdhsfhdskjhgf assert(start_run_abs_coords[r_level]<=end_run_abs_coords[r_level]);
                //shfdhsfhdskjhgf assert(l.grid().max_point_index(r_level)>=end_run_abs_coords[r_level]);


                if (end_run_abs_coords[r_level]-start_run_abs_coords[r_level]> end_run_rel_coords[r_level]-start_run_rel_coords[r_level]) {
                    //shfdhsfhdskjhgf assert(!(l.grid().is_boundary_periodic(r_level)));
                    //shfdhsfhdskjhgf assert(end_run_rel_coords[r_level]==l.grid().max_point_index(r_level) || start_run_rel_coords[r_level]==l.grid().min_point_index(r_level));
                }

                 //shfdhsfhdskjhgf assert(l.grid().global_index_2_local_index(r_level,start_run_abs_coords[r_level], offset[r_level])>=start_run_rel_coords[r_level]);
                 //shfdhsfhdskjhgf assert(l.grid().global_index_2_local_index(r_level,start_run_abs_coords[r_level], offset[r_level])<=end_run_rel_coords[r_level]);

                 //shfdhsfhdskjhgf assert(l.grid().global_index_2_local_index(r_level,end_run_abs_coords[r_level], offset[r_level])>=start_run_rel_coords[r_level]);
                 //shfdhsfhdskjhgf assert(l.grid().global_index_2_local_index(r_level,end_run_abs_coords[r_level], offset[r_level])<=end_run_rel_coords[r_level]);

                return true;
            } else {
                return false;
            }

        }


        bool go_previous_B() {

            //shfdhsfhdskjhgf assert(s_level==r_level+1);
            //shfdhsfhdskjhgf assert(end_run_abs_coords[r_level]<=l.grid().max_point_index(r_level));
            //shfdhsfhdskjhgf assert(start_run_abs_coords[r_level]>=l.grid().min_point_index(r_level));

            const sub_levelset_type & sl = l.sub_levelsets[sub];

            //shfdhsfhdskjhgf assert(run_type_pos[r_level]>=sl.GetStartIndex(r_level, start_indices_pos[s_level]));
            //shfdhsfhdskjhgf assert(run_type_pos[r_level]<sl.GetStartIndex(r_level, start_indices_pos[s_level]+1));

            if (start_run_abs_coords[r_level]!=l.grid().min_point_index(r_level)) {

                if (l.grid().is_boundary_periodic(r_level)) {       //periodic boundary conditions

                     if (start_run_rel_coords[r_level]==l.grid().min_point_index(r_level)) {
                        run_type_pos[r_level]=sl.GetStartIndex(r_level, start_indices_pos[s_level]+1)-1;
                        end_run_rel_coords[r_level]=l.grid().max_point_index(r_level);
                    } else {
                        --run_type_pos[r_level];
                        end_run_rel_coords[r_level]=start_run_rel_coords[r_level]-1;
                    }

                    start_run_rel_coords[r_level]=sl.GetRunStartCoord(r_level, start_indices_pos[s_level],run_type_pos[r_level]);

                    --start_run_abs_coords[r_level];
                    end_run_abs_coords[r_level]=start_run_abs_coords[r_level];

                } else {

                    if (start_run_rel_coords[r_level]==l.grid().min_point_index(r_level)) {
                        move_inverse.set(r_level);
                    } else if (end_run_rel_coords[r_level]==l.grid().max_point_index(r_level)) {
                        move_inverse.reset(r_level);
                    }

                    if (move_inverse.test(r_level)) {

                        //shfdhsfhdskjhgf assert(run_type_pos[r_level]<sl.GetStartIndex(r_level, start_indices_pos[s_level]+1)-1);

                        ++run_type_pos[r_level];

                        start_run_rel_coords[r_level]=end_run_rel_coords[r_level]+1;
                        end_run_rel_coords[r_level]=sl.GetRunEndCoord(r_level, start_indices_pos[s_level],run_type_pos[r_level]);

                        --start_run_abs_coords[r_level];
                        end_run_abs_coords[r_level]=start_run_abs_coords[r_level];

                        if (end_run_rel_coords[r_level]==l.grid().max_point_index(r_level)) {
                            start_run_abs_coords[r_level]-=(end_run_rel_coords[r_level]-start_run_rel_coords[r_level]);
                        }

                    } else {

                        //shfdhsfhdskjhgf assert(run_type_pos[r_level]>sl.GetStartIndex(r_level, start_indices_pos[s_level]));

                        --run_type_pos[r_level];

                        end_run_rel_coords[r_level]=start_run_rel_coords[r_level]-1;
                        start_run_rel_coords[r_level]=sl.GetRunStartCoord(r_level, start_indices_pos[s_level],run_type_pos[r_level]);

                        --start_run_abs_coords[r_level];
                        end_run_abs_coords[r_level]=start_run_abs_coords[r_level];

                        if (start_run_rel_coords[r_level]==l.grid().min_point_index(r_level)) {
                            start_run_abs_coords[r_level]-=(end_run_rel_coords[r_level]-start_run_rel_coords[r_level]);
                        }
                    }
                }

                start_run_abs_coords[r_level]-=(end_run_rel_coords[r_level]-start_run_rel_coords[r_level]);
                if (start_run_abs_coords[r_level]<l.grid().min_point_index(r_level)) start_run_abs_coords[r_level]=l.grid().min_point_index(r_level);


                //shfdhsfhdskjhgf assert(run_type_pos[r_level]>=sl.GetStartIndex(r_level, start_indices_pos[s_level]));
                //shfdhsfhdskjhgf assert(run_type_pos[r_level]<sl.GetStartIndex(r_level, start_indices_pos[s_level]+1));


                //shfdhsfhdskjhgf assert(start_run_rel_coords[r_level]>=l.grid().min_point_index(r_level));
                //shfdhsfhdskjhgf assert(start_run_rel_coords[r_level]<=end_run_rel_coords[r_level]);
                //shfdhsfhdskjhgf assert(l.grid().max_point_index(r_level)>=end_run_rel_coords[r_level]);

                //shfdhsfhdskjhgf assert(start_run_abs_coords[r_level]>=l.grid().min_point_index(r_level));
                //shfdhsfhdskjhgf assert(start_run_abs_coords[r_level]<=end_run_abs_coords[r_level]);
                //shfdhsfhdskjhgf assert(l.grid().max_point_index(r_level)>=end_run_abs_coords[r_level]);


                if (end_run_abs_coords[r_level]-start_run_abs_coords[r_level]> end_run_rel_coords[r_level]-start_run_rel_coords[r_level]) {
                    //shfdhsfhdskjhgf assert(!(l.grid().is_boundary_periodic(r_level)));
                    //shfdhsfhdskjhgf assert(end_run_rel_coords[r_level]==l.grid().max_point_index(r_level) || start_run_rel_coords[r_level]==l.grid().min_point_index(r_level));
                }

                 //shfdhsfhdskjhgf assert(l.grid().global_index_2_local_index(r_level,start_run_abs_coords[r_level], offset[r_level])>=start_run_rel_coords[r_level]);
                 //shfdhsfhdskjhgf assert(l.grid().global_index_2_local_index(r_level,start_run_abs_coords[r_level], offset[r_level])<=end_run_rel_coords[r_level]);

                 //shfdhsfhdskjhgf assert(l.grid().global_index_2_local_index(r_level,end_run_abs_coords[r_level], offset[r_level])>=start_run_rel_coords[r_level]);
                 //shfdhsfhdskjhgf assert(l.grid().global_index_2_local_index(r_level,end_run_abs_coords[r_level], offset[r_level])<=end_run_rel_coords[r_level]);

                return true;
            } else {
                return false;
            }

        }


    public:

        void print() const {
            levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_base::print();
            std::cout << "start_run_rel_coords: " << start_run_rel_coords << std::endl;
            std::cout << "end_run_rel_coords: " << end_run_rel_coords << std::endl;
            std::cout << "rel_coords: " << rel_coords << std::endl;
            std::cout << "offset: " << offset << std::endl;
            std::cout << "move_inverse: " << move_inverse << std::endl;
        }

        template <class V>
        const_iterator_runs_offset(const levelset& lx, const V & o, bool reverse=false):const_iterator_base(lx), offset(o) {
            if (reverse) {
                go_to_indices(l.grid().max_point_index());
            } else {
                go_to_indices(l.grid().min_point_index());
            }
        }

        template <class V1, class V2>
        const_iterator_runs_offset(const levelset& lx, const V1 & o, const V2& v):const_iterator_base(lx), offset(o) {
            go_to_indices(v);
        }

        void next() {
            l.next(*this);
        }

        void previous() {
            l.previous(*this);
        }

        template <class V>
        void go_to_indices(const V& v) {
            l.go_to_indices(v, *this);
        }

        template <class V>
        void go_to_indices_sequential(const V& v) {
             l.go_to_indices_sequential(v, *this);
        }

    };



    template <class GridTraitsType, class LevelSetTraitsType>
    void levelset<GridTraitsType, LevelSetTraitsType>::invert() {
        //this function inverts a level set function
        //by switching the sign of all level set values
        //and by interchanging for all undefined runs
        //POS_PT and NEG_PT

        #pragma omp parallel num_threads(sub_levelsets.size()) //use num_threads(sub_levelsets.size()) threads
        {
            {
                int p=0;
                #ifdef _OPENMP
                p=omp_get_thread_num();
                #endif

                sub_levelset_type & s=sub_levelsets[p];

                s.invert();
            }
        }
    }



    template <class GridTraitsType, class LevelSetTraitsType> class levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_neighbor {

         //this neighbor iterator consists of 2*Dimensions const_iterator_runs_offset for the neighbors and a
        //const_iterator_runs-iterator for the center
        //whenever one of these (2*Dimensions+1) iterators reach a defined grid point
        //the const_iterator_neighbor iterator stops
        //this iterator is internally used for rebuilding and expanding the level set function

        const levelset<GridTraitsType, LevelSetTraitsType> &l;
        vec<index_type,D> start_coords;
        const_iterator_runs it_mid;
        std::vector<const_iterator_runs_offset>  it_neighbors;


        /*bool is_defined() const {
            if (it_mid.is_defined()) return true;
            for (int i=0;i<2*D;i++) {
                if (it_neighbors[i].is_defined()) return true;
            }
            return false;
        }*/


    public:

        template <class V>
        const_iterator_neighbor(const levelset<GridTraitsType, LevelSetTraitsType>& lx,  const V& v) :
                                            l(lx),
                                            start_coords(v),
                                            it_mid(lx,v) {

            for (int i=0;i<2*D;i++) {
                vec<index_type,D> tv(index_type(0));
                if (i<D) tv[i]++; else tv[i-D]--;
                it_neighbors.push_back(const_iterator_runs_offset(lx, tv,v));
            }
        }

         const_iterator_neighbor(const levelset<GridTraitsType, LevelSetTraitsType>& lx) :
                                            l(lx),
                                            start_coords(lx.grid().min_point_index()),
                                            it_mid(lx) {

            for (int i=0;i<2*D;i++) {
                vec<index_type,D> tv(index_type(0));
                if (i<D) tv[i]++; else tv[i-D]--;
                it_neighbors.push_back(const_iterator_runs_offset(lx, tv));
            }
        }

        void next() {

            std::bitset<2*D+1> increment;
            increment.set(2*D);

            vec<index_type,D> end_coords=it_mid.end_indices();
            for (int i=0;i<2*D;i++) {
                switch (compare(end_coords, it_neighbors[i].end_indices())) {
                    case 1:
                        end_coords=it_neighbors[i].end_indices();
                        increment.reset();
                    case 0:
                        increment.set(i);
                }
            }

            if (increment.test(2*D)) it_mid.next();
            for (int i=0;i<2*D;i++) if (increment.test(i)) it_neighbors[i].next();


            start_coords=l.grid().increment_indices(end_coords);

        }

        const const_iterator_runs_offset& neighbor(int direction) const {
            return it_neighbors[direction];
        }

        const const_iterator_runs& center() const {
            return it_mid;
        }

        bool is_finished() const {
            return center().is_finished();
        }

        const vec<index_type,D>&  start_indices() const {
            return start_coords;
        }

        index_type  start_indices(int dir) const {
            return start_coords[dir];
        }

    };

    template <class GridTraitsType, class LevelSetTraitsType> template<class F> class levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_cells_filtered {

        //this iterator consist of 2^Dimensions const_iterator_runs_offset-iterators
        //one for each corner of a grid cell
        //all these 2^Dimensions iterators are synchronously moved
        //whenever the filter is fulfilled for one of
        //these iterators, the const_iterator_cells_filtered iterator stops

        //this iterator is used for surface extraction using the marching cubes algorithm

        const levelset<GridTraitsType, LevelSetTraitsType>& l;

        vec<index_type,D> coords;
        std::vector<const_iterator_runs_offset> it_corners;
        const F filter;

        bool is_defined() const {

            for (int i=0;i<D;i++) {
                if (
                    (!l.grid().is_boundary_periodic(i)) &&
                    (coords[i]==l.grid().max_point_index(i))
                ) return false;
            }

            for (int i=0;i<(1<<D);i++) if (filter(it_corners[i])) return true;

            return false;
        }

    public:


        void next() {

            do {

                std::bitset< (1<<D) > increment(1);

                coords=it_corners[0].end_indices();
                for (int i=1;i<(1<<D);++i) {
                    switch (compare(coords, it_corners[i].end_indices())) {
                        case 1:
                            coords=it_corners[i].end_indices();
                            increment.reset();
                        case 0:
                            increment.set(i);
                    }
                }

                for (int i=0;i<(1<<D);i++) if (increment.test(i)) it_corners[i].next();

                coords=l.grid().increment_indices(coords);

            } while ((!is_defined()) && (!is_finished()));

        }



        /*void next() {

            do {

                vec<index_type,D> tmp_coords=it_corners[0].end_indices();
                for (int i=1;i<(1<<D);i++) {
                    if (it_corners[i].end_indices()<tmp_coords) {
                        tmp_coords=it_corners[i].end_indices();
                    }
                }

                //go to next position
                for (int i=0;i<(1<<D);i++) {
                    if (tmp_coords==it_corners[i].end_indices()) {
                        it_corners[i].next();
                        coords=it_corners[i].start_indices();
                    }
                }


            } while ((!is_defined()) && (!is_finished()));

        }*/

        const_iterator_cells_filtered(const levelset<GridTraitsType, LevelSetTraitsType>& lx,  const F& f=F()) : l(lx), coords(l.grid().min_point_index()), filter(f) {

            //set all iterators = Min
            for (int i=0;i<(1<<D);i++) {
                it_corners.push_back(const_iterator_runs_offset(lx, get_corner<D, index_type>(i)));
            }

            if (!is_defined()) next();

        }

        const const_iterator_runs_offset& corner(int corner_id) const {
            //returns the corner given by the corner_id
            //0<=corner_id<2^Dimension
            //the corners are numbered in lexicographical order
            return it_corners[corner_id];
        }

        bool is_finished() const {
            //returns true if the end is reached
            return it_corners[0].is_finished();
        }

        const vec<index_type,D>&  indices() const {
            return coords;
        }

        index_type indices(int dir) const {
            return coords[dir];
        }

    };

}


#endif /*KERNEL_HPP_*/
