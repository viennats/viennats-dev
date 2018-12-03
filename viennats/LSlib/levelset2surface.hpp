#ifndef LEVELSET2SURFACE_HPP_
#define LEVELSET2SURFACE_HPP_

/* =========================================================================
   Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.

                            -----------------
                 ViennaTS - The Vienna Topography Simulator
                            -----------------

   Contact:         viennats@iue.tuwien.ac.at

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


#include <map>
#include <deque>
#include <cassert>
#include "marching_cubes.hpp"
#include "vector.hpp"
#include "math.hpp"

///Includes all level set operations and surface conersion tools.
namespace lvlset {

    namespace {

        class DefaultPointListType {
        public:
            template<class V>
            void push_back(const V& v) const {}
        };

    }


    template<class LevelSetType, class SurfaceInserterType>
    void extract(   const LevelSetType& l,                          //the level set function
                    SurfaceInserterType& srf,                       //the surface inserter type
                    typename LevelSetType::value_type eps=0.        //eps can be used to avoid degenerated triangles
                                                                    //which can occur for the marching cubes algorithm
                                                                    //the gridedge-surface intersection points are always set
                                                                    //at least a distance eps
                                                                    //away from both grid points of that grid edge
                     ) {
        DefaultPointListType tmp;
        extract(l,srf,eps,tmp);

    }







    template<class LevelSetType, class SurfaceInserterType, class ActivePointListType>
    void extract(   const LevelSetType& l,                          //the level set function
                    SurfaceInserterType& srf,                       //the surface inserter type
                    typename LevelSetType::value_type eps,         //eps can be used to avoid degenerated triangles
                    ActivePointListType& pt_lst                     //returns for each node of the surface mesh the nearest active grid point
                                                                  //which can occur for the marching cubes algorithm
                                                                    //the gridedge-surface intersection points are always set
                                                                    //at least a distance eps
                                                                    //away from both grid points of that grid edge
                ) {

        //this function extracts an explicit surface from the level set function
        //IMPORTANT: this function requires that the level set function consists of at least two layers of defined grid points
        //           this can be checked using the level set member function number_of_layers()
        //           the level set function can be expanded to 2 layers of defined grid points using the
        //           member function expand(new_number_of_layers) with a value for new_numer_of_layers>=2
        //
        //the SurfaceInserterType has to be a class with following type definitions and member functions
        //
        //class SurfaceInserterType {
        //public:
        //
        //  typedef unsigned int node_ref_type;     //here as example an unsigned int is used for referencing the nodes
        //
        //  node_ref_type insert_node(const LevelSetType::value_type*);  //the coordinates of each new node is passed to
        //                                                               //this function. So the user can decide how to store
        //                                                               //the new node. The function returns a reference, which can be for
        //                                                               //example a pointer or an index, to the new node
        //
        //  void insert_element(const node_ref_type*);                   //the node references of a new element are passed to this function
        //                                                               //In this function the user defines how the connectivity information
        //                                                               //should be stored
        //};

      typedef typename LevelSetType::index_type index_type;
    typedef typename SurfaceInserterType::node_ref_type node_ref_type;


      //const int       next[12]      ={-1, 3, 0,-1, 2, 7, 4, 1,-1, 8,11, 9};
      //const int       previous[12]  ={ 2, 7, 4, 1, 6,-1,-1, 5, 9,11,-1,10};
      /*const index_type    shift[12][3]={
          { 0, 0, 0},  // 0
          { 1, 0, 0},  // 1
          { 0, 1, 0},  // 2
          { 0, 0, 0},  // 3
          { 0,-1, 1},  // 4
          { 1, 0, 0}, // 5
          { 0, 1, 0}, // 6
          {-1, 0, 1}, // 7
          { 0, 0, 0}, // 8
          { 1, 0, 0}, // 9
          { 1, 0, 0}, //10
          {-1, 1, 0}  //11
      };*/


      const unsigned int corner0[12]   ={0,1,2,0,4,5,6,4,0,1,3,2};
      const unsigned int corner1[12]   ={1,3,3,2,5,7,7,6,4,5,7,6};
      const unsigned int direction[12] ={0,1,0,1,0,1,0,1,2,2,2,2};


        bool parity=l.grid().parity();

        assert(l.number_of_layers()>=2);       //test if level set function consists of at least 2 layers of defined grid points

        typedef typename LevelSetType::value_type value_type;

        const int D=LevelSetType::dimensions;

        typedef typename std::map<vec<index_type, D>, node_ref_type> node_container;

        //typedef typename std::vector<typename std::deque<node_ref_type> > node_container2;

        node_container nodes[D];

        //node_container2 nodes2((D==2)?4:12);

        typename node_container::iterator it;

//      SET BUT NOT USED WARNING:
//        index_type old_idx[D];
//        for (int i=0;i<D;++i) old_idx[i]=l.grid().min_point_index(i);

        //iterate over all active points
        for (typename LevelSetType::template const_iterator_cells_filtered<typename LevelSetType::filter_all> it_c(l);
            !it_c.is_finished();it_c.next()) {

          //assert(l.grid().is_cell_member(it_c.indices()));

          //std::cout << "cell_index = " << it_c.indices() << std::endl;


          for (int u=0;u<D;u++) {
        while (!nodes[u].empty() && nodes[u].begin()->first < vec<index_type,D>(it_c.indices())) nodes[u].erase(nodes[u].begin());
          }

            unsigned int signs=0;
            for(int i=0;i<(1<<D);i++) {
                if(it_c.corner(i).sign()==POS_SIGN) signs|=(1<<i);
            }

            if (signs==0) continue;
            if (signs==(1<<(1<<D))-1) continue;

            //std::vector<bool> point_on_edge((D==2)?4:12, false);

            //for each element
            for (  const int* Triangles =marching_cubes::polygonize<D>(signs);
                Triangles[0]!=-1;
                Triangles+=D) {

              vec<node_ref_type,D> nod_numbers;
              vec<node_ref_type,D> nod_numbers2;

                //for each node
                for (int n=0;n<D;n++) {

                  const int edge=Triangles[n];

                    //find grid points of corresponding edge
                    /*unsigned int p0,p1;
                    for (p0=0;((Triangles[e*D+n] & (1<<p0))==0);p0++);
                    for (p1=p0+1;((Triangles[e*D+n] & (1<<p1))==0);p1++);*/

                    unsigned int p0=corner0[edge];
                    unsigned int p1=corner1[edge];

                    //std::cout << "edge = " << p0 << " " << p1 << std::endl;

                    //determine direction of edge
                    int dir=direction[edge];
                    //while((p1-p0)!=(1u<<dir)) dir++;

                    //look for existing surface node
                    vec<index_type,D> d(it_c.indices());
                    d+=get_corner<D,index_type>(p0);

                    it=nodes[dir].find(d);
                    if (it!=nodes[dir].end()) {

                      //assert(!point_on_edge[edge]);

                        nod_numbers[n]=it->second;

            //int e=previous[edge];
            //vec<index_type,D> v=it_c.indices();

            /*while (e!=-1) {
              v-=vec<index_type,D>(shift[e]);
              if (l.grid().is_cell_member(v)) {
                assert(!nodes2[edge].empty());
                nod_numbers2[n]=nodes2[edge].front();
                nodes2[edge].pop_front();
                break;
              }
              e=previous[e];
            }
            assert(e!=-1);*/
            //assert(!nodes2[edge].empty());
            //nod_numbers2[n]=nodes2[edge].front();
            //nodes2[edge].pop_front();

            //std::cout << nod_numbers[n] << " " << nod_numbers2[n] << std::endl;

            //assert(nod_numbers[n]==nod_numbers2[n]);

            //point_on_edge[edge]=true;

                    } else {    //if node does not exist yet

                      //assert(point_on_edge[edge]==false);

                        //calculate coordinate of new node
                        vec<value_type, D> cc;
                        for (int z=0;z<D;z++) {
                            if (z!=dir) {
                                cc[z]=value_type(it_c.indices(z)+get_corner<D,index_type>(p0)[z]);
                            } else {
                                value_type d0, d1;

                                d0=it_c.corner(p0).value();
                                d1=it_c.corner(p1).value();

                                //calculate the surface-grid intersection point
                                if (d0==-d1) { //includes case where d0=d1=0
                                    pt_lst.push_back(it_c.corner(p0).active_pt_id());
                                    cc[z]=  static_cast<value_type>(it_c.indices(z))+0.5;
                                } else {
                                    if (math::abs(d0)<=math::abs(d1)) {
                                        pt_lst.push_back(it_c.corner(p0).active_pt_id());
                                        cc[z]=  static_cast<value_type>(it_c.indices(z))+(d0/(d0-d1));
                                    } else {
                                        pt_lst.push_back(it_c.corner(p1).active_pt_id());
                                        cc[z]=  static_cast<value_type>(it_c.indices(z)+1)-(d1/(d1-d0));
                                    }
                                }

                                cc[z]=std::max(cc[z], it_c.indices(z)+eps);
                                cc[z]=std::min(cc[z], (it_c.indices(z)+1)-eps);

                            }
                            cc[z]=l.grid().global_index_2_global_coordinate(z, cc[z]);
                        }

                        //insert new node
                        nod_numbers[n]=srf.insert_node(&cc[0]);             //insert new surface node
                        nodes[dir][d]=nod_numbers[n];

                        //nod_numbers2[n]=nod_numbers[n];
                        //nodes2[edge].push_front( nod_numbers2[n]);

                        //point_on_edge[edge]=true;


                    }
                    //assert(point_on_edge[edge]);
                    //assert(nod_numbers[n]==nodes2[edge].front());
                }

                if (parity) std::swap(nod_numbers[0], nod_numbers[1]);      //change the orientation of the elements in dependence of the grid parity

                srf.insert_element(&(nod_numbers[0]));                      //insert new surface element

            }

    }
    }

    ///
    ///     This function takes an empty LevelSet and fills it with a box with lower corner at start and upper corner at end
    ///
    template <class LevelSetType>
    void make_box(LevelSetType &LS,
            vec<typename LevelSetType::index_type,LevelSetType::dimensions> start,
            vec<typename LevelSetType::index_type,LevelSetType::dimensions> end){

        typedef typename LevelSetType::index_type index_type;
        typedef typename LevelSetType::value_type value_type;
        static constexpr int D=LevelSetType::dimensions;

        // Add borderpoints to levelset before extracting the surface
        std::vector< std::pair< vec<index_type, D>, value_type> > points;

        //allocate enough space for all points
        index_type num_points=0;
        for(int i=0; i<D; ++i) num_points+=4.2*std::abs(start[i]-end[i])* ((D==3)?std::abs(start[(i+1)%D]-end[(i+1)%D]):1);
        points.reserve(num_points);      //resize to stop many reallocs 4.2 instead of 4 for extra points

        // add points on simulation boundary
        std::vector< vec<index_type,D> > unity;
        for(unsigned i=0; i<D; ++i){ // fill with 1,0,0; 0,1,0; 0,0,1
          vec<index_type,D> temp = vec<index_type,D>(0);
          temp[i] = 1;
          unity.push_back(temp);
        }

        // slight offset for numerical stability
        double eps=1e-6;

        for(int i=0; i<D; ++i){
            vec<index_type,D> index;
            int y=(i+1)%D, z=(i+2)%D;   //permutation of other dimensions

            int jmin = start[y], jmax = end[y];
            int kmin = (D==3)?start[z]:0, kmax = (D==3)?end[z]:2; //in 2D, make loop run only once

            for(int j=jmin+1; j<jmax; ++j){
                for(int k=kmin+1; k<kmax; ++k){
                    index[y] = j;
                    index[z] = k;

                    index[i] = start[i];
                    points.push_back(std::make_pair(index+unity[i], -1.+eps));
                    points.push_back(std::make_pair(index, 0.+eps));

                    index[i] = end[i];
                    points.push_back(std::make_pair(index-unity[i], -1.+eps));
                    points.push_back(std::make_pair(index, 0.+eps));
                }
            }
        }

        //sort and remove degenerate points(-1 points on inner corner)
        std::sort(points.begin(), points.end());
        typename std::vector< std::pair< vec<index_type, D>, value_type> >::iterator it = std::unique(points.begin(), points.end());

        points.resize(std::distance(points.begin(), it));

        //put into levelset
        LS.insert_points(points);
    }

}

#endif /*LEVELSET2SURFACE_HPP_*/
