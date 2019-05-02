#ifndef DEF_NARROWBAND_H
#define DEF_NARROWBAND_H

/* =========================================================================
   Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.

                            -----------------
                 ViennaTS - The Vienna Topography Simulator
                            -----------------

   Contact:         viennats@iue.tuwien.ac.at

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */

//IF PARALLEL_MODE IS ACTIVE

#include "Cells.h"
#include "Time.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Math.h"
#include <algorithm>
#include <vector>
#include <numeric>
#include "LSlib/vector.hpp"
#include "message.h"

#include "boundaries.h"

#include <cmath>
#include <cassert>
#include "LSlib/math.hpp"

///Namespace for calculation helpers.
namespace calc {

    template <int Dimensions> class Make3DVector {
        const double* v;
    public:
        inline double operator[](int i) const {
            return (i<Dimensions)?v[i]:0.;
        }
        Make3DVector(const double* v2):v(v2) {}

    };


    template <class ParameterType> class PartitionTraits {
    public:
        typedef int IntCoordType;
        typedef unsigned int IntLinkType;
        //typedef const geom::cell<D>* CellRefType;
        typedef unsigned int CellRefType;

        static inline CellRefType UndefinedCellRef() {  return std::numeric_limits<CellRefType>::max();}//return 0;}
        static inline IntLinkType UndefinedLink() { return std::numeric_limits<IntLinkType>::max();}

        //static const double SurfaceAreaHeuristicLambda=ParameterType::SurfaceAreaHeuristicLambda;
        //static const partition::SplittingModeType PartitionMode=ParameterType::PartitionMode;

        static const int Dimension=ParameterType::Dimension;

    };

  template<int D, class LS, class NormalVectorVectorClass, class DistancesVectorClass>
  void SetupCells(   const LS& l,geom::cells<D> &Cells,
                        std::vector<lvlset::vec<int,D> > &CellCoordinates,
                        const NormalVectorVectorClass& NormalVectors,
                        const DistancesVectorClass& DistancesToReceiver,
                        double RecepterRadius
                     ) {

        Cells.clear();

        //int cell_counter=0;
        //int cell_equal_signs_counter=0;
        //int cell_contains_disk_counter=0;
        //int cell_inserted_counter=0;

        for (typename LS::template const_iterator_cells_filtered<typename LS::filter_active> it(l);!it.is_finished();it.next()) {

            //cell_counter++;

            bool cell_contains_disk=false;

            int sgn_count=0;
            for (int i=0;i<(1<<D);i++) {
                sgn_count+=it.corner(i).sign();
            }

            if ((sgn_count!=(1<<D)) && (sgn_count!=0)) {
                cell_contains_disk=true;
//                //std::cout << "cell_cont_disk! ";
            } else {

                //cell_equal_signs_counter++;

                //check all corners
                for (int i=0;i<(1<<D);++i) {

                    //if corner is active
                    if (it.corner(i).is_active()) {

                        unsigned int id=it.corner(i).active_pt_id();

                        const double &d=DistancesToReceiver[id];

                        //check if disk of active grid point intersects corresponding cell

                        //check for all dimensions
                        bool cell_contains_disk2=true;
                        for (int dir=0;dir < D;++dir) {

                            const double &n=NormalVectors[id*D+dir];

                            double min=n*d+((i>>dir) & 1);
                            double max=min;
                            double tmp=std::sqrt(std::max(0.,1-n*n))*RecepterRadius;
                            min-=tmp;
                            max+=tmp;
                            if ((max<0.) || (min>1.) ) {
                                cell_contains_disk2=false;
                                break;
                            }
                        }

                        if (cell_contains_disk2) {
//                            //std::cout << "cell_cont_disk2! ";
                            cell_contains_disk=true;
                            //cell_contains_disk_counter++;
                            break;
                        }
                    }
                }
            }



            if (cell_contains_disk) {

                //cell_inserted_counter++;

                lvlset::vec<unsigned int, (1<<D)> points;

                CellCoordinates.push_back(lvlset::vec<int,D>(it.indices()));

                Cells.push_back(geom::cell<D>());

                for (unsigned int i=0;i<(1<<D); i++ ) {
//                  if (it.corner(i).pt_id()==1028) //std::cout << "it.corner("<<i<<").pt_id():" <<it.corner(i).pt_id() <<std::endl;
                    Cells.back().Points[i]=it.corner(i).pt_id();
                    assert(it.corner(i).is_defined());
                }
            }
        }

        ////std::cout << "num_cells=" << cell_counter << std::endl;
        ////std::cout << "num_cells_equal_signs=" << cell_equal_signs_counter << std::endl;
        ////std::cout << "num_cells_contains_disk=" << cell_contains_disk_counter << std::endl;
        ////std::cout << "num_cells_inserted=" << cell_inserted_counter << std::endl;
    }

  template <class LS> void CalculateNormalVectors(
      const LS& l,
      std::vector<double>& NormalVectors,
      std::vector<double>& DistancesToReceiver,
      int open_boundary_direction,
      bool is_open_boundary_negative,
      double ReceptorRadius,
      const lvlset::vec<double,LS::dimensions> & default_directions=lvlset::vec<double,LS::dimensions>(0)) {

        const int D=LS::dimensions;

        lvlset::vec<double,D> t=default_directions;
        double tmp=Norm(t);
        if (tmp==0) {
          t[open_boundary_direction]=(is_open_boundary_negative)?-1.:1.;
        } else {
          t/=tmp;
        }

        NormalVectors.clear();
        NormalVectors.resize(l.num_active_pts()*D);

        DistancesToReceiver.clear();
        DistancesToReceiver.resize(l.num_active_pts(),0.);

        //!Calculate Normalvectors

        typename LS::points_type segmentation=l.get_new_segmentation();

        #pragma omp for schedule(static, 1) // parallelization - Iterations divided into chunks of size 1. Each chunk is assigned to a thread
        for (int p=0;p<= static_cast<int>(segmentation.size());++p) {
            typename LS::point_type  begin_v=(p==0)?l.grid().min_point_index():segmentation[p-1];
            typename LS::point_type  end_v=(p!=static_cast<int>(segmentation.size()))?segmentation[p]:l.grid().increment_indices(l.grid().max_point_index());

            for (typename LS::template const_iterator_neighbor_filtered<typename LS::filter_active,  1> it(l, typename LS::filter_active(), begin_v);it.indices()<end_v;it.next()) {

                double* n=&NormalVectors[it.center().active_pt_id2()*D];
                double& dist=DistancesToReceiver[it.center().active_pt_id2()];

                for (int i=0;i<D;i++) {
                  double pos=it.neighbor(i).value()-it.center().value();
                  double neg=it.center().value()-it.neighbor(i+D).value();
                  n[i]=0;
                  if ((pos > 0 && neg < 0) || (pos < 0 && neg > 0)) {
                    if (default_directions[i]<0) {
                      n[i]=std::min(neg,pos);
                    } else if (default_directions[i]>0) {
                      n[i]=std::max(neg,pos);
                    } else {
                      n[i]=(pos+neg)*0.5;
                    }
                  } else {
                    n[i]=(pos+neg)*0.5;
                  }
                }

                //for (int i=0;i<D;i++) n[i]=it.gradient2(i);

                double tmp_max=std::fabs(n[0]);
                for (int uu=1;uu<D;uu++) tmp_max=std::max(tmp_max,std::fabs(n[uu]));

                if (tmp_max==0.) {
                  for (int uu=0;uu<D;uu++) n[uu]=t[uu];
                  dist=0.;
                } else {
                    double no2=0.;
                    for (int uu=0;uu<D;uu++) no2+=my::math::pow2(n[uu]/tmp_max);
                    double no=tmp_max*std::sqrt(no2);
                    for (int uu=0;uu<D;uu++) n[uu]/=no;
                    dist=-it.center().value()/no;
                }

                //reduce distance if receptor disk would not completely lie inside of attached voxels

                for (int i=0;i<D;i++) {
                    dist=my::math::Sign(dist)*std::max(
                            0.,
                            std::min(
                                    std::fabs(dist),
                                    (1.-ReceptorRadius*std::sqrt(
                                                        std::max(0.,1-n[i]*n[i])
                                                    ))/std::fabs(n[i]))
                        );
                    assert(!std::isnan(dist));
                }

                //for (int i=0;i<D;++i) NormalVectors.push_back(n[i]);
                //DistancesToReceiver.push_back(dist);
            }
        }
    }

  template <class LS> void CalculateCurvatureVectors(const LS& l, std::vector<double>& CurvatureVectors, bool initialized) {

    ////std::cout << "here!\n";
        const int D=LS::dimensions;
        typedef typename LS::index_type index_type;

        CurvatureVectors.clear();
        CurvatureVectors.resize(l.num_active_pts());

        typename LS::points_type segmentation=l.get_new_segmentation();
        std::vector<typename LS::const_iterator_runs_offset> it_neighbors;

    #pragma omp for schedule(static, 1) // parallelization - Iterations divided into chunks of size 1. Each chunk is assigned to a thread
    for (int p=0;p<= static_cast<int>(segmentation.size());++p) {
      typename LS::point_type  begin_v=(p==0)?l.grid().min_point_index():segmentation[p-1];
      typename LS::point_type  end_v=(p!=static_cast<int>(segmentation.size()))?segmentation[p]:l.grid().increment_indices(l.grid().max_point_index());

      for (typename LS::template const_iterator_neighbor_filtered<typename LS::filter_active,  1> it(l, typename LS::filter_active(), begin_v);it.indices()<end_v;it.next()) {

              double* curv=&CurvatureVectors[it.center().active_pt_id2()];

        if (initialized) {
          for (unsigned int i=0;i<it_neighbors.size();i++) it_neighbors[i].go_to_indices_sequential(it.indices());
        } else {
          for (int i=-1;i<=1;++i) {
            for (int j=-1;j<=1;++j) {
              for (int k=-1;k<=1;++k) {
                if (((i!=0) || (j!=0) || (k!=0)) && ((i==0) || (j==0) || (k==0))) {
                  lvlset::vec<index_type,D> v(i,j,k);
                  ////std::cout << "v(" << v[0] << ", " << v[1] << ", " << v[2] << ")\n";
                  it_neighbors.push_back(typename LS::const_iterator_runs_offset(l, v,it.indices()));
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

              double PhiXX=it_neighbors[XpY0Z0].value()+it_neighbors[XmY0Z0].value()-2*it.center().value();
              double PhiYY=it_neighbors[X0YpZ0].value()+it_neighbors[X0YmZ0].value()-2*it.center().value();
              double PhiZZ=it_neighbors[X0Y0Zp].value()+it_neighbors[X0Y0Zm].value()-2*it.center().value();

              double PhiXY=(it_neighbors[XpYpZ0].value()+it_neighbors[XmYmZ0].value()-it_neighbors[XpYmZ0].value()-it_neighbors[XmYpZ0].value())*0.25;
              double PhiXZ=(it_neighbors[XpY0Zp].value()+it_neighbors[XmY0Zm].value()-it_neighbors[XpY0Zm].value()-it_neighbors[XmY0Zp].value())*0.25;
              double PhiYZ=(it_neighbors[X0YpZp].value()+it_neighbors[X0YmZm].value()-it_neighbors[X0YmZp].value()-it_neighbors[X0YpZm].value())*0.25;

              //const int mode=0;

              double denom=PhiX*PhiX+PhiY*PhiY+PhiZ*PhiZ;

              double num=     0.5*PhiX*PhiX*(PhiYY+PhiZZ)-PhiY*PhiZ*PhiYZ+        //mean curvature
                              0.5*PhiY*PhiY*(PhiXX+PhiZZ)-PhiX*PhiZ*PhiXZ+
                              0.5*PhiZ*PhiZ*(PhiXX+PhiYY)-PhiX*PhiY*PhiXY;

//              double s=0.;

              //if (material<material_level) {

                  if (denom!=0) {
                    *curv=num/(denom*std::sqrt(denom));
                      ////std::cout << "*curv: " << *curv << std::endl;
//                      if ((k>max_curvature) || (k<min_curvature)) s= -num/denom;

                  } else {
                      //std::cout << "warning!!!dlkajf" << std::endl;
                      if (num>0) {
                        *curv=-std::numeric_limits<double>::max();
                      } else {
                        *curv=std::numeric_limits<double>::max();
                      }
                  }
              //}
//              return s;
//                  //std::cout << "curv: " << *curv << "\n";
      }
        }

        //!Calculate Curvature vectors

    }

  namespace {

        template<int D, class PartitionType>
        class ClusterPositionType {
        public:
            double X[D];
            typename PartitionType::subbox Subbox;
        };

  }

  // [Josef] Main function where particles are tracked and collision with interface is checked!
  template <class ModelType, class ParameterType, class PartitionType, class LevelSetType> void CalculateRates(
    const ModelType &Model,
    const ParameterType &Parameter,
    const PartitionType &Partition,
    const LevelSetType &SurfaceLevelSet,
    const std::vector<double>& NormalVectors,
    const std::vector<double>& DistancesToReceiver,
    const std::vector<double>& Coverages,
    std::vector<double>& Rates,
    const std::vector<unsigned int>& PointMaterials,
    const geom::cells<ParameterType::Dimension>& Cells,
    double ProcessTime) {
//      std::cout << "1\n";
    const int D=ParameterType::Dimension;
    typedef ClusterPositionType<D, PartitionType> ClusterPositionType;


    const double ReceptorRadius=Parameter.receptor_radius;
    const double ReceptorRadius2=ReceptorRadius*ReceptorRadius;

    const double further_tracking_distance=Parameter.further_tracking_distance; //default is 3

    //Initialize Rates
    unsigned int  num_active_points=SurfaceLevelSet.num_active_pts();
    Rates.clear();
    if(NormalVectors.size()!=num_active_points*D){
        std::cout << "Assert normal vector size: " << NormalVectors.size() << " = " << num_active_points*D << "\n";
        assert(0);
    }
    if(Coverages.size()<num_active_points*Model.CoverageStorageSize){
        std::cout << "Assert Coverage size: " << Coverages.size() << " >= " << num_active_points*Model.CoverageStorageSize << "\n";
        assert(0);
    }

        #ifdef _OPENMP
    const int max_threads=omp_get_max_threads();
        #else
    const int max_threads=1;
        #endif


    std::vector<std::vector<double > > all_tmp_Rates(
                        max_threads,
                        std::vector<double>(num_active_points*Model.RatesStorageSize,0.)
                      );

    double RecepterArea=(D==3)?(Parameter.receptor_radius*Parameter.receptor_radius*my::math::Pi):(2.*Parameter.receptor_radius);
    if (!ModelType::SpatiallyEqualDistributedFlux) {
        RecepterArea*=Parameter.grid_delta;
        if (D==3) RecepterArea*=Parameter.grid_delta;
    }
    ////std::cout << "Recepter Area: " << RecepterArea << endl;


    #pragma omp parallel
    {

        //determine the number of different starting locations
        //in case of equally distributed flux the number of starting places is equal to the open surface area measured in grid spacings
        //in case of non-equally distributed flux the number of starting places is set to the number of threads

        #ifdef _OPENMP
            const int my_num_threads=omp_get_num_threads();
            const int my_thread_num=omp_get_thread_num();
            #else
            const int my_num_threads=1;
            const int my_thread_num=0;
            #endif

        const int NumStartingPlaces=(ModelType::SpatiallyEqualDistributedFlux)?
                                            static_cast<int>(Partition.AreaSize(Parameter.open_boundary)):   //AXIS
                        my_num_threads;


            //for each thread a vector is defined, where the rates are stored
      std::vector<double>& tmp_Rates=all_tmp_Rates[my_thread_num];

      assert(tmp_Rates.size()==num_active_points*Model.RatesStorageSize);
      //std::cout << "assert temp rates size \n";

      //stacks to store the particles and their positions
      std::stack<typename ModelType::ParticleType> ParticleStack;
      std::stack<ClusterPositionType> ParticlePositionsStack;


      //beginning of parallel section with dynamic scheduling
      #pragma omp for schedule(dynamic) //Chunks are dynamically assigned to threads on a first-come, first-serve basis as threads become available.
      for (int StartingPlace=0;StartingPlace<NumStartingPlaces;++StartingPlace) {  //for each starting place do

          //used to store the partition subbox the particles starts from
          typename PartitionType::subbox starting_subbox;

          //the start position of the particle (in global coordinates)
          double StartPosition[3];

          //if spatially equal distributed flux, determine the start_box
          if (ModelType::SpatiallyEqualDistributedFlux) {

          unsigned int tmp_s=StartingPlace;
          int tmp_dim=Parameter.open_boundary;    //AXIS
          if (tmp_dim==0) tmp_dim=D;
          --tmp_dim;

          for (int i=0;i<D-2;++i) {
            StartPosition[tmp_dim]=tmp_s%Partition.Extension(tmp_dim);
            tmp_s/=Partition.Extension(tmp_dim);
            if (tmp_dim==0) tmp_dim=D;
            --tmp_dim;
          }
          StartPosition[tmp_dim]=tmp_s;

          starting_subbox=Partition.Access(StartPosition, Parameter.open_boundary, Parameter.open_boundary_negative);
          //std::cout << "equaldistributed \n";

          }


           //for each involved particle type do
          for (unsigned int ParticleType=0;ParticleType<Model.NumberOfParticleTypes;++ParticleType) {
              //if (ParticleType==1) //std::cout << "AH!!!\n";

              //determine the number of particles which have to be simulated
             const unsigned int NumOfParticles=(ModelType::SpatiallyEqualDistributedFlux)?
                            Model.NumberOfParticleClusters[ParticleType]:
                            Model.NumberOfParticleClusters[ParticleType]/my_num_threads
                                +(my_thread_num < static_cast<int>(Model.NumberOfParticleClusters[ParticleType]%my_num_threads)?1:0);

            //for each particle do
            for (unsigned int ParticleCounter=0;ParticleCounter<NumOfParticles;++ParticleCounter) {
                      //std::cout << "\nparticles\n";


            //generate cluster energy and direction
            typename ModelType::ParticleType p;
            //typename ModelType::TipHeightType dist;

            Model.ParticleGeneration(p,ParticleType,ProcessTime, StartPosition);
            //std::cout << "\nparticlegeneration\n";
            //if particle is not moving downwards
            if (Parameter.open_boundary_negative) {
                if(p.Direction[Parameter.open_boundary]<=0.) continue;
            } else {
                if(p.Direction[Parameter.open_boundary]>=0.) continue;
            }

            //calculate represented flux by that particle
            p.Flux/=Model.NumberOfParticleClusters[ParticleType];
//            //std::cout<<"p.Flux1="<<p.Flux<<"\n";
            p.Flux/=RecepterArea;
//            //std::cout<<"p.Flux2="<<p.Flux<<"\n";

            //determine starting position and starting subbox
            ClusterPositionType cp;

            if (ModelType::SpatiallyEqualDistributedFlux) {
                //if flux is equal distributed

                //chose random start position
                for (int i=0;i<D;++i) {
                    cp.X[i]=StartPosition[i];
                    if (i!=Parameter.open_boundary) cp.X[i]+=my::stat::RandomNumber();
                }

                cp.Subbox=starting_subbox;


                //determine additional particles, which are necessary to account for extended boundaries
                int zmax[D-1];

                int dir=Parameter.open_boundary;

                for (int i=0;i<D-1;++i) {
                    dir=(Parameter.open_boundary+i+1)%D;
                    zmax[i]=0;
                    if (dir!=Parameter.open_boundary) {
                        if ((Parameter.boundary_conditions[dir].min==bnc::EXTENDED_BOUNDARY) && (p.Direction[dir]>0)) {
                            zmax[i]=static_cast<int>(
                                                                    std::ceil(
                                                                                (
                                                                                        -std::min(   std::fabs((Partition.Extension(Parameter.open_boundary)*p.Direction[dir])/p.Direction[Parameter.open_boundary]),
                                                   static_cast<double>(Parameter.max_extended_starting_position)
                                                )
                                                                                       -(cp.X[dir]+(cp.Subbox.Min(dir)-Partition.Min(dir)))
                                                                                )/Partition.Extension(dir)
                                                                             )
                                                                );
                            assert(zmax[i]<=0);
                                    } else if ((Parameter.boundary_conditions[dir].max==bnc::EXTENDED_BOUNDARY) && (p.Direction[dir]<0)) {
                                        zmax[i]=static_cast<int>(
                                                                    std::floor(
                                                                                (
                                                                                        std::min(   std::fabs((Partition.Extension(Parameter.open_boundary)*p.Direction[dir])/p.Direction[Parameter.open_boundary]),
                                                  static_cast<double>(Parameter.max_extended_starting_position)
                                                )
                                                                                        -(cp.X[dir]+(cp.Subbox.Min(dir)-Partition.Min(dir)))
                                                                                )/Partition.Extension(dir)
                                                                             )
                                                                )+1;
                                        assert(zmax[i]>=0);
                                    }
                                }
                }

                int counter[D-1];
                for (int k=0;k<D-1;++k) counter[k]=0;

                //add additional particles to the stack
                while (true) {

                    int h=0;
                    for (;h<D-1;++h) {
                        if (counter[h]!=zmax[h]) {
                            if (zmax[h]>0) ++counter[h]; else --counter[h];
                            break;
                        } else {
                            counter[h]=0;
                        }
                    }
                    if (h==D-1) break;


                    ClusterPositionType new_cp;

                    for (int g=0;g<D-1;++g) {
                        int dir=(g+Parameter.open_boundary+1)%D;
                        new_cp.X[dir]=  cp.X[dir]+
                                                    static_cast<double>(cp.Subbox.Min(dir)-Partition.Min(dir))+                         //TODO check!!!
                                                    static_cast<double>(counter[g])*static_cast<double>(Partition.Extension(dir));
                    }
                    new_cp.Subbox=Partition.Access(new_cp.X, Parameter.open_boundary, Parameter.open_boundary_negative);
                    ParticlePositionsStack.push(new_cp);
                    ParticleStack.push(p);
                }
                //std::cout << "again EDF\n";
            }
            else {

                for (int i=0;i<D;++i) cp.X[i]=StartPosition[i]/Parameter.grid_delta;       //scale starting position

                double t=-(  cp.X[Parameter.open_boundary]-
                    ((Parameter.open_boundary_negative)?Partition.Min(Parameter.open_boundary):Partition.Max(Parameter.open_boundary))
                  )/p.Direction[Parameter.open_boundary];

                //Move cp.X to the top LS surface and update horizontal axis values (not open boundary value)
                for (int dir=0;dir<D;++dir) {
                    if (dir!=Parameter.open_boundary) {
                        bool ReverseSign;
                        cp.X[dir]=Parameter.boundary_conditions[dir].map_coordinate(Partition.Min(dir), Partition.Max(dir),cp.X[dir]+p.Direction[dir]*t, ReverseSign);
                        if (ReverseSign) p.Direction[dir]=-p.Direction[dir];
                    }
                }
                cp.Subbox=Partition.Access(cp.X, Parameter.open_boundary, Parameter.open_boundary_negative);
                            //cp.X is now position within subbox after removing the "global components"
            }

            //loop until particle stack is empty
            while (true) {
              //initialize the travelled distance from the intersection with -oo
              double travelled_distance_from_intersection(-std::numeric_limits<double>::max());
//std::cout << "DTFI\n";
              //the indices of the surface grid cell which was previously visited
              int last_surface_cell_indices[D];
              for (int r=0;r<D;++r) last_surface_cell_indices[r]=Partition.Min(r)-2;   //initialize with invalid indices


              //Iterate through the cells between LS.Max() and LS.Min() until surface reached or particle exits environment
              while (true) {
//std::cout << "particleIteration\n";
                //get reference to actual cluster
                const typename PartitionType::subbox &Subbox= cp.Subbox;

                //#######################################################
                //# find max distance within box                        #
                //#######################################################
                double max_distance_in_box=std::numeric_limits<double>::max();

                int LeavingDirection=-1;                  //LeavingDirection : 0,1,2 particle leaves box in x,y,z direction respectively

                std::bitset<D> PositionStatusPos;       //for each direction the bit is set if a particle is outside (in positive direction) of the regular simulation domain
                std::bitset<D> PositionStatusNeg;       //for each direction the bit is set if a particle is outside (in negative direction) of the regular simulation domain


                //for each dimension do
                int i;
                for (i=0;i<D;i++) {
//std::cout << "Dimension\n";
                  double t_temp=std::numeric_limits<double>::max();

                  //Subbox.Extension(dir) is the length of the subbox in the dir direction
                  //Subbox.Min(dir) is the global coordinate (grid points) of the Subbox edge in the min dir direction
                  //Subbox.Max(dir) is the global coordinate (grid points) of the Subbox edge in the max dir direction
                  //Global coordinate is then found by Subbox.Min(dir)+cp.X[dir]

                  //When outside the min extended boundary
                  if ((cp.X[i]<=0) && (Parameter.boundary_conditions[i].min==bnc::EXTENDED_BOUNDARY)) {
                      if (p.Direction[i]>0.) {
                          if (cp.X[i]==0) {
                              t_temp=(Subbox.Extension(i)-cp.X[i])/p.Direction[i];
                          } else {
                              t_temp=-cp.X[i]/p.Direction[i];
                              PositionStatusNeg.set(i);
                          }
                      } else {
                          if (cp.X[i]<0/*-Parameter.DomainExtension*/) break;
                          PositionStatusNeg.set(i);
                      }
                  //When outside the max extended boundary
                  } else if ((cp.X[i]>=Subbox.Extension(i)) && (Parameter.boundary_conditions[i].max==bnc::EXTENDED_BOUNDARY)) {
                      if (p.Direction[i]<0.) {
                          if (cp.X[i]==Subbox.Extension(i)) {
                              t_temp=-cp.X[i]/p.Direction[i];
                          } else {
                              t_temp=(Subbox.Extension(i)-cp.X[i])/p.Direction[i];
                              PositionStatusPos.set(i);
                          }
                      } else {
                          if (cp.X[i]>Subbox.Extension(i)/*+Parameter.DomainExtension*/) break;
                          PositionStatusPos.set(i);
                      }
                  } else {
                      if (p.Direction[i]>0.) {
                        //t_temp is the variable to determine time to reach Subbox.Extension(i) from cp.X in p.Direction[i]
                          t_temp=(Subbox.Extension(i)-cp.X[i])/p.Direction[i];
                      } else if (p.Direction[i]<0.) {
                          t_temp=-cp.X[i]/p.Direction[i];
                      }
                  }
                  //Determine which axis is the leaving direction of the particle
                  if (t_temp<max_distance_in_box) {
                      max_distance_in_box=t_temp;
                      LeavingDirection=i;
                  }
                }

                //cp.X remains unchanged at this point, only max_distance_in_box and LeavingDirection are found
//                              //std::cout << "cp.X4 " << "(" << cp.X[0] << "," << cp.X[1] << "," << cp.X[2] << ")\n";
//                //std::cout << "max_dinstance_in_box: " << max_distance_in_box << "\n";
//                //std::cout << "LeavingDirection: " << LeavingDirection << "\n";
                //Now have max_distance_in_box = distance from cp.X to point of exit in the i direction
                //where i is the leaving direction (x,y,z)=(0,1,2)

                if ((i!=D) || (LeavingDirection==-1)) break;

                //When the subbox which received the particle also contains within it the surface boundary
                if (Subbox.ContainsCell()) {
//                  std::cout << "ContainsCell()\n";
                    //if subbox is a surface grid cell

                  const geom::cell<D> &Cell=Cells[Subbox.Cell()];
//std::cout << "containsCell\n";
                  //Calculate the exit direction and distance as before to see if surface is intersected
                  //#######################################################
                  //# check for surface intersection                      #
                  //#######################################################

                  //[Josef] This is where the particle is tracked within a subbox for surface intersection

                  //Check if Surface is intersected between position and position+max_distance_in_box*direction
                  if (travelled_distance_from_intersection==-std::numeric_limits<double>::max()) {

                    //get distances at corners
                    double Rho[1<<D];
                    int sgn_count=0;


                    for (int i=0;i<(1<<D);i++) {
                      Rho[i]=SurfaceLevelSet.value2(Cell.Points[((std::bitset<D>(i) | PositionStatusPos) & (~PositionStatusNeg)).to_ulong()]);
//std::cout << "Rho[" << i << "] = " << Rho[i] << "\n";
                      if (Rho[i]>0) sgn_count++;
                    }
//                    std::cout << "HERE!\n";

                    if (sgn_count!=(1<<D)) {

                      my::math::TransformRho2<D>::exec(Rho);

                                            double relative_distance_to_intersection;

                      if (sgn_count==0) {

                          relative_distance_to_intersection=0.;

                      } else {

                        my::math::Polynom<double, D> poly;

                        my::math::DetermineCoefficientsForImplicitRayTracing<D>(
                          cp.X,
                          p.Direction,
                          Rho,
                          &(poly.Coefficients()[0])
                        );
//                        std::cout << "cp.X5: (" << cp.X[0] << "," << cp.X[1] << "," << cp.X[2] << ")\n";
                        relative_distance_to_intersection=my::math::FindFirstTransitionFromPosToNegOfPolynomNewton(0., max_distance_in_box, poly,1e-6);
                      }

                      if (relative_distance_to_intersection < std::numeric_limits<double>::max()) {  //if particle hits surface
//std::cout << "hits\n";
                          travelled_distance_from_intersection=-relative_distance_to_intersection;

                          ClusterPositionType new_cp=cp;

                          for (int kk=0;kk<D;kk++) new_cp.X[kk]+=relative_distance_to_intersection*p.Direction[kk];
//                          std::cout << "RDI\n";
                          //p.Direction[kk] is unchanged at this point
                          //new_cp.X contains the coordinates within the Subbox at the subbox exit point

                          //determine normal vector
                          double tmp_normalvec[3];
                          my::math::CalculateNormal<D>(tmp_normalvec,new_cp.X,Rho);
                          if (D==2) tmp_normalvec[2]=0.;
//                                                std::cout << "tmp_normalvec(): " << tmp_normalvec[0] << ", " << tmp_normalvec[1] << ", " << tmp_normalvec[2] << "\n";

                          double dot=tmp_normalvec[0]*p.Direction[0];
                          for (int w=1;w<D;++w) dot+=tmp_normalvec[w]*p.Direction[w];

                          if (dot>=0.) {
                              msg::print_warning("Particle hits negative side of surface! Particle is skipped.");
                              break;
                          }

                          //calculate nearest active grid point to determine coverages and material (using Manhattan distance for speedup)
                          unsigned int gp=0;
                          unsigned int mat=0;
                          if ((ModelType::CoverageStorageSize>0) || (ModelType::ReemissionIsMaterialDependent)) {
//                              std::cout << "CSS>0 RMD\n";
                              double dist=std::numeric_limits<double>::max();
//                              std::cout << "D = " << D << "\n";
                              for (int g=0;g<(1<<D);g++) {
//                                  std::cout << "g = " << g << "\n";
                                  unsigned int tmp_gp= SurfaceLevelSet.active_pt_id(Cell.Points[g]);
//                                  std::cout << "tmp_gp = " << tmp_gp << "\n";
//                                  std::cout << "new_cp = " << new_cp.X[0] << ", " << new_cp.X[1] << ", " << new_cp.X[2] << "\n";
//                                  std::cout << "LevelSetType::INACTIVE = " << LevelSetType::INACTIVE << "\n";
                                  if (tmp_gp!=LevelSetType::INACTIVE) {
//                                      std::cout << "4\n";
                                      double tmp_dist=0;
//                                      std::cout << "5\n";
                                      for (int iii=0;iii<D;iii++) tmp_dist+=(((g & (1<<iii))==0)?(new_cp.X[iii]):(1.-new_cp.X[iii]));
//                                      std::cout << "6\n";
//                                          std::cout << "dist = " << dist << "\n";
//                                          std::cout << "tmp_dist = " << tmp_dist << "\n";
                                      if (tmp_dist<dist) {
//                                          std::cout << "7\n";
                                          dist=tmp_dist;
                                          gp=tmp_gp;
//                                          std::cout << "gp = " << gp << "\n";
                                      }
                                  }
                              }
//                              std::cout << "PointMaterials["<<gp<<"] = "<<PointMaterials[gp]<<"\n";
                              mat=PointMaterials[gp];
                          }

                          //perform particle reemission
//                          std::cout << "Reflection\n";
                          Model.ParticleReflexion(  p,
                                                      ParticleStack,
                                                      tmp_normalvec,
                                                      &Coverages[gp*Model.CoverageStorageSize],
                                                      mat//, D, dot
                                              );
                          while (ParticleStack.size()>ParticlePositionsStack.size()) ParticlePositionsStack.push(new_cp);
                      }
                    }
                  }
//std::cout << "positionStatusPos\n";
                  if (PositionStatusPos.none() && PositionStatusNeg.none()) {

                    //#####################################################################
                    //# determine corners which have to be checked for disk intersections #
                    //#####################################################################

                    std::bitset<(1<<D)> corners;

                    for (int dir=0;dir<D;++dir) {

                        switch(Subbox.Min(dir)-last_surface_cell_indices[dir]) {
                            case 0:
                                break;
                            case 1:
                                corners>>= (1<<dir);
                                corners|= ((dir<2)?((dir<1)?0xAA:0xCC):0xF0);
                                break;
                            case -1:
                                corners<<= (1<<dir);
                                corners|= ((dir<2)?((dir<1)?0x55:0x33):0x0F);
                                break;
                            default:
                                corners.set();
                        }
                    }

                    for (int s=0;s<D;++s) last_surface_cell_indices[s]=Subbox.Min(s);

                    //#######################################################
                    //# check for disk intersections                        #
                    //#######################################################


//[Josef] This is where the four corners of the box containing the particle are checked for intersection
                    //all 8 neighbors have to be checked if they are active and their disks are hit
                    for (int g=0;g<(1<<D);g++) {

                        if(corners.test(g)) {

                            unsigned int gp= SurfaceLevelSet.active_pt_id(Cell.Points[g]);
                            if (gp!=LevelSetType::INACTIVE) {

                                unsigned int gpD=gp*D;

                                double cos=-NormalVectors[gpD]*p.Direction[0];
                                for (int kk=1;kk<D;kk++) cos-=NormalVectors[gpD+kk]*p.Direction[kk];

                                if (cos > 0.) {

                                    //calculate relative position to disk midpoint
                                    double rel_pos[D];
                                    for (int kk=0;kk<D;kk++) rel_pos[kk]=cp.X[kk]-((g>>kk) & 1)-NormalVectors[gpD+kk]*DistancesToReceiver[gp];
                                    //rel_pos holds cp.X not cp_new.X

                                    //calculate rel_pos*disk_normal
                                    double rel_pos_dot_normal=rel_pos[0]*NormalVectors[gpD];
                                    for (int kk=1;kk<D;kk++) rel_pos_dot_normal+=rel_pos[kk]*NormalVectors[gpD+kk];


                                    if ( rel_pos_dot_normal <= (further_tracking_distance-travelled_distance_from_intersection)*cos ) {

                                        double tmpx=my::math::pow2(rel_pos[0]*cos+p.Direction[0]*rel_pos_dot_normal);
                                        for (int kk=1;kk<D;kk++) tmpx+=my::math::pow2(rel_pos[kk]*cos+p.Direction[kk]*rel_pos_dot_normal);

                                        if (tmpx<=cos*cos*ReceptorRadius2) {

                                            int Factor=1;
                                            for (int kk=0;kk<D;kk++) {
                                                if (kk!=Parameter.open_boundary) {
                                                    if ((Parameter.boundary_conditions[kk].min==bnc::REFLECTIVE_BOUNDARY) || (Parameter.boundary_conditions[kk].min==bnc::EXTENDED_BOUNDARY)) {
                                                        if ((g & (1<<kk))==0) {
                                                            if (Partition.Min(kk)==Subbox.Min(kk)) {
                                                                if (cp.X[kk]*cos+p.Direction[kk]*rel_pos_dot_normal>=0.) {
                                                                    Factor<<=1;
                                                                } else {
                                                                    Factor=0;
                                                                    break;
                                                                }
                                                            }
                                                        } else {
                                                            if (Partition.Max(kk)==Subbox.Max(kk)) {
                                                                if (cp.X[kk]*cos+p.Direction[kk]*rel_pos_dot_normal<=cos) {
                                                                    Factor<<=1;
                                                                } else {
                                                                    Factor=0;
                                                                    break;
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }

                                            for (;Factor>0;--Factor) {
//[Josef] Here, the particle has collided with the surface so the model's function to deal with this is called.


//                                                int mat = 0;
//                                                if ((ModelType::CoverageStorageSize>0) || (ModelType::ReemissionIsMaterialDependent))
//                                                  mat = PointMaterials[gp];
//std::cout << "Collision\n";
                                                Model.ParticleCollision(  p,
                                                                            Make3DVector<D>(&NormalVectors[gpD]),
                                                                            //&NormalVectors[gpD],
                                                                            &(tmp_Rates[gp*Model.RatesStorageSize]),
                                                                            &(Coverages[gp*Model.CoverageStorageSize]),
                                                                            ProcessTime//,
//                                                                            mat
                                                );
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                  }
                }

                //Dealing with particles beyond simulation boundaries and boundary conditions:

                //######################################################################
                //# check if calculation of particle cluster trajectory can be stopped #
                //######################################################################

                if (travelled_distance_from_intersection!=-std::numeric_limits<double>::max()) {
                  travelled_distance_from_intersection+=max_distance_in_box;
                  if (travelled_distance_from_intersection>=further_tracking_distance) break;
                }

                //#######################################################
                //# calculate exit point                                #
                //#######################################################

                for (int kk=0;kk<D;kk++) {
                  if (kk!=LeavingDirection) cp.X[kk]+=p.Direction[kk]*max_distance_in_box;
                }

                if (PositionStatusNeg.test(LeavingDirection)) {                     //particle enters regular simulation domain from negative side
                    cp.X[LeavingDirection]=0;
                } else if (PositionStatusPos.test(LeavingDirection)) {              //particle enters regular simulation domain from positive side
                    cp.X[LeavingDirection]=cp.Subbox.Extension(LeavingDirection);
                } else {

                                    //#######################################################
                                    //# get next box                                        #
                                    //#######################################################

                                    int old_min=cp.Subbox.Min(LeavingDirection);

                                    bool IsDirectionPositive=(p.Direction[LeavingDirection]>=0.);

                                    if (Partition.GoToNeighborBox(cp.Subbox,cp.X,LeavingDirection, IsDirectionPositive)) {

                                        if (IsDirectionPositive) {
                                            if (Parameter.boundary_conditions[LeavingDirection].max==bnc::INFINITE_BOUNDARY) break;

                                            if (Parameter.boundary_conditions[LeavingDirection].max==bnc::REFLECTIVE_BOUNDARY) p.Direction[LeavingDirection]=-p.Direction[LeavingDirection];
                                            cp.X[LeavingDirection]=cp.Subbox.Extension(LeavingDirection);

                                        } else {
                                            if (Parameter.boundary_conditions[LeavingDirection].min==bnc::INFINITE_BOUNDARY) break;

                                            if (Parameter.boundary_conditions[LeavingDirection].min==bnc::REFLECTIVE_BOUNDARY) p.Direction[LeavingDirection]=-p.Direction[LeavingDirection];
                                            cp.X[LeavingDirection]=0;
                                        }

                                        last_surface_cell_indices[LeavingDirection]=Partition.Min(LeavingDirection)-2;

                                    } else {

                                        cp.X[LeavingDirection]=(IsDirectionPositive)?0:cp.Subbox.Extension(LeavingDirection);
                                        if ((Parameter.boundary_conditions[LeavingDirection].min==bnc::PERIODIC_BOUNDARY) && PositionStatusPos.none() && PositionStatusNeg.none()) {
                                            if (IsDirectionPositive) {
                                                if (old_min>=cp.Subbox.Min(LeavingDirection)) last_surface_cell_indices[LeavingDirection]-=Partition.Extension(LeavingDirection);
                                            } else {
                                                if (old_min<=cp.Subbox.Min(LeavingDirection)) last_surface_cell_indices[LeavingDirection]+=Partition.Extension(LeavingDirection);
                                            }
                                        }
                                    }
                }

              }

//              std::cout << ParticleStack.size() << std::endl;
              if (ParticleStack.empty()) break;

              //#######################################################
              //# retrieve particle from stack                        #
              //#######################################################

              p=ParticleStack.top();
              ParticleStack.pop();

              cp=ParticlePositionsStack.top();
              ParticlePositionsStack.pop();

            } // end while loop: until particle stack is empty
          }//end of particle loop
        }//end of particle type loop
      }

      #pragma omp single //run by a single available thread.
      {
        Rates.swap(all_tmp_Rates[0]);
      }

      #pragma omp for
      for (int i=0;i<static_cast<int>(Rates.size());i++) {
        for (int j=1;j<my_num_threads;j++) {
          Rates[i]+=all_tmp_Rates[j][i];
        }
      }

      // [josef] now that all thead-exclusive thread rates have been merged, we can output them
    if (Model.OutputFluxes) {
        {
        std::ofstream outputfile("rates.csv");
        for (typename LevelSetType::const_iterator_runs it(SurfaceLevelSet); !it.is_finished(); it.next())
        {
        if(it.active_pt_id() != LevelSetType::INACTIVE)
        {
          for (int j=0;j<D;j++) outputfile << (it.start_indices()[j]) << " ";
          outputfile << Rates[it.active_pt_id()] << std::endl;
        }
        }
        outputfile.close();
        }
        {
        std::ofstream outputfile("rates_griddelta.csv");
        for (typename LevelSetType::const_iterator_runs it(SurfaceLevelSet); !it.is_finished(); it.next())
        {
        if(it.active_pt_id() != LevelSetType::INACTIVE)
        {
          for (int j=0;j<D;j++) outputfile << (it.start_indices()[j])*Parameter.grid_delta << " ";
          outputfile << Rates[it.active_pt_id()] << std::endl;
        }
        }
        outputfile.close();
        }
    }
    }

    //local_time=my::time::GetTime()-StartTime;
  }

//  template<class ModelType> void UpdateCoverages(const std::vector<double>& Rates, std::vector<double>& Coverages, const ModelType& Model) {
  template<class ModelType> void UpdateCoverages(const std::vector<double>& Rates, std::vector<double>& Coverages,
          const ModelType& Model, double &time_step) {//, double &current_time) {
    double* c=&Coverages[0];
    const double* r=&Rates[0];
    while (r!=&(*(Rates.end()))) {
      Model.UpdateCoverage(c, r, time_step);//, current_time);
//                    //std::cout << "time_step = " << time_step << "\n";
//      else Model.UpdateCoverage(c, r);
      c+=Model.CoverageStorageSize;
      r+=Model.RatesStorageSize;
    }
  }

  template<class ModelType> void UpdateCoverages(const std::vector<double>& Rates, std::vector<double>& Coverages, const ModelType& Model) {
    double* c=&Coverages[0];
    const double* r=&Rates[0];
    while (r!=&(*(Rates.end()))) {
//      if (time_step != 0) Model.UpdateCoverage(c, r, time_step);
      Model.UpdateCoverage(c, r);
      c+=Model.CoverageStorageSize;
      r+=Model.RatesStorageSize;
    }
  }
}




#endif //DEF_NARROWBAND_H
