#ifndef PROCESS_H_
#define PROCESS_H_

/* =========================================================================
   Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.

                            -----------------
                 ViennaTS - The Vienna Topography Simulator
                            -----------------

   Contact:         viennats@iue.tuwien.ac.at

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


#include "Time.h"
#include "calc.h"
#include <vector>
#include <fstream>
#include <string>
#include <list>
#include <algorithm>
#include <iostream>

#define BOOST_NO_HASH

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
//#include </home/filipov/boost/boost_1_49_0/boost/graph/adjacency_list.hpp>
//#include </home/filipov/boost/boost_1_49_0/boost/graph/connected_components.hpp>

#include "message.h"
//#include "Model/ModelPlanarization.h"

#include "Partition/PartitionNeighborLinksArrays.h"
#include "Partition/PartitionUpDownLinkTree.h"
#include "Partition/PartitionFullGrid.h"

#include "../LSlib/vector.hpp"
#include "boundaries.h"

namespace proc {

        template <class LevelSetType> void AddLayer(std::list<LevelSetType>& LS, int num_layers) {

            for (int i=0;i<num_layers;++i) {
            LS.push_back(LS.back());
        }

            for (int i=0;i>num_layers;--i) {
            assert(LS.size()>=2);
            LS.erase((LS.end()--)--);
        }
        }

	template <class LevelSetsType> void DetermineTopMostLayer(
			const LevelSetsType& LS,
			std::vector<unsigned int>& PointMaterials) {

	    //this function determines the materials of the most top levelset

		typedef typename LevelSetsType::value_type LevelSetType;

		PointMaterials.clear();
		PointMaterials.resize(LS.back().num_active_pts());

		typename LevelSetType::points_type segmentation=LS.back().get_new_segmentation();

//		std::cout << "segmentation.size(): " << segmentation.size() << std::endl;
		#pragma omp for schedule(static, 1) // parallelization - Iterations divided into chunks of size 1. Each chunk is assigned to a thread
		for (int p=0;p<= static_cast<int>(segmentation.size());++p) {

			typename LevelSetType::point_type  begin_v=(p==0)?LS.back().grid().min_point_index():segmentation[p-1];
			typename LevelSetType::point_type  end_v=(p!=static_cast<int>(segmentation.size()))?segmentation[p]:LS.back().grid().increment_indices(LS.back().grid().max_point_index());

			//iterator necessary to access
			std::vector< typename LevelSetType::const_iterator_runs> ITs;
			for (typename LevelSetsType::const_iterator it=LS.begin();&(*it)!=&(LS.back());++it)  ITs.push_back(typename LevelSetType::const_iterator_runs(*it,begin_v));

			for (typename LevelSetType::const_iterator_runs it(LS.back(),begin_v );it.start_indices()<end_v;it.next()) {
				if (!it.is_active()) continue;

				const typename LevelSetType::value_type d=it.value2();

				int z=LS.size()-1;
				for (;z>0;z--) {
					ITs[z-1].go_to_indices_sequential(it.start_indices());
					if (d<ITs[z-1].value()) break;
				}

				PointMaterials[it.active_pt_id2()]=LS.size()-1-z;
			}
		}

/*
        #pragma omp for schedule(static, 1) // parallelization - Iterations divided into chunks of size 1. Each chunk is assigned to a thread
        for (int p=0;p<= static_cast<int>(segmentation.size());++p) {


            typename LevelSetType::point_type  begin_v=(p==0)?LS.back().grid().min_point_index():segmentation[p-1];
            typename LevelSetType::point_type  end_v=(p!=static_cast<int>(segmentation.size()))?segmentation[p]:LS.back().grid().increment_indices(LS.back().grid().max_point_index());

            //iterator necessary to access
            std::vector< typename LevelSetType::const_iterator_runs> ITs;
            for (typename LevelSetsType::const_iterator it=LS.begin();&(*it)!=&(LS.back());++it)  ITs.push_back(typename LevelSetType::const_iterator_runs(*it,begin_v));

            for (typename LevelSetType::const_iterator_runs it(LS.back(),begin_v );it.start_indices()<end_v;it.next()) {
                if (!it.is_active()) continue;

                const typename LevelSetType::value_type d=it.value2();

                int z=LS.size()-1;
                for (;z>0;z--) {
                    ITs[z-1].go_to_indices_sequential(it.start_indices());
                    if (d<ITs[z-1].value()) break;
                }

                PointMaterials[it.active_pt_id2()]=LS.size()-1-z;
            }
        } */
        }

        namespace {

        template <class I1, class I2>
            bool connected(const I1& it1, const I2& it2) {
            return (it1.sign()==it2.sign());
        }

        }


        template <class LStype> std::pair<unsigned int, unsigned int> CalculateConnectivities(
        const LStype& l,
        std::vector<bool>& Connectivities,
        bool is_open_boundary_negative) {

        const int D=LStype::dimensions;

        boost::adjacency_list <boost::setS, boost::vecS, boost::undirectedS> Graph;

        unsigned int num_components=0;
        //unsigned int total_number_of_runs=0;

        //allocate memory for component list
        std::vector<int> comp_lst[l.number_of_segments()][D+1];
        for (unsigned int sub=0;sub<l.number_of_segments();++sub) {
            for (int i = -1;i<D;++i) {
                comp_lst[sub][i+1].resize(l.number_of_runs(i,sub),-1);
                //total_number_of_runs+=l.number_of_runs(i,sub);
            }
        }


        bool is_first_run=true;
        int node_of_first_run=0;
        int node_of_last_run=0;

        //cycle through
        for (typename LStype::template const_iterator_neighbor_filtered<typename LStype::filter_all,1> it(l);!it.is_finished();it.next()) {

            int & tc = comp_lst[it.center().get_segment_num()][it.center().get_level()][it.center().run_type_position()];

            if (tc==-1) {
                for (int k=0;k<2*D;++k) {
                    const int & tn= comp_lst[it.neighbor(k).get_segment_num()][it.neighbor(k).get_level()][it.neighbor(k).run_type_position()];
                    if (tn!=-1) {
                        if (connected(it.center(),it.neighbor(k))) {
                            tc=tn;
                            break;
                        }
                    }
                }
            }
            if (tc==-1) {
                tc=num_components;
                boost::add_vertex(Graph);
                ++num_components;
            }

           for (int k=0;k<2*D;++k) {
                int & tn= comp_lst[it.neighbor(k).get_segment_num()][it.neighbor(k).get_level()][it.neighbor(k).run_type_position()];
                if (connected(it.center(),it.neighbor(k))) {
                    if (tn!=-1) {
                        if (tc!=tn) boost::add_edge(tc,tn,Graph);
                    } else {
                        tn=tc;
                    }
                }
            }

           if (is_first_run) {
               is_first_run=false;
               node_of_first_run=tc;
           }
           node_of_last_run=tc;
        }

        assert(boost::num_vertices(Graph)==num_components);
        std::vector<int> component(boost::num_vertices(Graph));

        unsigned int num_components_after = connected_components(Graph, &component[0]);

        //determine component number of source region

        int source_node=(is_open_boundary_negative)?component[node_of_first_run]:component[node_of_last_run];

        Connectivities.clear();
        for (typename LStype::template const_iterator_neighbor_filtered<typename LStype::filter_active,1> it(l);!it.is_finished();it.next()) {
            if (it.center().sign()==lvlset::POS_SIGN) {
                assert(it.center().get_level()==0);
                assert(it.center().get_segment_num()<l.number_of_segments());
                Connectivities.push_back(component[comp_lst[it.center().get_segment_num()][0][it.center().run_type_position()]]==source_node);          //TODO
            } else {
                int k;
                for (k=0;k<2*D;++k) {
                    if (component[comp_lst[it.neighbor(k).get_segment_num()][it.neighbor(k).get_level()][it.neighbor(k).run_type_position()]]==source_node) break;
                }
                Connectivities.push_back(k!=2*D);
            }
        }

        return std::make_pair(num_components, num_components_after);

        }

        template <class LStype> void CalculateVisibilities(
        const LStype& l,
        std::vector<bool>& Visibilities,
        int open_boundary_direction,
        bool is_open_boundary_negative) {

            const int D=LStype::dimensions;

            const typename LStype::value_type max=std::numeric_limits<typename LStype::value_type>::max();

            Visibilities.resize(l.num_active_pts());

            std::vector<typename LStype::index_type> old_indices(D-1-open_boundary_direction, std::numeric_limits<typename LStype::index_type>::max());

	    unsigned int size=1;
	    for (int i=0;i<open_boundary_direction;++i) {
		assert(!l.grid().is_pos_boundary_infinite(i));
		assert(!l.grid().is_neg_boundary_infinite(i));

		size*=(l.grid().max_point_index(i)-l.grid().min_point_index(i)+1);
	    }

	    std::vector<typename LStype::value_type> min_values(size, max);

	    typename LStype::size_type id=0;

	    typename LStype::const_iterator_runs it(l,!is_open_boundary_negative);
	    while (!it.is_finished()) {

		for (int i=0;i<D-1-open_boundary_direction;++i) {
		    bool b=false;
		    if (old_indices[i]!=it.start_indices(i+open_boundary_direction+1)) {
			old_indices[i]=it.start_indices(i+open_boundary_direction+1);
			b=true;
		    }
		    if (b) min_values.assign(size,max);
		}

		unsigned int pos_begin=0;
		unsigned int pos_end=0;

		for (int i=open_boundary_direction-1;i>=0;--i) {
		    pos_begin*=(l.grid().max_point_index(i)-l.grid().min_point_index(i)+1);
		    pos_end*=(l.grid().max_point_index(i)-l.grid().min_point_index(i)+1);
		    pos_begin+=(it.start_indices(i)-l.grid().min_point_index(i));
		    pos_end+=(it.end_indices(i)-l.grid().min_point_index(i));
		}

		if (it.is_active()) {
		    Visibilities[is_open_boundary_negative?id:(l.num_active_pts()-1-id)]=(it.value()<min_values.at(pos_begin));
		    ++id;
		}

		for (unsigned int i=pos_begin; i<=pos_end;++i) min_values.at(i)=std::min(min_values.at(i), it.value());

		if (is_open_boundary_negative) {
		    it.next();
		} else {
		    it.previous();
		}
	}

            assert(id==l.num_active_pts());
    }

        namespace {

                template <class ModelType, int Dimensions> class VelocityClass {
            const ModelType& Model;
            const double* NormalVector;
            const double* Coverages;
            const double* Rates;
            const std::vector<bool>& Connectivities;
            const std::vector<bool>& Visibilities;
        public:

            VelocityClass(	const ModelType& m,
                            const double * n,
                            const double * c,
                            const double * r,
                            const std::vector<bool>& co,
                            const std::vector<bool>& vi

                            ) : Model(m), NormalVector(n), Coverages(c), Rates(r), Connectivities(co), Visibilities(vi)  {}

                        double operator()(unsigned int active_pt,int matnum) const {
                double v;
               // double s=0.;

                Model.CalculateVelocity(
                        v,
                        calc::Make3DVector<Dimensions>(NormalVector+active_pt*Dimensions),
                        Coverages+active_pt*Model.CoverageStorageSize,
                        Rates+active_pt*Model.RatesStorageSize,
                        matnum,
                        (Model.CalculateConnectivities)?Connectivities[active_pt]:true,
                        (Model.CalculateVisibilities)?Visibilities[active_pt]:true//,
                        //s
                );
                return v;
            }
        };

		template <class ModelType, int Dimensions> class VelocityClass1 {
			const ModelType& Model;
	    const double* NormalVector;
	    const double* Coverages;
	    const double* Rates;
	    const std::vector<bool>& Connectivities;
	    const std::vector<bool>& Visibilities;
			double NextTime;
			double RelTime;
	    double *Thick;
	    double *Vel;
	    bool Top;
	public:

            VelocityClass1(	const ModelType& m,
                            const double * n,
                            const double * c,
                            const double * r,
                            const std::vector<bool>& co,
                            const std::vector<bool>& vi,
                                                        double t,
                                                        double s,
                            double *th,
                            double *ve,
                                                        bool tp

                            ) : Model(m), NormalVector(n), Coverages(c), Rates(r), Connectivities(co), Visibilities(vi), NextTime(t), RelTime(s), Thick(th), Vel(ve), Top(tp)  {}

            double operator()(unsigned int active_pt,int matnum) const {
                double v=0;
                bool Mask=false;

				Model.CalculateVelocity(
						//v_top,
						//v_bottom,
						v,
						calc::Make3DVector<Dimensions>(NormalVector+active_pt*Dimensions),
						Coverages+active_pt*Model.CoverageStorageSize,
						Rates+active_pt*Model.RatesStorageSize,
						matnum,
						(Model.CalculateConnectivities)?Connectivities[active_pt]:true,
						(Model.CalculateVisibilities)?Visibilities[active_pt]:true,
						NextTime,
						RelTime,
						Thick,
						Vel,
						Top,
						Mask
				);
		return v;
	    }
	};


		template <class ModelType, int Dimensions> class VelocityClass2 {
			const ModelType& Model;
			const double* NormalVector;
			const double* Coverages;
			const double* Rates;
			const std::vector<bool>& Connectivities;
			const std::vector<bool>& Visibilities;
		public:

			VelocityClass2(	const ModelType& m,
							const double * n,
							const double * c,
							const double * r,
							const std::vector<bool>& co,
							const std::vector<bool>& vi

							) : Model(m), NormalVector(n), Coverages(c), Rates(r), /*RelTime(s),*/ Connectivities(co), Visibilities(vi)  {}

			void scalar_velocity(double & v, unsigned int active_pt,int matnum) const {

				Model.CalculateVelocity(
						v,
						calc::Make3DVector<Dimensions>(NormalVector+active_pt*Dimensions),
						Coverages+active_pt*Model.CoverageStorageSize,
						Rates+active_pt*Model.RatesStorageSize,
						matnum,
						(Model.CalculateConnectivities)?Connectivities[active_pt]:true,
						(Model.CalculateVisibilities)?Visibilities[active_pt]:true//,
				);
			}

			void vector_velocity(double* v, unsigned int active_pt, double location, int matnum) const {
				Model.CalculateVectorVelocity(
						v,
						location,
						calc::Make3DVector<Dimensions>(NormalVector+active_pt*Dimensions),
						Coverages+active_pt*Model.CoverageStorageSize,
						Rates+active_pt*Model.RatesStorageSize,
						matnum,
						(Model.CalculateConnectivities)?Connectivities[active_pt]:true,
						(Model.CalculateVisibilities)?Visibilities[active_pt]:true,
						0.0,
						false
				);
			}
		};

		template <class ModelType, int Dimensions> class VelocityClass3 {
//			const model::Oxidation& Model;
			const ModelType& Model;
			const double* NormalVector;
			const double* Coverages;
			const double* Rates;
			const std::vector<bool>& Connectivities;
			const std::vector<bool>& Visibilities;
			double NextTime;
			double RelTime;
			double *Thick;
	    double *Vel;
	    bool Top;
	    bool Mask;
		public:

			VelocityClass3(	const ModelType& m,//const model::Oxidation& m,
							const double * n,
							const double * c,
							const double * r,
							const std::vector<bool>& co,
							const std::vector<bool>& vi,
							double t,
							double s,
			    double *th,
			    double *ve,
			    bool tp,
			    bool msk

							) : Model(m), NormalVector(n), Coverages(c), Rates(r), Connectivities(co), Visibilities(vi), NextTime(t), RelTime(s), Thick(th), Vel(ve), Top(tp), Mask(msk) {

//				if (tp){
//					double velocity=*ve;
//					velocity *=(-55./45.);
//					*Vel=velocity;
//				} else {
//					Vel = ve;
//				}

			}

//			if (Mask) std::cout << "Mask!\n";
			void scalar_velocity(double & v, unsigned int active_pt, int matnum) const {
//				double v_top=0;
//				double v_bottom=0;
				Model.CalculateVelocity(
						//v_top,
						//v_bottom,
						v,
						calc::Make3DVector<Dimensions>(NormalVector+active_pt*Dimensions),
						Coverages+active_pt*Model.CoverageStorageSize,
						Rates+active_pt*Model.RatesStorageSize,
						matnum,
						(Model.CalculateConnectivities)?Connectivities[active_pt]:true,
						(Model.CalculateVisibilities)?Visibilities[active_pt]:true,
						NextTime,
						RelTime,
						Thick,
						Vel,
						Top,
						Mask
				);
				//std::cout << "Vel: " << Vel << std::endl;
//				if (Top) {
//					v = v_top;
//					std::cout << "v_top: " << v << std::endl;
//				} else {
//					v = v_bottom;
//					std::cout << "v_bottom: " << v << std::endl;
//				}
//				if (Mask) std::cout << "Vel: " << *Vel << "\n";
			}

			void vector_velocity(double* v, unsigned int active_pt, double location, int matnum) const {
				Model.CalculateVectorVelocity(
						v,
						location,
						calc::Make3DVector<Dimensions>(NormalVector+active_pt*Dimensions),
						Coverages+active_pt*Model.CoverageStorageSize,
						Rates+active_pt*Model.RatesStorageSize,
						matnum,
						(Model.CalculateConnectivities)?Connectivities[active_pt]:true,
						(Model.CalculateVisibilities)?Visibilities[active_pt]:true,
						*Vel,
						Top,
						Mask
				);
			}
		};

        template <class ModelType, int Dimensions>
        class DataAccessClass {
            const ModelType& Model;
            const double* Coverages;
            const double* Rates;
            const double* NormalVector;
            const std::vector<unsigned int>& Materials;
            const std::vector<bool>& Connectivities;
            const std::vector<bool>& Visibilities;
            bool OutputVelocities;
            bool OutputCoverages;
            bool OutputRates;
            bool OutputMaterials;
         public:

            DataAccessClass(    const ModelType& m,
                                const double * c,
                                const double * r,
                                const double * n,
                                const std::vector<unsigned int>& ma,
                                const std::vector<bool>& co,
                                const std::vector<bool>& vi,
                                bool out_v=false,
                                bool out_c=false,
                                bool out_r=false,
                                bool out_m=false
                                ) : Model(m), Coverages(c), Rates(r), NormalVector(n), Materials(ma), Connectivities(co), Visibilities(vi), OutputVelocities(out_v), OutputCoverages(out_c), OutputRates(out_r), OutputMaterials(out_m)  {}

             int number_of_series() const {
                return (1+ModelType::CoverageStorageSize+ModelType::RatesStorageSize+1);
            }

            template <class PT_ID_TYPE>
            std::string get_series_data(PT_ID_TYPE active_pt_id, int series) const {

                std::ostringstream out;

                if (series==0) {

                    double v=0.;
                   // double s=0.;

                    unsigned int mat=0;
                    bool connected=true;
                    bool visible=true;

                    if (Materials.size()>0) mat= Materials[active_pt_id];
                    if (Connectivities.size()>0) connected=Connectivities[active_pt_id];
                    if (Visibilities.size()>0) visible=Visibilities[active_pt_id];

                        Model.CalculateVelocity(
                                        v,
                                        calc::Make3DVector<Dimensions>(NormalVector+active_pt_id*Dimensions),
                                        Coverages+active_pt_id*Model.CoverageStorageSize,
                                        Rates+active_pt_id*Model.RatesStorageSize,
                                        mat,
                                        connected,
                                        visible//,
                                        //s
                                                );

                    out << static_cast<float>(v);

                } else if (series<=ModelType::CoverageStorageSize) {

                   out <<  static_cast<float>(Coverages[active_pt_id*ModelType::CoverageStorageSize+series-1]);

                } else if (series<=ModelType::CoverageStorageSize+ModelType::RatesStorageSize) {

                    out << static_cast<float>(Rates[active_pt_id*ModelType::RatesStorageSize+series-ModelType::CoverageStorageSize-1]);

                } else {

                        unsigned int mat=0;
                        if (Materials.size()>0) mat= Materials[active_pt_id];

                    out << mat;

                }
                return out.str();
            }

            std::string get_series_label(int series) const {
                if (series==0) {
                    return std::string("Velocities");
                } else if (series<=ModelType::CoverageStorageSize) {
                    std::ostringstream out;
                    out << "Coverage" << series-1;
                    return out.str();
                } else if (series<=ModelType::CoverageStorageSize+ModelType::RatesStorageSize) {
                    std::ostringstream out;
                    out << "Rate" << series-ModelType::CoverageStorageSize-1;
                    return out.str();
                } else {
                    return std::string("Material");
                }
            }

            std::string get_series_type(int series) const {
                if (series<=ModelType::CoverageStorageSize+ModelType::RatesStorageSize) {
                    return std::string("float");
                } else {
                    return std::string("int");
                }
            }

            bool get_series_output(int series) const {
               if (series==0) {
                    return OutputVelocities;
                } else if (series<=ModelType::CoverageStorageSize) {
                    return OutputCoverages;
                } else if (series<=ModelType::CoverageStorageSize+ModelType::RatesStorageSize) {
                    return OutputRates;
                } else {
                    return OutputMaterials;
                }
            }

        };


        template <class ModelType, int Dimensions>
        class DataAccessClass2 {
            const ModelType& Model;
            const double* Coverages;
            const double* Rates;
            const double* NormalVector;
            const std::vector<unsigned int>& Materials;
            const std::vector<bool>& Connectivities;
            const std::vector<bool>& Visibilities;
            bool OutputVelocities;
            bool OutputCoverages;
            bool OutputRates;
            bool OutputMaterials;
         public:

            DataAccessClass2(   const ModelType& m,// const model::Oxidation& m,
                                const double * c,
                                const double * r,
                                const double * n,
                                const std::vector<unsigned int>& ma,
                                const std::vector<bool>& co,
                                const std::vector<bool>& vi,
                                bool out_v=false,
                                bool out_c=false,
                                bool out_r=false,
                                bool out_m=false
                                ) : Model(m), Coverages(c), Rates(r), NormalVector(n), Materials(ma), Connectivities(co), Visibilities(vi), OutputVelocities(out_v), OutputCoverages(out_c), OutputRates(out_r), OutputMaterials(out_m)/*, NextTime(t), TimeStep(s), Thick(th), Vel(ve)*/ {}

             int number_of_series() const {
                 return (1+ModelType::CoverageStorageSize+ModelType::RatesStorageSize+1);
            }

            template <class PT_ID_TYPE>
            std::string get_series_data(PT_ID_TYPE active_pt_id, int series) const {

                std::ostringstream out;

                if (series==0) {

                    double v;

                    unsigned int mat=0;
                    bool connected=true;
                    bool visible=true;
                    double NextTime=0;
                    double TimeStep=0;
                    double *Thick=0;
                    double *Vel=0;

                    if (Materials.size()>0) mat= Materials[active_pt_id];
                    if (Connectivities.size()>0) connected=Connectivities[active_pt_id];
                    if (Visibilities.size()>0) visible=Visibilities[active_pt_id];

                    Model.CalculateVelocity(
                                v,
                                calc::Make3DVector<Dimensions>(NormalVector+active_pt_id*Dimensions),
                                Coverages+active_pt_id*Model.CoverageStorageSize,
                                Rates+active_pt_id*Model.RatesStorageSize,
                                mat,
                                connected,
                                visible,
                                                NextTime,
                                                TimeStep,
                                Thick,
                                Vel,
                                true
                        );

                    out << static_cast<float>(v);

                } else if (series<=ModelType::CoverageStorageSize) {

                   out <<  static_cast<float>(Coverages[active_pt_id*ModelType::CoverageStorageSize+series-1]);

                } else if (series<=ModelType::CoverageStorageSize+ModelType::RatesStorageSize) {

                    out << static_cast<float>(Rates[active_pt_id*ModelType::RatesStorageSize+series-ModelType::CoverageStorageSize-1]);

                } else {

                        unsigned int mat=0;
                        if (Materials.size()>0) mat= Materials[active_pt_id];

                    out << mat;

                }
                return out.str();
            }

            std::string get_series_label(int series) const {
                if (series==0) {
                    return std::string("Velocities");
                } else if (series<=ModelType::CoverageStorageSize) {
                    std::ostringstream out;
                    out << "Coverage" << series-1;
                    return out.str();
                } else if (series<=ModelType::CoverageStorageSize+ModelType::RatesStorageSize) {
                    std::ostringstream out;
                    out << "Rate" << series-ModelType::CoverageStorageSize-1;
                    return out.str();
                } else {
                    return std::string("Material");
                }
            }

            std::string get_series_type(int series) const {
                if (series<=ModelType::CoverageStorageSize+ModelType::RatesStorageSize) {
                    return std::string("float");
                } else {
                    return std::string("int");
                }
            }

            bool get_series_output(int series) const {
               if (series==0) {
                    return OutputVelocities;
                } else if (series<=ModelType::CoverageStorageSize) {
                    return OutputCoverages;
                } else if (series<=ModelType::CoverageStorageSize+ModelType::RatesStorageSize) {
                    return OutputRates;
                } else {
                    return OutputMaterials;
                }
            }

	};
	}

	// Normal original ExecuteProcess2 from Ertl
	template <class LevelSetsType, class ModelType, class ParameterType, class ProcessParameterType, class OutputInfoType> void ExecuteProcess2(
				LevelSetsType& LevelSets,
				const ModelType& Model,
				const ParameterType& Parameter,
				const ProcessParameterType& ProcessParameter,
				OutputInfoType & output_info
		) {

		const int D=LevelSetsType::value_type::dimensions;


	    const std::vector<double> & OutputTimes=ProcessParameter.output_times;
	    //int OutputTimesSize = OutputTimes.size();
	    //msg::print_message("OutputTimes.size(): ", OutputTimesSize);
	    std::vector<double>::const_iterator OutputTimesIter = OutputTimes.begin();//std::lower_bound(OutputTimes.begin(), OutputTimes.end(), AbsoluteTime);











		int init_cycles=ProcessParameter.StartIterationCycles;
		int rec_cycles=ProcessParameter.IterationCycles;


		geom::cells<ParameterType::Dimension> Cells;

		std::vector<double> Coverages(std::max(LevelSets.back().num_active_pts()* Model.CoverageStorageSize,1u),0.);
		std::vector<double> Rates(1,0);
		std::vector<double> NormalVectors;
		std::vector<double> DistancesToReceiver;
		std::vector<unsigned int> PointMaterials;
		std::vector<bool> Connectivities;
		std::vector<bool> Visibilities;

		//time statistics
		const std::string TimeStatFileName=Parameter.GetCompleteOutputFileName("StatisticsTimes.cvs");
		std::ofstream f;

		//unsigned int LineNumber;
		if (Parameter.PrintStatistics) {
			if(!std::ifstream(TimeStatFileName.c_str())) {

#ifdef VERBOSE
				msg::print_message("Print Header in StatisticsTimes.cvs");
#endif

				f.open(TimeStatFileName.c_str());
				f << "Time for expansion"				<<";";
				f << "Time for normal vector calc."		<<";";
				f << "Determining materials"			<<";";
				f << "Determining connectivities"       <<";";
				f << "Reduced graph num vertices"       <<";";
				f << "num componenets"                  <<";";
				f << "Time for smoothing"               <<";";
				f << "Determining visibilities"         <<";";
				f << "Setup active cells"				<<";";
				f << "Setup partition"					<<";";
				f << "Rate calculation"					<<";";
				f << "Memory Ray Tracing Data Structure"<<";";
				f << "Level set time integration"    	<<";";
				f << "Output"							<<";";
				f << "Time for Output"					<<";";
				f << "Total time step excl. Output"		<<";";
				f << "Total time step incl. Output"		<<";";         //TODO
				f << "Chosen time step"					<<";";          //TODO
				f << "Time"								<<";";           //TODO
				f << "Left Time"						<<std::endl;
				f.close();
			}
		}

		const double & ProcessTime = ProcessParameter.ProcessTime;
		double RelativeTime=0;

		//while ((OutputTimesIter!=OutputTimes.end()) && (RelativeTime>*OutputTimesIter)) ++OutputTimesIter;

#ifdef VERBOSE
				msg::print_message("Start loop over time");
#endif


                while(true) {

                    double TimeTotalExclOutput=-my::time::GetTime();
            double TimeTotalInclOutput=-my::time::GetTime();
            double TimeExpansion=0;
            double TimeNormals=0;
            double TimeMaterials=0;
            double TimeCells=0;
            double TimePartition=0;
            double TimeRates=0;
            double TimeTimeIntegration=0;
            double TimeOutput=0;
            double TimeConnectivities=0;
            double TimeVisibilities=0;
            double TimeSmoothing=0;
            double ray_tracing_memory=0;

            unsigned int graph_size=0;
            unsigned int num_components=0;

		    bool MakeOutput=false;
		    if (OutputTimesIter!=OutputTimes.end()) {
			assert(RelativeTime<=*OutputTimesIter);
			if (RelativeTime==*OutputTimesIter) {
			    MakeOutput=true;
			    OutputTimesIter++;
			}
		    }

		    //if ((RelativeTime==EndTime) && (ProcessParameter.final_output)) MakeOutput=true;
		    //if ((RelativeTime==StartTime) && (ProcessParameter.initial_output)) MakeOutput=true;
		    if (!MakeOutput) if (RelativeTime==ProcessTime) break; // Reached the end of processing time

                    //###########################
            // smooth surface level set
            //###########################

                    if (ProcessParameter.smoothing_material_level>0) {
#ifdef VERBOSE
                                msg::print_message("smoothing");
#endif

			TimeSmoothing-=my::time::GetTime();

			double time_step;

			int dummy;

			int counter=0;

                        do {
                    time_step=lvlset::time_integrate(
                                LevelSets,
                                dummy,
                                lvlset::SMOOTHING_SCHEME(ProcessParameter.smoothing_material_level, ProcessParameter.smoothing_max_curvature, ProcessParameter.smoothing_min_curvature),
                                Parameter.TimeStepRatio,
                                std::numeric_limits<double>::max(),
                                Coverages,
                                Model.CoverageStorageSize);
                                        counter++;
                        } while (time_step!=std::numeric_limits<double>::max() && counter < ProcessParameter.smoothing_max_iterations);

			if (time_step!=std::numeric_limits<double>::max()) {
				msg::print_message("maximum number of iterations reached during smoothing operation");
			}


			TimeSmoothing+=my::time::GetTime();
		    }

			/*
			//Output statistics for level sets
			if (Parameter.PrintStatistics) {
			    TimeTotalExclOutput+=my::time::GetTime();
				int i=0;
				for (typename LevelSetsType::iterator it=LevelSets.begin();it!=LevelSets.end();++it) {
					std::ostringstream tmp;
					tmp << Parameter.OutputPath << "StatisticsLevelSet" << i << ".cvs";
					lvlset::misc::PrintStatistics(*it, tmp.str());
					i++;
				}
				TimeTotalExclOutput-=my::time::GetTime();
			}
	    */

	    if (Model.ReemissionIsMaterialDependent) {
#ifdef VERBOSE
				msg::print_message("determine top most layer");
#endif
                TimeMaterials-=my::time::GetTime();
                DetermineTopMostLayer(LevelSets, PointMaterials);
                TimeMaterials+=my::time::GetTime();
            }

            if (Model.CalculateConnectivities) {
#ifdef VERBOSE
                                msg::print_message("calculate connectivities");
#endif
                TimeConnectivities-=my::time::GetTime();
                std::pair<unsigned int, unsigned int> x=CalculateConnectivities(LevelSets.back(), Connectivities, Parameter.is_open_boundary_negative);
                graph_size=x.first;
                num_components=x.second;
                TimeConnectivities+=my::time::GetTime();
            }
            if (Model.CalculateVisibilities) {
#ifdef VERBOSE
                                msg::print_message("calculate visibilities");
#endif
                TimeVisibilities-=my::time::GetTime();
                CalculateVisibilities(LevelSets.back(), Visibilities, Parameter.open_boundary_direction, Parameter.is_open_boundary_negative);
                TimeVisibilities+=my::time::GetTime();
            }

            if ((Model.CalculateNormalVectors) || (Model.NumberOfParticleTypes>0)){
#ifdef VERBOSE
                                msg::print_message("calculate normal vectors");
#endif

#ifdef VERBOSE
                                msg::print_message("expansion");
#endif
                                TimeExpansion-=my::time::GetTime();
                LevelSets.back().expand(3);
                TimeExpansion+=my::time::GetTime();
#ifdef VERBOSE
                                msg::print_message("normal vector calculation");
#endif
                TimeNormals-=my::time::GetTime();
                calc::CalculateNormalVectors(LevelSets.back(), NormalVectors, DistancesToReceiver, Parameter.open_boundary_direction, Parameter.is_open_boundary_negative, Parameter.ReceptorRadius, lvlset::vec<double,D>(Parameter.DefaultDiskOrientation));
                TimeNormals+=my::time::GetTime();
            }

            if (Model.NumberOfParticleTypes>0) {

#ifdef VERBOSE
                                msg::print_message("start monte carlo");
#endif

                std::vector<lvlset::vec<int,ParameterType::Dimension > > CellCoordinates; //Vector of vec types to set up the cell coordinates

                TimeExpansion-=my::time::GetTime();
                LevelSets.back().add_voxel_corners();
                TimeExpansion+=my::time::GetTime();

                TimeCells-=my::time::GetTime();
                calc::SetupCells(LevelSets.back(),Cells, CellCoordinates, NormalVectors, DistancesToReceiver, Parameter.ReceptorRadius);
                TimeCells+=my::time::GetTime();

                typedef typename calc::PartitionTraits<ParameterType> tmp_type;

#ifdef COMPILE_PARTITION_NEIGHBOR_LINKS_ARRAYS
                if (ProcessParameter.partition_data_structure==partition::NEIGHBOR_LINKS_ARRAYS) {
                    partition::NeighborLinksArrays<tmp_type> Partition;
                    TimePartition-=my::time::GetTime();
                    Partition.Setup(0, Cells.size(), CellCoordinates, LevelSets.back().grid().boundary_conditions(),ProcessParameter.partition_splitting_strategy,ProcessParameter.partition_surface_area_heuristic_lambda);
                    TimePartition+=my::time::GetTime();
                    ray_tracing_memory=Partition.get_memory();
                    if (Parameter.PrintStatistics) {
                        TimeTotalExclOutput+=my::time::GetTime();
                        Partition.PrintStatistics(Parameter.GetCompleteOutputFileName("StatisiticsPartition.cvs"));
                        TimeTotalExclOutput-=my::time::GetTime();
                    }

                    TimeRates-=my::time::GetTime();
                    do {
                        //[Josef] This is where the CalculateRates function found in calc.h is called.
                        calc::CalculateRates(Model,Parameter,Partition,LevelSets.back(),NormalVectors,DistancesToReceiver,Coverages,Rates,PointMaterials,Cells,RelativeTime);
                        calc::UpdateCoverages(Rates, Coverages, Model);
                        init_cycles--;
                    } while (init_cycles>=0);
                    init_cycles=rec_cycles;

                    TimeRates+=my::time::GetTime();
                }
#endif
#ifdef COMPILE_PARTITION_FULL_GRID
                if (ProcessParameter.partition_data_structure==partition::FULL_GRID) {
                        msg::print_message("FULL_GRID");
                    partition::FullGrid<tmp_type> Partition;

                    TimePartition-=my::time::GetTime();
                    Partition.Setup(0, Cells.size(), CellCoordinates, LevelSets.back().grid().boundary_conditions(),ProcessParameter.partition_splitting_strategy,ProcessParameter.partition_surface_area_heuristic_lambda);
                    TimePartition+=my::time::GetTime();
                    ray_tracing_memory=Partition.get_memory();
                    if (Parameter.PrintStatistics) {
                        TimeTotalExclOutput+=my::time::GetTime();
                        Partition.PrintStatistics(Parameter.GetCompleteOutputFileName("StatisiticsPartition.cvs"));
                        TimeTotalExclOutput-=my::time::GetTime();
                    }

                    TimeRates-=my::time::GetTime();
                    do {
                        calc::CalculateRates(Model,Parameter,Partition,LevelSets.back(),NormalVectors,DistancesToReceiver,Coverages,Rates,PointMaterials,Cells,RelativeTime);
                        calc::UpdateCoverages(Rates, Coverages, Model);
                        init_cycles--;
                    } while (init_cycles>=0);
                    init_cycles=rec_cycles;

                    TimeRates+=my::time::GetTime();
                }
#endif
#ifdef COMPILE_UP_DOWN_LINKED_TREE
                if (ProcessParameter.partition_data_structure==partition::UP_DOWN_LINKED_TREE) {
                        msg::print_message("UP_DOWN_LINKED_TREE");
                    partition::UpDownLinkTree<tmp_type> Partition;

                    TimePartition-=my::time::GetTime();
                    Partition.Setup(0, Cells.size(), CellCoordinates, LevelSets.back().grid().boundary_conditions(),ProcessParameter.partition_splitting_strategy,ProcessParameter.partition_surface_area_heuristic_lambda);
                    TimePartition+=my::time::GetTime();
                    ray_tracing_memory=Partition.get_memory();
                    if (Parameter.PrintStatistics) {
                        TimeTotalExclOutput+=my::time::GetTime();
                        Partition.PrintStatistics(Parameter.GetCompleteOutputFileName("StatisiticsPartition.cvs"));
                        TimeTotalExclOutput-=my::time::GetTime();
                    }

                    TimeRates-=my::time::GetTime();
                    do {
                        calc::CalculateRates(Model,Parameter,Partition,LevelSets.back(),NormalVectors,DistancesToReceiver,Coverages,Rates,PointMaterials,Cells,RelativeTime);
                        calc::UpdateCoverages(Rates, Coverages, Model);
                        init_cycles--;
                    } while (init_cycles>=0);
                    init_cycles=rec_cycles;

                    TimeRates+=my::time::GetTime();
                }
#endif

                        }

            //#######################################
            // output
            //#######################################

            TimeTotalExclOutput+=my::time::GetTime();
            TimeOutput-=my::time::GetTime();

            if (MakeOutput) {

#ifdef VERBOSE
                                msg::print_message("make output");
#endif


                DataAccessClass<ModelType, ParameterType::Dimension> Data(  Model,
                                                                            &Coverages[0],
                                                                            &Rates[0],
                                                                            &NormalVectors[0],
                                                                            PointMaterials,
                                                                            Connectivities,
                                                                            Visibilities,
                                                                            ProcessParameter.print_velocities || Parameter.print_velocities,
                                                                            ProcessParameter.print_coverages || Parameter.print_coverages,
                                                                            ProcessParameter.print_rates || Parameter.print_rates,
                                                                            ProcessParameter.print_materials || Parameter.print_materials
                                                                        );

                {
                    std::ostringstream oss;
                    oss<<"Write " << output_info.output_counter << "-th output (time = " << RelativeTime << ")...";
                    msg::print_start(oss.str());
                }

                typename LevelSetsType::iterator it=LevelSets.begin();
                for (unsigned int i=0;i<LevelSets.size();i++) {

                    if (Parameter.print_dx) {
                        std::ostringstream oss;
                        oss << Parameter.OutputPath<< output_info.file_name <<"_" << i << "_" << output_info.output_counter << ".dx";
#ifdef VERBOSE
                        msg::print_message("print dx");
#endif

                        if (i!=LevelSets.size()-1) {
                            write_explicit_surface_opendx(*it,oss.str());
                        } else {
                            write_explicit_surface_opendx(*it,oss.str(), Data);
                        }
                    }
                    if (Parameter.print_vtk) {
                        std::ostringstream oss;
                        oss << Parameter.OutputPath<< output_info.file_name <<"_" << i << "_" << output_info.output_counter << ".vtk";
#ifdef VERBOSE
                        msg::print_message("print vtk");
#endif

                        if (i!=LevelSets.size()-1) {
                                write_explicit_surface_vtk(*it,oss.str());
                        } else {
                            write_explicit_surface_vtk(*it,oss.str(), Data);
                        }
                    }
                    it++;
                }

                output_info.output_counter++;

                msg::print_done();
            }

            TimeOutput+=my::time::GetTime();
            TimeTotalExclOutput-=my::time::GetTime();


            bool is_finished=(RelativeTime==ProcessTime);

//            std::cout << "RelativeTime: " << RelativeTime << std::endl;
            //#######################################
            // time integration
            //#######################################
#ifdef VERBOSE
                        msg::print_message("time integration");
#endif

            double time_step=0;
            if (!is_finished) {

                //determine next time stop
                double NextTimeStop=std::min(ProcessTime, RelativeTime+ProcessParameter.MaxTimeStep);
                if (OutputTimesIter!=OutputTimes.end()) NextTimeStop=std::min(NextTimeStop, *OutputTimesIter);

                double MaxTimeStep=NextTimeStop-RelativeTime;

                if (ProcessParameter.FiniteDifferenceScheme==ProcessParameter.ENGQUIST_OSHER_1ST_ORDER) {

                        VelocityClass2<ModelType, ParameterType::Dimension> Velocities(Model, &NormalVectors[0], &Coverages[0], &Rates[0], Connectivities, Visibilities);

                        TimeExpansion-=my::time::GetTime();
                        LevelSets.back().expand(3);
                        TimeExpansion+=my::time::GetTime();

                        TimeTimeIntegration-=my::time::GetTime();
//                		if (Parameter.separate_materials) {
//                            for (typename LevelSetsType::iterator it=LevelSets.begin();it!=LevelSets.end();it++) {
//                            	time_step=lvlset::time_integrate(
//                            			*it,
//                            			Velocities,
//                            			lvlset::ENGQUIST_OSHER_SV_1ST_ORDER,
//                            			Parameter.TimeStepRatio,
//                            			MaxTimeStep,
//                            			Coverages,
//                            			Model.CoverageStorageSize);
//                            }
//                		} else {
					time_step=lvlset::time_integrate(
							LevelSets,
							Velocities,
							lvlset::ENGQUIST_OSHER_SV_1ST_ORDER,
							Parameter.TimeStepRatio,
							MaxTimeStep,
							Coverages,
							Model.CoverageStorageSize);
//                		}
			TimeTimeIntegration+=my::time::GetTime();

		} else if (ProcessParameter.FiniteDifferenceScheme==ProcessParameter.ENGQUIST_OSHER_2ND_ORDER) {

			VelocityClass2<ModelType, ParameterType::Dimension> Velocities(Model, &NormalVectors[0], &Coverages[0], &Rates[0], Connectivities, Visibilities);

                    TimeExpansion-=my::time::GetTime();
                    LevelSets.back().expand(5);
                    TimeExpansion+=my::time::GetTime();

                    TimeTimeIntegration-=my::time::GetTime();
                    time_step=lvlset::time_integrate(
                            LevelSets,
                            Velocities,
                            lvlset::ENGQUIST_OSHER_SV_2ND_ORDER,
                            Parameter.TimeStepRatio,
                            MaxTimeStep,
                            Coverages,
                            Model.CoverageStorageSize);
                    TimeTimeIntegration+=my::time::GetTime();

                } else if (ProcessParameter.FiniteDifferenceScheme==ProcessParameter.LAX_FRIEDRICHS_1ST_ORDER) {                  //TODO

                        VelocityClass<ModelType, ParameterType::Dimension> Velocities(Model, &NormalVectors[0], &Coverages[0], &Rates[0], Connectivities, Visibilities);

                    TimeExpansion-=my::time::GetTime();
                    LevelSets.back().expand(3);
                    TimeExpansion+=my::time::GetTime();

                    TimeTimeIntegration-=my::time::GetTime();
                    time_step=lvlset::time_integrate(
                            LevelSets,
                            Velocities,
                            lvlset::LAX_FRIEDRICHS_SCALAR_1ST_ORDER(ProcessParameter.LaxFriedrichsDissipationCoefficient),
                            Parameter.TimeStepRatio,
                            MaxTimeStep,
                            Coverages,
                            Model.CoverageStorageSize);
                    TimeTimeIntegration+=my::time::GetTime();

                } else assert(0);

                if (time_step>=MaxTimeStep) {
                    assert(time_step==MaxTimeStep);
                    time_step=MaxTimeStep;
                    RelativeTime=NextTimeStop;
                } else {
                    RelativeTime+=time_step;
                }


            }

            TimeTotalExclOutput+=my::time::GetTime();
            TimeTotalInclOutput+=my::time::GetTime();

            //#######################################
            // print statistics
            //#######################################
                        if (Parameter.PrintStatistics) {
#ifdef VERBOSE
                                msg::print_message("print statistics");
#endif

				f.open(TimeStatFileName.c_str(),std::ios_base::app);
				f<<TimeExpansion			<<";";
				f<<TimeNormals				<<";";
				f<<TimeMaterials			<<";";
				f<<TimeConnectivities       <<";";
				f<<graph_size               <<";";
				f<<num_components           <<";";
				f<<TimeSmoothing            <<";";
				f<<TimeVisibilities         <<";";
				f<<TimeCells				<<";";
				f<<TimePartition			<<";";
				f<<TimeRates				<<";";
				f<<ray_tracing_memory		<<";";
				f<<TimeTimeIntegration		<<";";
				f<<MakeOutput				<<";";
				f<<TimeOutput				<<";";
				f<<TimeTotalExclOutput		<<";";
				f<<TimeTotalInclOutput		<<";";
				f<<time_step				<<";";
				f<<RelativeTime             <<";";
				f<<(ProcessTime-RelativeTime)	<< std::endl;
				f.close();
			}

			if (is_finished) break;
		}
	}

}


#endif /*PROCESS_H_*/