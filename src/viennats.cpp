/* =========================================================================
   Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.

                            -----------------
                 ViennaTS - The Vienna Topography Simulator
                            -----------------

   Contact:         viennats@iue.tuwien.ac.at

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */



//!COMPILE OPTIONS#####################################
//#define TEST_MODE
//#define VERBOSE
//#define USE_HDF5

//Dimensions
#define DIMENSION_3
#define DIMENSION_2

//Processes
#define PROCESS_CONSTANT_RATES
#define PROCESS_SIMPLE_DEPOSITION
#define PROCESS_SF6_O2_PLASMA_ETCHING
#define PROCESS_SiO2_PLASMA_ETCHING
#define PROCESS_HBr_O2_PLASMA_ETCHING
#define PROCESS_NONLINEAR_DEPOSITION
#define PROCESS_WET_ETCHING
#define PROCESS_FIB

//LS Processes
#define PROCESS_PLANARIZATION
#define PROCESS_MASK
#define PROCESS_BOOLEANOPS

//Flux calculation
#define PROCESS_CALCULATEFLUX

#define COMPILE_PARTITION_NEIGHBOR_LINKS_ARRAYS
#define COMPILE_PARTITION_FULL_GRID
#define COMPILE_UP_DOWN_LINKED_TREE

#define MAX_NUM_THREADS 110

//!##################################################

#include "Time.h"

#include <cassert>
#include <sstream>
#include <vector>
#include <list>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Parameters.h"

#include "Model/ModelConstantRates.h"
#include "Model/ModelSimpleDeposition.h"
#include "Model/ModelSF6_O2PlasmaEtching.h"
#include "Model/ModelSiO2_PlasmaEtching.h"
#include "Model/ModelHBr_O2PlasmaEtching.h"
#include "Model/ModelNonlinearDeposition.h"
#include "Model/ModelWetEtching.h"
#include "Model/ModelFIB.h"

#include "Model/ModelCalculateFlux.h"

#include "Model/ModelPlanarization.h"
#include "Model/ModelMask.h"
#include "Model/ModelBooleanOps.h"

#include "Geometry.h"
#include "Process.h"
#include "LSlib/levelset.hpp"
#include "message.h"
#include "boundaries.h"

template<int D>
class GridTraitsType {
public:
	typedef int index_type;
	typedef double coord_type;
private:
	index_type Min_[D], Max_[D];
	lvlset::boundary_type BoundaryConditions_[D];
	coord_type GridDelta;
public:
	const static int dimensions = D;
	std::vector<std::string> Materials;

	template<class A, class B, class C>
	GridTraitsType(const A& min, const A& max, const B& b, C grid_delta) ://, int open_boundary) :
		GridDelta(grid_delta) {//, open_boundary_direction(open_boundary) {

		for (int i = 0; i < D; i++) {
			Min_[i] = min[i];
			Max_[i] = max[i];
			BoundaryConditions_[i] = b[i];
		}
	}

	index_type min_index(int dir) const {
		return Min_[dir];
	}

	index_type max_index(int dir) const {
		return Max_[dir];
	}

	coord_type grid_position(int dir, index_type Index) const {
		return Index * GridDelta;
	}

	lvlset::boundary_type boundary_condition(int dir) const {
		return BoundaryConditions_[dir];
	}
};

class LevelSetTraitsType {
public:
	typedef unsigned int size_type;
	typedef double value_type;
};

template<class ParameterType, int D>
class ParameterDimType: public ParameterType {
public:
	static const int Dimension = D;
	ParameterDimType(const ParameterType& p) :
		ParameterType(p) {
	}
};

class OutputInfoType {
public:

	std::string file_name;
	unsigned int output_counter;
	unsigned int process_counter;
	double start_time;
	double end_time;

	OutputInfoType() :
		file_name("Interface"), output_counter(0), process_counter(0),
				start_time(0), end_time(0) {
	}
};

template<int D, class ParameterType2>
void main_(const ParameterType2& p2) {

	const ParameterDimType<ParameterType2, D> p = p2;

	int grid_min[D]={ };
	int grid_max[D]={ };

	//!Read Geometry
	geometry::geometry<D> g;
	int num_surfaces=p.InputFiles.size();
//	geometry::surface<D> s[num_surfaces];
	geometry::surface<D> *s = new geometry::surface<D> [num_surfaces];

	if (p.surface_geometry) {
		//!Read .vtk surface geometries
		std::cout << "The geometry consists of " << p.InputFiles.size() <<" input surfaces. \n";
		for(int cs=0;cs<num_surfaces;cs++) {
			msg::print_start("Read surface input file "+p.InputFiles[cs]+"...");
			s[cs].ReadVTK(p.InputFiles[num_surfaces-cs-1], p.InputScale, p.InputTransformationDirections,
					      p.InputTransformationSigns, p.change_input_parity, p.InputShift);

			for (int h = 0; h < D; ++h) {
				grid_min[h] = std::min(grid_min[h],int(std::ceil(s[cs].Min[h] / p.GridDelta - p.snap_to_boundary_eps)));
				grid_max[h] = std::max(grid_max[h],int(std::floor(s[cs].Max[h] / p.GridDelta + p.snap_to_boundary_eps)));
			}
			msg::print_done();
#ifdef VERBOSE
	std::cout << "min = " << (s[cs].Min) << "   " << "max = " << (s[cs].Max)
			<< std::endl;
	std::cout << "min = " << (s[cs].Min / p.GridDelta) << "   " << "max = "
			<< (s[cs].Max / p.GridDelta) << std::endl;
#endif
		}
	} else {
		//!Read Geometry

		// g.Read reads a geometry file and modifies it according to the user-set parameters
		msg::print_start("Read geometry input file...");
		g.Read(p.InputFiles[0], p.InputScale, p.InputTransformationDirections,
			p.InputTransformationSigns, p.change_input_parity, p.MapMaterials,
			p.InputShift, p.IgnoreMaterials);
		{
			// output to a .vtk file the modified initial geometry
			std::ostringstream oss;
			oss << p.OutputPath << "Initial_Volume_Mesh.vtk";
			g.Write(oss.str());
		}
		for (int h = 0; h < D; ++h) {
			grid_min[h]	= std::ceil(g.Min[h] / p.GridDelta - p.snap_to_boundary_eps);
			grid_max[h] = std::floor(g.Max[h] / p.GridDelta
						+ p.snap_to_boundary_eps);
		}
#ifdef VERBOSE
	std::cout << "min = " << (g.Min) << "   " << "max = " << (g.Max)
			<< std::endl;
	std::cout << "min = " << (g.Min / p.GridDelta) << "   " << "max = "
			<< (g.Max / p.GridDelta) << std::endl;
#endif
		msg::print_done();
	}

	//determine boundary conditions for level set domain
	lvlset::boundary_type bnc[D];
	for (int hh = 0; hh < D; ++hh) {
		if ((p.boundary_conditions[hh].min == bnc::PERIODIC_BOUNDARY)
				&& (p.boundary_conditions[hh].max == bnc::PERIODIC_BOUNDARY)) {
			bnc[hh] = lvlset::PERIODIC_BOUNDARY;
		} else if ((p.boundary_conditions[hh].min == bnc::INFINITE_BOUNDARY)
				&& (p.boundary_conditions[hh].max == bnc::INFINITE_BOUNDARY)) {
			bnc[hh] = lvlset::INFINITE_BOUNDARY;
		} else if (p.boundary_conditions[hh].min == bnc::INFINITE_BOUNDARY) {
			bnc[hh] = lvlset::NEG_INFINITE_BOUNDARY;
		} else if (p.boundary_conditions[hh].max == bnc::INFINITE_BOUNDARY) {
			bnc[hh] = lvlset::POS_INFINITE_BOUNDARY;
		} else {
			bnc[hh] = lvlset::SYMMETRIC_BOUNDARY;
		}
	}

	//level set grid

#ifdef VERBOSE
		std::cout << "dim " << h << " min=" << grid_min[h] << " max="
				<< grid_max[h] << std::endl;
#endif
//	}

	GridTraitsType<D> GridProperties(grid_min, grid_max, bnc, p.GridDelta);

	lvlset::grid_type<GridTraitsType<D> > grid(GridProperties);

	//Transform to surfaces
	msg::print_start("Extract surface and interfaces...");
	typedef std::list<geometry::surface<D> > SurfacesType;
	SurfacesType Surfaces;
	{
		std::bitset<2 * D> remove_flags;

		for (int i = 0; i < D; ++i) {
			if (p.boundary_conditions[i].min == bnc::PERIODIC_BOUNDARY
					|| p.boundary_conditions[i].min == bnc::REFLECTIVE_BOUNDARY
					|| p.boundary_conditions[i].min == bnc::EXTENDED_BOUNDARY) {
				remove_flags.set(i);
			} else {
				if (i == p.open_boundary_direction
						&& !p.is_open_boundary_negative && p.remove_bottom)
					remove_flags.set(i);
			}
			if (p.boundary_conditions[i].max == bnc::PERIODIC_BOUNDARY
					|| p.boundary_conditions[i].max == bnc::REFLECTIVE_BOUNDARY
					|| p.boundary_conditions[i].max == bnc::EXTENDED_BOUNDARY) {
				remove_flags.set(i + D);
			} else {
				if (i == p.open_boundary_direction
						&& p.is_open_boundary_negative && p.remove_bottom)
					remove_flags.set(i + D);
			}
		}

		if (p.surface_geometry) {
			for(int cs=num_surfaces-1;cs>=0;cs--) Surfaces.push_back(s[cs]);
		} else {
			geometry::TransformGeometryToSurfaces(g, Surfaces, remove_flags,
												  p.GridDelta * p.snap_to_boundary_eps);
		}
	}
	msg::print_done();

	//Output of initial surfaces
	int SurfaceCounter = 0;
	for (typename SurfacesType::const_iterator it = Surfaces.begin(); it != Surfaces.end(); ++it) {//unsigned int i=0;i<Surfaces.size();++i) {
		std::ostringstream oss, oss2;
		oss << p.OutputPath << "Interface" << "Initial" << "_"
				<< SurfaceCounter << ".dx";
		oss2 << p.OutputPath << "Interface" << "Initial" << "_"
				<< SurfaceCounter << ".vtk";
		it->Write(oss.str());
		it->WriteVTK(oss2.str());
		++SurfaceCounter;
	}
	//Create levelsets

	typedef std::list<lvlset::levelset<GridTraitsType<D> , LevelSetTraitsType> > LevelSetsType;
	LevelSetsType LevelSets; //list of all level set functions

	msg::print_start("Distance transformation...");

	for (typename SurfacesType::const_iterator it = Surfaces.begin(); it != Surfaces.end(); ++it) {
		LevelSets.push_back(lvlset::levelset<GridTraitsType<D> , LevelSetTraitsType>(grid));
		lvlset::init(LevelSets.back(), *it);
	}

	msg::print_done();

	msg::print_start("Add Initial Layers...");
	proc::AddLayer(LevelSets, p.AddLayer);
	msg::print_done();

	//organization of the output information by initiation of required models
	OutputInfoType output_info;

	for (typename std::list<typename ParameterType2::ProcessParameterType>::const_iterator
			pIter = p.ProcessParameters.begin(); pIter
			!= p.ProcessParameters.end(); ++pIter) {
		{
			std::ostringstream oss;
			oss << "Start execution of process \"" + pIter->ModelName + "\""
					<< std::endl << "(processing time = " << pIter->ProcessTime
					<< ")";
			msg::print_message(oss.str());
		}
		output_info.end_time += pIter->ProcessTime;
		std::cout << "AddLayer = " << pIter->AddLayer << "\n";
		proc::AddLayer(LevelSets, pIter->AddLayer);

#ifdef PROCESS_CONSTANT_RATES
		if (pIter->ModelName == "ConstantRates") {
			model::ConstantRates m(pIter->ModelParameters, pIter->MaskLayer);
			proc::ExecuteProcess(LevelSets, m, p, *pIter, output_info);
		}
#endif

#ifdef PROCESS_SIMPLE_DEPOSITION
		if (pIter->ModelName == "SimpleDeposition") {
			model::SimpleDeposition m(pIter->ModelParameters);
			proc::ExecuteProcess(LevelSets, m, p, *pIter, output_info);
		}
#endif

#ifdef PROCESS_SF6_O2_PLASMA_ETCHING
		if (pIter->ModelName == "SF6_O2PlasmaEtching") {
			model::SF6_O2PlasmaEtching m(pIter->ModelParameters);
			proc::ExecuteProcess(LevelSets, m, p, *pIter, output_info);
		}
#endif

#ifdef PROCESS_SiO2_PLASMA_ETCHING
		if (pIter->ModelName == "SiO2_PlasmaEtching") {
			model::SiO2_PlasmaEtching m(pIter->ModelParameters);
			proc::ExecuteProcess(LevelSets, m, p, *pIter, output_info);
		}
#endif

#ifdef PROCESS_HBr_O2_PLASMA_ETCHING
		if (pIter->ModelName == "HBr_O2PlasmaEtching") {
			model::HBr_O2PlasmaEtching m(pIter->ModelParameters);
			proc::ExecuteProcess(LevelSets, m, p, *pIter, output_info);
		}
#endif

#ifdef PROCESS_NONLINEAR_DEPOSITION
		if (pIter->ModelName == "NonlinearDeposition") {
			model::NonlinearDeposition m(pIter->ModelParameters);
			proc::ExecuteProcess(LevelSets, m, p, *pIter, output_info);
		}
#endif

#ifdef PROCESS_WET_ETCHING
		if (pIter->ModelName == "WetEtching") {
			model::WetEtching m(pIter->ModelParameters);
			proc::ExecuteProcess(LevelSets, m, p, *pIter, output_info);
		}
#endif

#ifdef PROCESS_FIB
		if (pIter->ModelName == "FIB") {
			model::FIB m(pIter->ModelParameters);
			proc::ExecuteProcess(LevelSets, m, p, *pIter, output_info);
		}
#endif

#ifdef PROCESS_CALCULATEFLUX
		if (pIter->ModelName == "CalculateFlux") {
			model::CalculateFlux m(pIter->ModelParameters);
			proc::ExecuteProcess(LevelSets, m, p, *pIter, output_info);
		}
#endif

#ifdef PROCESS_PLANARIZATION
    if (pIter->ModelName=="Planarization") {
        model::Planarization m(pIter->ModelParameters);
        proc::ExecuteProcess(LevelSets, m, p, *pIter,output_info);
    }
#endif

#ifdef PROCESS_MASK
    if (pIter->ModelName=="Mask") {
        model::Mask m(pIter->ModelParameters);
        proc::ExecuteProcess(LevelSets, m, p, *pIter,output_info);
    }
#endif

#ifdef PROCESS_BOOLEANOPS
    if (pIter->ModelName=="BooleanOperation") {
        model::BooleanOps m(pIter->ModelParameters);
        proc::ExecuteProcess(LevelSets, m, p, *pIter,output_info);
    }
#endif

		output_info.start_time = output_info.end_time;
		output_info.process_counter++;
	}
	delete [] s;
}

int main(int argc, char *argv[]) {

  double timer = my::time::GetTime();

	msg::print_welcome();

	//check intrinsinc double-type
	assert(std::numeric_limits<double>::is_iec559);

	//!Read Parameters-File
	par::Parameters p(argv[1]);

	//!set maximum number of threads
#ifdef _OPENMP
	if (p.OpenMP_threads>0) omp_set_num_threads(p.OpenMP_threads);
#endif

	//!Initialize Random-Generators
	int my_rank = 0;
	int num_nodes = 1;
	my::stat::InitRandomGenerator(my_rank, num_nodes, p.RNG_Seed, p.RNG_Type,
			p.RNG_Par);

#ifdef DIMENSION_2
	if (p.Dimensions == 2)
		main_<2, par::Parameters> (p);
#endif

#ifdef DIMENSION_3
	if (p.Dimensions == 3)
		main_<3, par::Parameters> (p);
#endif

	//!Finalize Random-Generators
	my::stat::FreeRandomGenerator();

  double exec_time = my::time::GetTime()-timer;
  std::stringstream ss;
  ss << exec_time;
	msg::print_message("Finished - exec-time: "+ss.str()+" s");

	return 0;

}

