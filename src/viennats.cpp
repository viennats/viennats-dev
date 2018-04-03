/* =========================================================================
   Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.

                            -----------------
                 ViennaTS - The Vienna Topography Simulator
                            -----------------

   Contact:         viennats@iue.tuwien.ac.at

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */



//COMPILE OPTIONS#####################################
#define VIENNATS_VERSION 0.11
#define TEST_MODE
//#define VERBOSE

//Dimensions
#define DIMENSION_3
#define DIMENSION_2

//Processes
#define PROCESS_CONSTANT_RATES
#define PROCESS_SIMPLE_DEPOSITION
#define PROCESS_TiN_ALD
#define PROCESS_TiN_PEALD
#define PROCESS_TiO2_ALD
#define PROCESS_SF6_O2_PLASMA_ETCHING
#define PROCESS_Cl2_CH4_PLASMA_ETCHING
#define PROCESS_BCl3_PLASMA_ETCHING
#define PROCESS_SiO2_PLASMA_ETCHING
#define PROCESS_SF6_CH2F2_PLASMA_ETCHING
#define PROCESS_Cl2_N2_ETCHING
#define PROCESS_CFx_DEPOSITION
#define PROCESS_HfO2_DEPOSITION
#define PROCESS_HBr_O2_PLASMA_ETCHING
#define PROCESS_N2_PLASMA_ETCHING
#define PROCESS_NONLINEAR_DEPOSITION
#define PROCESS_TWOSPECIES_DEPOSITION
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

//##################################################

#include "Time.h"

#include <cassert>
#include <sstream>
#include <vector>
#include <list>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Parameters.h"

#ifdef PROCESS_CONSTANT_RATES
#include "Model/ModelConstantRates.h"
#endif
#ifdef PROCESS_SIMPLE_DEPOSITION
#include "Model/ModelSimpleDeposition.h"
#endif
#ifdef PROCESS_TiN_ALD
#include "Model/ModelTiNAtomicLayerDeposition.h"
#endif
#ifdef PROCESS_TiN_PEALD
#include "Model/ModelTiNPlasmaEnhancedAtomicLayerDeposition.h"
#endif
#ifdef PROCESS_TiO2_ALD
#include "Model/ModelTiO2AtomicLayerDeposition.h"
#endif
#ifdef PROCESS_SF6_O2_PLASMA_ETCHING
#include "Model/ModelSF6_O2PlasmaEtching.h"
#endif
#ifdef PROCESS_Cl2_CH4_PLASMA_ETCHING
#include "Model/ModelCl2_CH4PlasmaEtching.h"
#endif
#ifdef PROCESS_BCl3_PLASMA_ETCHING
#include "Model/ModelBCl3PlasmaEtching.h"
#endif
#ifdef PROCESS_SiO2_PLASMA_ETCHING
#include "Model/ModelSiO2_PlasmaEtching.h"
#endif
#ifdef PROCESS_SF6_CH2F2_PLASMA_ETCHING
#include "Model/ModelSF6_CH2F2PlasmaEtching.h"
#endif
#ifdef PROCESS_Cl2_N2_ETCHING
#include "Model/ModelCl2_N2Etching.h"
#endif
#ifdef PROCESS_N2_PLASMA_ETCHING
#include "Model/ModelN2_Flash.h"
#endif
#ifdef PROCESS_HBr_O2_PLASMA_ETCHING
#include "Model/ModelHBr_O2PlasmaEtching.h"
#endif
#ifdef PROCESS_CFx_DEPOSITION
#include "Model/ModelCFx_Deposition.h"
#endif
#ifdef PROCESS_HfO2_DEPOSITION
#include "Model/ModelHfO2_Deposition.h"
#endif
#ifdef PROCESS_NONLINEAR_DEPOSITION
#include "Model/ModelNonlinearDeposition.h"
#endif
#ifdef PROCESS_TWOSPECIES_DEPOSITION
#include "Model/ModelTwoSpeciesDeposition.h"
#endif
#ifdef PROCESS_WET_ETCHING
#include "Model/ModelWetEtching.h"
#endif
#ifdef PROCESS_FIB
#include "Model/ModelFIB.h"
#endif

#ifdef PROCESS_CALCULATEFLUX
#include "Model/ModelCalculateFlux.h"
#endif

#include "Model/ModelPlanarization.h"
#include "Model/ModelMask.h"
#include "Model/ModelBooleanOps.h"

#include "Geometry.h"
#include "Process.h"
#include "LSlib/levelset.hpp"
#include "message.h"
#include "boundaries.h"


///GridTraitsType Contains minimum and maximum indexes in each dimensional direction, boundary conditions, and the grid spacing
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
  ///Can be 2 or 3-dimensional
  const static int dimensions = D;
  const index_type* minima() const {return Min_;}
  const index_type* maxima() const {return Max_;}
  const lvlset::boundary_type* boundary_conditions() const {return BoundaryConditions_;}


  ///GridTraitsType constructor to populate the minimum and maximum values and the boundary conditions
  template<class A, class B, class C>
  GridTraitsType(const A& min, const A& max, const B& b, C grid_delta) ://, int open_boundary) :
    GridDelta(grid_delta) {//, open_boundary_direction(open_boundary) {
    for (int i = 0; i < D; i++) {
      Min_[i] = min[i];
      Max_[i] = max[i];
      BoundaryConditions_[i] = b[i];
    }
  }
  //copy constructor
  GridTraitsType(const GridTraitsType& g) :
    GridDelta(g.GridDelta) {
    for (int i = 0; i < D; i++) {
      Min_[i] = g.Min_[i];
      Max_[i] = g.Max_[i];
      BoundaryConditions_[i] = g.BoundaryConditions_[i];
    }
  }
  //empty constructor
  GridTraitsType(){}

  void print() const{
    std::ostringstream oss;
    for(int i=D;i--;){
      oss << "Dimension " << i << ":" << std::endl
          << "    Grid min: " << Min_[i] << std::endl
          << "    Grid max: " << Max_[i] << std::endl
          << "    Boundary condition: " << (int)BoundaryConditions_[i];
    }
    oss << "Grid delta: " << GridDelta;
    msg::print_message(oss.str());
  }

  ///Returns the minimum index at in the "dir" axial direction
  index_type min_index(int dir) const {
    return Min_[dir];
  }

  ///Returns the maximum index at in the "dir" axial direction
  index_type max_index(int dir) const {
    return Max_[dir];
  }

  ///Returns the grid position in the "dir" axial direction (Index*GridDelta)
  coord_type grid_position(int dir, index_type Index) const {
    return Index * GridDelta;
  }

  ///Returns the boundary conditions in the "dir" axial direction
  lvlset::boundary_type boundary_condition(int dir) const {
    return BoundaryConditions_[dir];
  }
};

///Defines the size_type (unsigned int) and value_type (double) for the level set function
class LevelSetTraitsType {
public:
  typedef unsigned int size_type;
  typedef double value_type;
};

///ParameterType is the initial parameter type, which is inherited by a dimension-dependent ParameterDimType
template<class ParameterType, int D>
class ParameterDimType:
  public ParameterType {
    public:
    static const int Dimension = D;
    ParameterDimType(const ParameterType& p) :
      ParameterType(p) {}
};

///OutputInfoType stores the information necessary to output a geometry
class OutputInfoType {
public:

  ///Name of the output file
  std::string file_name;
  ///Counter keeping track of the outputs
  unsigned int output_counter;
  ///Counter keeping track of the processes
  unsigned int process_counter;
  ///Start simulation time
  double start_time;
  ///End simulation time
  double end_time;

  ///OutputInfoType constructor - only set file name to begin with "Interface"
  OutputInfoType() :
    file_name("Interface"), output_counter(0), process_counter(0),
        start_time(0), end_time(0) {
  }
};

/// Dimension specific main function reading input and starting correct process
template<int D, class ParameterType2>
void main_(ParameterType2& p2) {          //TODO changed from const to not const

  ParameterDimType<ParameterType2, D> p = p2;    //TODO changed to not const

  GridTraitsType<D> GridProperties;
  lvlset::grid_type<GridTraitsType<D> > grid;

  //Create the LevelSets - a list of all level set functions
  typedef lvlset::levelset<GridTraitsType<D> , LevelSetTraitsType> LevelSetType;
  typedef std::list<LevelSetType> LevelSetsType;
  LevelSetsType LevelSets; //list of all level set functions

  //if the first input is a lvst file, then it is assumed all are a lvst files
  if(p.geometry_files[0].find(".lvst") != std::string::npos){
    //import levelsets
    GridProperties = lvlset::get_grid_traits_from_lvst_file<GridTraitsType<D>>(p.geometry_files[0]);
    grid = lvlset::grid_type<GridTraitsType<D>>(GridProperties);
    for(unsigned int i=0; i<p.geometry_files.size(); i++){
      msg::print_start("Read levelset input file " + p.geometry_files[i] + "...");
      LevelSets.push_back(LevelSetType(grid));
      LevelSets.back().import_levelset(p.geometry_files[i]);
      msg::print_done();
    }
  }
  else {
    //read in surfaces or volume and transform them to levelsets
    if(p.surface_geometry){
      for(unsigned int i=0; i<p.geometry_files.size(); i++)
        geometry::import_levelset_from_surface<D, GridTraitsType<D>, ParameterType2, LevelSetType>(GridProperties, grid, p, LevelSets, i);
    } else {
      geometry::import_levelsets_from_volume<D, GridTraitsType<D>, ParameterType2, LevelSetType>(GridProperties, grid, p, LevelSets);
    }
  }

  //output initial LevelSets
  {
    int LevelsetCounter = 0;
    for(auto ls: LevelSets){
      std::ostringstream oss;
      oss << p.output_path << "InterfaceInitial" << LevelsetCounter++ << ".lvst";
      ls.export_levelset(oss.str(), p.bits_per_distance);
    }
  }

  if(p.add_layer>0){
    std::stringstream oss;
    oss << "Adding " << p.add_layer << " initial layer" << ((p.add_layer==1)?"...":"s...");
    msg::print_start(oss.str());
    proc::AddLayer(LevelSets, p.add_layer);
    msg::print_done();
  }

  //organization of the output information by initiation of required models
  OutputInfoType output_info;

  //!Initialize the required models and call "proc::ExecuteProcess(...)"
  //!    Possible models are: ConstantRates, SimpleDeposition, SF6_O2PlasmaEtching, SiO2_PlasmaEtching,
  //!    HBr_O2PlasmaEtching, NonlinearDeposition, WetEtching, FIB, CalculateFlux, Planarization, Mask,
  //!    and BooleanOperation

#ifdef PROCESS_TiO2_ALD
            std::vector<double> CoveragesALD_TiO2(4*LevelSets.back().num_active_pts(),0.);
#endif
#ifdef PROCESS_TiN_ALD
            std::vector<double> CoveragesALD_TiN(12*LevelSets.back().num_active_pts(),0.);
            for (unsigned int i=0;i<CoveragesALD_TiN.size();i++) CoveragesALD_TiN[i]=(i%12==10)?1.:0.;
#endif
#ifdef PROCESS_TiN_PEALD
            std::vector<double> CoveragesPEALD_TiN(12*LevelSets.back().num_active_pts(),0.);
            for (unsigned int i=0;i<CoveragesPEALD_TiN.size();i++) CoveragesPEALD_TiN[i]=(i%12==10)?1.:0.;
#endif

  for (typename std::list<typename ParameterType2::ProcessParameterType>::iterator
      pIter = p.process_parameters.begin(); pIter
      != p.process_parameters.end(); ++pIter) {
    {
      std::ostringstream oss;
      oss << "Start execution of process \"" + pIter->ModelName + "\""
          << std::endl << "(processing time = " << pIter->ProcessTime
          << ")";
      msg::print_message(oss.str());
    }
    output_info.end_time += pIter->ProcessTime;

    //Reassign Active layers to correspond to kernel layer numbering
    if(pIter->ActiveLayers.size()>LevelSets.size()) assert(0);

    LevelSetsType temp_levelSets;

    std::vector<int> layer_order;
    for(unsigned int i=0; i<pIter->ActiveLayers.size(); i++){
      layer_order.push_back(pIter->ActiveLayers[i]-1);
      pIter->ActiveLayers[i] = LevelSets.size() - pIter->ActiveLayers[i];  //reorder for model use
    }
    //if(pIter->ActiveLayers.empty()) temp_levelSets = LevelSets;  //if no materials specified, etch highest one

    //put inactive layers to new levelset ordering
    typename LevelSetsType::iterator LSIter = LevelSets.begin(), LSIter_old;

    std::cout << "Inactive/Mask/Active: ";
    for(unsigned int i=0; i<LevelSets.size(); ++i){
      if(!my::stat::AnyElement<int>(layer_order, i) && (pIter->MaskLayers.empty() || pIter->MaskLayers[0] != int(i+1))){  //neither mask nor active
        std::cout << i << ",";
        temp_levelSets.push_back(*LSIter);
      }
      ++LSIter;
    }
    std::cout << '\b' << " ";

    //add mask layers on top of inactive
    if(!pIter->MaskLayers.empty()){
      for(unsigned int i=0; i<pIter->MaskLayers.size(); ++i){
        pIter->MaskLayers[i] -= 1;  //kernel numbering
        assert(unsigned(pIter->MaskLayers[i]) < LevelSets.size());  //check if numbering is correct
      }
      LSIter = LevelSets.begin();
      for(int a=0; a<pIter->MaskLayers[0]; a++)  LSIter++;  // advance iterator to first mask layer
      temp_levelSets.push_back(*LSIter);     //this is now the only mask layer, all the other ones are AND'ed onto it
      std::cout << "\b/" << pIter->MaskLayers[0];
      for(unsigned int i=1; i<pIter->MaskLayers.size(); i++){
        std::cout << "," << pIter->MaskLayers[i];
        LSIter = LevelSets.begin();
        for(int a=0; a<pIter->MaskLayers[i]; a++)  LSIter++;    //Advance iterator to corresponding levelset
        temp_levelSets.back().min(*LSIter);      // Union second mask levelset with first
        temp_levelSets.back().prune();    //remove unnecessary points
      }
      temp_levelSets.back().segment();    //segment levelset to balance load
    }

    if(!layer_order.empty()){
      std::cout << "/";
      for(unsigned int i=0; i<layer_order.size(); i++){    //Reorder active Levelsets for next step
        std::cout << layer_order[i] << ",";
        assert(unsigned(layer_order[i]) < LevelSets.size());
        LSIter = LevelSets.begin();
        for(int a=0; a<layer_order[i]; ++a)  LSIter++;  //Advance iterator to corresponding levelset
        temp_levelSets.push_back(*LSIter);      //push levelset to temporary list
      }
      std::cout << '\b' << " ";
    }

    //wrap new top levelset around lower layers
    LSIter = temp_levelSets.begin();  //Advance iterator to first reassigned LS
    for(unsigned i=0; i<temp_levelSets.size()-1; ++i){
      temp_levelSets.back().min(*LSIter);
      temp_levelSets.back().prune();
      ++LSIter;
    }
    temp_levelSets.back().segment();

    std::swap(LevelSets, temp_levelSets);

    std::cout << std::endl;
    if(pIter->AddLayer>0){
      std::cout << "Add Layer = " << pIter->AddLayer << "\n";
      proc::AddLayer(LevelSets, pIter->AddLayer);
    }

    for(int i=0; i<pIter->AddLayer; i++) pIter->ActiveLayers.push_back(i+1);
    std::cout << "Active/Total Layers: " << pIter->ActiveLayers.size() << "/" << LevelSets.size() << "\n\n";

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

#ifdef PROCESS_Cl2_N2_ETCHING
    if(pIter->ModelName == "Cl2_N2Etching"){
      model::Cl2_N2Etching<typename ParameterType2::ProcessParameterType> m(pIter);
      proc::ExecuteProcess(LevelSets, m, p,*pIter, output_info);
    }
#endif

#ifdef PROCESS_Cl2_CH4_PLASMA_ETCHING
    if (pIter->ModelName == "Cl2_CH4PlasmaEtching") {
      model::Cl2_CH4PlasmaEtching<typename ParameterType2::ProcessParameterType> m(pIter);
      proc::ExecuteProcess(LevelSets, m, p, *pIter, output_info);
    }
#endif

#ifdef PROCESS_BCl3_PLASMA_ETCHING
    if (pIter->ModelName == "BCl3PlasmaEtching") {
      model::BCl3PlasmaEtching m(pIter->ModelParameters);
      proc::ExecuteProcess(LevelSets, m, p, *pIter, output_info);
    }
#endif

#ifdef PROCESS_SiO2_PLASMA_ETCHING
    if (pIter->ModelName == "SiO2_PlasmaEtching") {
      model::SiO2_PlasmaEtching m(pIter->ModelParameters);
      proc::ExecuteProcess(LevelSets, m, p, *pIter, output_info);
    }
#endif

#ifdef PROCESS_SF6_CH2F2_PLASMA_ETCHING
    if (pIter->ModelName == "SF6_CH2F2PlasmaEtching") {
      model::SF6_CH2F2_PlasmaEtching<typename ParameterType2::ProcessParameterType> m(pIter);
      proc::ExecuteProcess(LevelSets, m, p, *pIter, output_info);
    }
#endif

#ifdef PROCESS_N2_PLASMA_ETCHING
    if(pIter->ModelName == "N2_Flash"){
      model::N2_FLASH m(pIter->ModelParameters);
      proc::ExecuteProcess(LevelSets, m, p, *pIter, output_info);
    }
#endif

#ifdef PROCESS_HBr_O2_PLASMA_ETCHING
    if (pIter->ModelName == "HBr_O2PlasmaEtching") {
      model::HBr_O2_PlasmaEtching<typename ParameterType2::ProcessParameterType> m(pIter);
      proc::ExecuteProcess(LevelSets, m, p, *pIter, output_info);
    }
#endif

#ifdef PROCESS_CFx_DEPOSITION
    if (pIter->ModelName == "CFx_Deposition") {
      model::CFx_Deposition m(pIter->ModelParameters);
      proc::ExecuteProcess(LevelSets, m, p, *pIter, output_info);
    }
#endif

#ifdef PROCESS_HfO2_DEPOSITION
    if (pIter->ModelName == "HfO2_Deposition") {
      model::HfO2_Deposition m(pIter->ModelParameters);
      proc::ExecuteProcess(LevelSets, m, p, *pIter, output_info);
    }
#endif

#ifdef PROCESS_NONLINEAR_DEPOSITION
    if (pIter->ModelName == "NonlinearDeposition") {
      model::NonlinearDeposition m(pIter->ModelParameters);
      proc::ExecuteProcess(LevelSets, m, p, *pIter, output_info);
    }
#endif

#ifdef PROCESS_TWOSPECIES_DEPOSITION
    if (pIter->ModelName == "TwoSpeciesDeposition") {
      model::TwoSpeciesDeposition m(pIter->ModelParameters);
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
      model::CalculateFlux m(pIter->ModelParameters, D);
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

#ifdef PROCESS_TiN_ALD
    if (pIter->ModelName == "TiN_ALD") {
                    for (unsigned int i=0;i<CoveragesALD_TiN.size();i++)
                        if ((i%12==6*(pIter->ALDStep-1))||(i%12==6*(pIter->ALDStep-1)+1)||(i%12==5)||(i%12==11))
                            CoveragesALD_TiN[i]=0.;
                    model::TiN_ALD m(pIter->ModelParameters, pIter->ALDStep);
                    proc::ExecuteProcess(LevelSets, m, p, *pIter, output_info, CoveragesALD_TiN);
        }
#endif

#ifdef PROCESS_TiN_PEALD
    if (pIter->ModelName == "TiN_PEALD") {
                    for (unsigned int i=0;i<CoveragesPEALD_TiN.size();i++)
                        if ((i%12==6*(pIter->ALDStep-1))||(i%12==6*(pIter->ALDStep-1)+1)||(i%12==5)||(i%12==11))
                            CoveragesPEALD_TiN[i]=0.;
                    model::TiN_PEALD m(pIter->ModelParameters, pIter->ALDStep);
                    proc::ExecuteProcess(LevelSets, m, p, *pIter, output_info, CoveragesPEALD_TiN);
          }
#endif

#ifdef PROCESS_TiO2_ALD
    if (pIter->ModelName == "TiO2_ALD") {
      model::TiO2_ALD m(pIter->ModelParameters);
      proc::ExecuteProcess(LevelSets, m, p, *pIter, output_info, CoveragesALD_TiO2);
          }
#endif

    output_info.start_time = output_info.end_time;
    output_info.process_counter++;
  }
}

/**
*  Main function: Times the whole program and passes parameters to main_()
*/

int main(int argc, char *argv[]) {
  if(argc > 1){
    if(argv[1][0] == '-'){//if option was passed
      std::string option(argv[1]);
      std::ostringstream oss;
      if(option == "-ls2vtk" || option == "-ls2dx"){
        //assume all remaining argvs are lvst files, convert them to vtk or dx
        for(int i=2; i<argc; i++){
          std::string file(argv[i]);
          char buff[8] = {};
          std::ifstream fin(file);
          fin.read(buff, 8);//read file header
          fin.close();
          if(std::string("LvSt").compare(std::string(buff).substr(0,4))) msg::print_error(file + " is not a lvst file.");
          const int D = buff[6]-48;
          oss.str("");
          oss << file.substr(0, file.find(".lvst")) << (option == "-ls2vtk" ? ".vtk" : ".dx");

          msg::print_start("Converting " + file + std::string(" to a ") + std::string((option == "-ls2vtk" ? "vtk" : "dx")) + std::string(" file....")  );
          if(D == 2){
            GridTraitsType<2> gridP = lvlset::get_grid_traits_from_lvst_file<GridTraitsType<2>>(file);
            lvlset::grid_type<GridTraitsType<2>> grid(gridP);
            lvlset::levelset<GridTraitsType<2>> ls(grid);
            ls.import_levelset(file);
            if(option == "-ls2vtk") write_explicit_surface_vtk(ls, oss.str());
            else if(option == "-ls2dx") write_explicit_surface_opendx(ls, oss.str());
          }
          else if(D == 3){
            GridTraitsType<3> gridP = lvlset::get_grid_traits_from_lvst_file<GridTraitsType<3>>(file);
            lvlset::grid_type<GridTraitsType<3>> grid(gridP);
            lvlset::levelset<GridTraitsType<3>> ls(grid);
            ls.import_levelset(file);
            if(option == "-ls2vtk") write_explicit_surface_vtk(ls, oss.str());
            else if(option == "-ls2dx") write_explicit_surface_opendx(ls, oss.str());
          }
          msg::print_done();
        }
      }
      else if(option == "-p" || option == "-p2o"){
        for(int i=2; i<argc; i++){
          std::string file(argv[i]);
          char buff[8] = {};
          std::ifstream fin(file);
          fin.read(buff, 8);//read file header
          fin.close();
          if(std::string("LvSt").compare(std::string(buff).substr(0,4))) msg::print_error(file + " is not a lvst file.");
          const int D = buff[6]-48;
          std::ofstream fout;
          if(option == "-p2o"){
            oss.str("");
            oss << file.substr(0, file.find(".lvst")) << ".txt";
            fout = std::ofstream(oss.str());
            msg::print_start("Converting " + file + " to a txt file...");
          }

          if(D == 2){
            GridTraitsType<2> gridP = lvlset::get_grid_traits_from_lvst_file<GridTraitsType<2>>(file);
            lvlset::grid_type<GridTraitsType<2>> grid(gridP);
            lvlset::levelset<GridTraitsType<2>> ls(grid);
            ls.import_levelset(file);
            if(option == "-p2o") ls.print_without_segmentation(fout);
            else if(option == "-p") ls.print_without_segmentation();
          }
          else if(D == 3){
            GridTraitsType<3> gridP = lvlset::get_grid_traits_from_lvst_file<GridTraitsType<3>>(file);
            lvlset::grid_type<GridTraitsType<3>> grid(gridP);
            lvlset::levelset<GridTraitsType<3>> ls(grid);
            ls.import_levelset(file);
            if(option == "-p2o") ls.print_without_segmentation(fout);
            else if(option == "-p") ls.print_without_segmentation();
          }
          if(option == "-p2o"){
            fout.close();
            msg::print_done();
          }
        }
      }
      else if(option == "-h"){
        msg::print_help();
      }
      else if(option == "--help"){
        msg::print_help_extended();
      }
      else if(option == "--version"){
        msg::print_version();
      }
      else {
        oss << argv[1] << " is not an available option.\nSee --help for more information.";
        msg::print_message_2(oss.str());
      }
    }
    else{
      double timer = my::time::GetTime();

      msg::print_welcome();

      //check intrinsic double-type
      assert(std::numeric_limits<double>::is_iec559);

      //!Read Parameters-File and populate Parameters class
      client::Parameters p(argv[1]);

      //!Set maximum number of threads
      #ifdef _OPENMP
        if (p.omp_threads>0) omp_set_num_threads(p.omp_threads);
      #endif

      //!Initialize number of dimensions and execute main_(const ParameterType2) accordingly
      #ifdef DIMENSION_2
        if (p.num_dimensions == 2)
          main_<2, client::Parameters> (p);
      #endif

      #ifdef DIMENSION_3
        if (p.num_dimensions == 3)
          main_<3, client::Parameters> (p);
      #endif

      double exec_time = my::time::GetTime()-timer;
      std::stringstream ss;
      ss << exec_time;
      msg::print_message_2("Finished - exec-time: "+ss.str()+" s");
    }
  }
  else{
    std::cout << "No option or parameter file passed.\nSee -h or --help for more information." << std::endl;
  }

  return 0;

}
