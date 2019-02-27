#ifndef LEVELSET2VOLUME_HPP_
#define LEVELSET2VOLUME_HPP_


#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkXMLRectilinearGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkPolyData.h>
#include <vtkGeometryFilter.h>
#include <vtkBooleanOperationPolyDataFilter.h>
#include <vtkTableBasedClipDataSet.h>
#include <vtkDataSetTriangleFilter.h>
#include <vtkAppendFilter.h>
#include <vtkProbeFilter.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include "vector.hpp"

// TODO remove
#include <vtkXMLPolyDataWriter.h>


namespace lvlset{

  // This function takes a levelset and converts it into a vtkRectilinearGrid
  // The full domain contains values, which are capped at numLayers * gridDelta
  // the offset just increases the size of the grid so it is mapped correctly to the mesh
  template<class LevelSetType>
  vtkSmartPointer<vtkRectilinearGrid> LS2RectiLinearGrid(const LevelSetType& LevelSet, const double gridScaleOffset=0.){
    typedef typename LevelSetType::grid_type2 GridType;
    static const int D=GridType::dimensions;

    constexpr bool debugOutput = true;

    // check if enough information
    constexpr int minLayers = 2;
    const int numLayers = LevelSet.number_of_layers();

    assert(numLayers >= minLayers);

    // get extent of rectilinear grid needed
    double gridDelta = LevelSet.grid().grid_delta();
    lvlset::vec<int, D> gridMin(0), gridMax(0); // intialise as origin vector

    for(unsigned i=0; i<D; ++i){
      if(LevelSet.grid().boundary_conditions(i)==lvlset::INFINITE_BOUNDARY){
        gridMin[i] = LevelSet.get_min_runbreak(i);
        gridMax[i] = LevelSet.get_max_runbreak(i);
      }else{
        gridMin[i] = LevelSet.grid().min_grid_index(i);
        gridMax[i] = LevelSet.grid().max_grid_index(i);
      }
    }

    // make grid
    vtkSmartPointer<vtkFloatArray> xCoords =
    vtkSmartPointer<vtkFloatArray>::New();
    for(int x=gridMin[0]; x <= gridMax[0]; ++x){
      xCoords->InsertNextValue(x * gridDelta);
    }
    // y
    vtkSmartPointer<vtkFloatArray> yCoords =
    vtkSmartPointer<vtkFloatArray>::New();
    for(int x=gridMin[1]; x <= gridMax[1]; ++x){
      yCoords->InsertNextValue(x * gridDelta);
    }
    // z
    vtkSmartPointer<vtkFloatArray> zCoords =
    vtkSmartPointer<vtkFloatArray>::New();
    if(D==3){
      for(int x=gridMin[2]; x <= gridMax[2]; ++x){
        zCoords->InsertNextValue(x * gridDelta);
      }
    }else{
      zCoords->InsertNextValue(0);
    }


    vtkSmartPointer<vtkRectilinearGrid> rgrid =
      vtkSmartPointer<vtkRectilinearGrid>::New();

    rgrid->SetDimensions(xCoords->GetNumberOfTuples(),
                     yCoords->GetNumberOfTuples(),
                     zCoords->GetNumberOfTuples());
    rgrid->SetXCoordinates(xCoords);
    rgrid->SetYCoordinates(yCoords);
    rgrid->SetZCoordinates(zCoords);

    // now iterate over grid and fill with LS values
    // initialise iterator over levelset and pointId to start at gridMin
    vtkIdType pointId=0;
    typename LevelSetType::const_iterator_runs it_l(LevelSet, gridMin);


    // Make array to store signed distance function
    // Create an array to hold distance information
    vtkSmartPointer<vtkFloatArray> signedDistances =
      vtkSmartPointer<vtkFloatArray>::New();
    signedDistances->SetNumberOfComponents(1);
    signedDistances->SetName("SignedDistances");

    //iterate until both ends have been reached
    // while( (pointId < rgrid->GetNumberOfPoints()) || !it_l.is_finished())
    while( (pointId < rgrid->GetNumberOfPoints()) ){
      // std::cout << "ID<Num: " << (pointId < rgrid->GetNumberOfPoints()) << ", !is finished: " << !it_l.is_finished() << std::endl;
      double p[3];
      rgrid->GetPoint(pointId, p);
      // create index vector
      vec<typename LevelSetType::index_type, D> indices(LevelSet.grid().global_coordinates_2_global_indices(p));

      // if(debugOutput){
      //   std::cout << "Grid point: " << indices << " / " << it_l.end_indices() << " value: " << it_l.value() << std::endl;
      // }

      typename LevelSetType::value_type LSValue = (it_l.is_finished())?(LevelSetType::POS_VALUE):it_l.value();
      if(LSValue == LevelSetType::POS_VALUE){
        LSValue = numLayers;
      }else if(LSValue == LevelSetType::NEG_VALUE){
        LSValue = - numLayers;
      }
      signedDistances->InsertNextValue(LSValue * gridDelta);


      // move iterator if point was visited
      switch(compare(it_l.end_indices(), indices)) {
          case -1:
              it_l.next();
              break;
          case 0:
              it_l.next();
          default:
              ++pointId;
      }
    }

    // Add the SignedDistances to the grid
    rgrid->GetPointData()->SetScalars(signedDistances);

    if(debugOutput){
      vtkSmartPointer<vtkXMLRectilinearGridWriter> gwriter =
        vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
      gwriter->SetFileName("grid.vtr");
      gwriter->SetInputData(rgrid);
      gwriter->Write();
    }

    return rgrid;
  }







  // This function creates a vtkRectilinearGrid and fills it with the level set values
  // Then it uses ClipDataSetWithPolyData to create tetra meshes for all wrapped materials
  // Explicit boolean operations are then used to create the final layers in a single mesh
  template<class LevelSetsType>
  void extract_volume(const LevelSetsType& LevelSets, vtkSmartPointer<vtkUnstructuredGrid>& volumeMesh){

    // typedef typename LevelSetsType::value_type LevelSetType;
    // typedef typename LevelSetType::grid_type2 GridType;
    // static const int D=GridType::dimensions;

    constexpr bool debugOutput=true;

    // number of levels required to output
    constexpr int numLayers = 2;



    // store volume for each material
    std::vector< vtkSmartPointer<vtkUnstructuredGrid> > materialMeshes;

    unsigned counter = 1; // TODO remove

    // create volume mesh for largest LS
    // Use vtkClipDataSet to slice the grid
    vtkSmartPointer<vtkTableBasedClipDataSet> clipper =
      vtkSmartPointer<vtkTableBasedClipDataSet>::New();
    clipper->SetInputData(LS2RectiLinearGrid(*(LevelSets.rbegin())));  // last LS
    clipper->InsideOutOn();
    clipper->Update();

    materialMeshes.push_back(clipper->GetOutput());

    // now cut large volume mesh with all the smaller ones
    for(typename LevelSetsType::const_reverse_iterator it=++LevelSets.rbegin(); it!=LevelSets.rend(); ++it){
      if(it->num_active_pts() == 0) continue; //ignore empty levelsets

      if(debugOutput) std::cout << "Cutting levelset " << counter << std::endl;

      assert(it->number_of_layers()>=numLayers);  // check if enough information


      // create grid of next LS with slight offset and project into current mesh
      vtkSmartPointer<vtkRectilinearGrid> rgrid = vtkSmartPointer<vtkRectilinearGrid>::New();
      rgrid = LS2RectiLinearGrid(*it, 1e-3);  // relative offset wrt gridDelta

      // now transfer implicit values to mesh points
      vtkSmartPointer<vtkProbeFilter> probeFilter = vtkSmartPointer<vtkProbeFilter>::New();
      probeFilter->SetInputData(*(materialMeshes.rbegin()));  // last element
      probeFilter->SetSourceData(rgrid);
      probeFilter->CategoricalDataOn();
      //probeFilter->ComputeToleranceOn();
      probeFilter->Update();

      if(debugOutput){
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> owriter =
          vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        owriter->SetFileName(("probed_" + std::to_string(counter) + ".vtu").c_str());
        owriter->SetInputData(probeFilter->GetOutput());
        owriter->Write();
      }

      // now clip the mesh and save the clipped as the 1st layer and use the inverse for the next layer clipping
      // Use vtkTabelBasedClipDataSet to slice the grid
      vtkSmartPointer<vtkTableBasedClipDataSet> insideClipper =
        vtkSmartPointer<vtkTableBasedClipDataSet>::New();
      insideClipper->SetInputConnection(probeFilter->GetOutputPort());
      insideClipper->GenerateClippedOutputOn();
      insideClipper->Update();

      materialMeshes.rbegin()[0] = insideClipper->GetOutput();
      materialMeshes.push_back(insideClipper->GetClippedOutput());

      if(debugOutput){
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> owriter =
          vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        owriter->SetFileName(("layer_" + std::to_string(counter) + ".vtu").c_str());
        owriter->SetInputData(insideClipper->GetOutput());
        owriter->Write();

        owriter->SetFileName(("clipped_away_" + std::to_string(counter) + ".vtu").c_str());
        owriter->SetInputData(insideClipper->GetClippedOutput());
        owriter->Write();
      }


      // TODO remove; increment counter
      ++counter;
    }

    vtkSmartPointer<vtkAppendFilter> appendFilter =
      vtkSmartPointer<vtkAppendFilter>::New();
    //appendFilter->MergePointsOn();

    for(int i=materialMeshes.size()-1; i>=0; --i){
      if(debugOutput){
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> owriter =
          vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        owriter->SetFileName(("material_" + std::to_string(i) + ".vtu").c_str());
        owriter->SetInputData(materialMeshes[i]);
        owriter->Write();
      }

      appendFilter->AddInputData(materialMeshes[i]);
    }

    appendFilter->Update();
    volumeMesh = appendFilter->GetOutput();
  }
}

#endif
