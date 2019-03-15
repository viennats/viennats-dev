#ifndef LEVELSET2VOLUME_HPP_
#define LEVELSET2VOLUME_HPP_


#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkCellData.h>
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
  vtkSmartPointer<vtkRectilinearGrid> LS2RectiLinearGrid(const LevelSetType& LevelSet, const int gridExtraPoints=0, const double LSOffset=0.){
    typedef typename LevelSetType::grid_type2 GridType;
    static const int D=GridType::dimensions;

    // check if enough information
    constexpr int minLayers = 2;
    const int numLayers = LevelSet.number_of_layers();

    assert(numLayers >= minLayers);

    // get extent of rectilinear grid needed
    double gridDelta = LevelSet.grid().grid_delta();

    vtkSmartPointer<vtkFloatArray> coords[D];
    lvlset::vec<int, D> VgridMin(0); // intialise as origin vector
    int gridMin=0, gridMax=0;


    // fill grid with offset depending on orientation
    for(unsigned i=0; i<D; ++i){
      coords[i] = vtkSmartPointer<vtkFloatArray>::New();

      if(LevelSet.grid().boundary_conditions(i)==lvlset::INFINITE_BOUNDARY){
        gridMin = LevelSet.get_min_runbreak(i)-1;  // make bigger by one for numerical stability
        gridMax = LevelSet.get_max_runbreak(i)+1;
        VgridMin[i] = gridMin;

        for(int x=gridMin; x <= gridMax; ++x){
          coords[i]->InsertNextValue(x * gridDelta);
        }

      }else{
        gridMin = LevelSet.grid().min_grid_index(i)-gridExtraPoints;
        gridMax = LevelSet.grid().max_grid_index(i)+gridExtraPoints;
        VgridMin[i] = gridMin;

        for(int x=gridMin; x<=gridMax;++x){
          coords[i]->InsertNextValue(x * gridDelta);
        }
      }
    }


    vtkSmartPointer<vtkRectilinearGrid> rgrid =
      vtkSmartPointer<vtkRectilinearGrid>::New();

    rgrid->SetDimensions(coords[0]->GetNumberOfTuples(),
                     coords[1]->GetNumberOfTuples(),
                     coords[2]->GetNumberOfTuples());
    rgrid->SetXCoordinates(coords[0]);
    rgrid->SetYCoordinates(coords[1]);
    rgrid->SetZCoordinates(coords[2]);


    // now iterate over grid and fill with LS values
    // initialise iterator over levelset and pointId to start at gridMin
    vtkIdType pointId=0;
    typename LevelSetType::const_iterator_runs it_l(LevelSet, VgridMin);

    // Make array to store signed distance function
    // Create an array to hold distance information
    vtkSmartPointer<vtkFloatArray> signedDistances =
      vtkSmartPointer<vtkFloatArray>::New();
    signedDistances->SetNumberOfComponents(1);
    signedDistances->SetName("SignedDistances");

    //iterate until both ends have been reached
    while( (pointId < rgrid->GetNumberOfPoints()) ){
      double p[3];
      rgrid->GetPoint(pointId, p);
      // create index vector
      vec<typename LevelSetType::index_type, D> indices(LevelSet.grid().global_coordinates_2_global_indices(p));

      // write the corresponding LSValue
      typename LevelSetType::value_type LSValue;

      // if indices are outside of domain mark point with max value type
      if(LevelSet.grid().is_outside_domain(indices)){
        signedDistances->InsertNextValue(signedDistances->GetDataTypeValueMax());
      }else{
        // if inside domain just write the correct value
        LSValue = (it_l.is_finished())?(LevelSetType::POS_VALUE):it_l.value()+LSOffset;
        if(LSValue == LevelSetType::POS_VALUE){
          LSValue = numLayers;
        }else if(LSValue == LevelSetType::NEG_VALUE){
          LSValue = - numLayers;
        }

        signedDistances->InsertNextValue(LSValue * gridDelta);
      }

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

    // now need to go through again to fix border points, this is done by mapping existing points onto the points outside of the domain according to the correct boundary conditions
    if(gridExtraPoints!=0){
      pointId = 0;
      while( (pointId < rgrid->GetNumberOfPoints()) ){

        if(signedDistances->GetValue(pointId) == signedDistances->GetDataTypeValueMax()){
          double p[3];
          rgrid->GetPoint(pointId, p);
          // create index vector
          vec<typename LevelSetType::index_type, D> indices(LevelSet.grid().global_coordinates_2_global_indices(p));

          vec<typename LevelSetType::index_type,D> localIndices = LevelSet.grid().global_indices_2_local_indices(indices);

          // now find Id of point we need to take value from ( always assume reflective for now)
          int originalPointId=0;
          for(int i=D-1; i>=0; --i){
            originalPointId*=(LevelSet.grid().grid_extent(i)+2*gridExtraPoints+1);
            originalPointId += localIndices[i]-indices[i];
          }
          originalPointId+=pointId;

          // now put value of mapped point in global point
          signedDistances->SetValue(pointId, signedDistances->GetValue(originalPointId));
        }
        ++pointId;
      }
    }

    // Add the SignedDistances to the grid
    rgrid->GetPointData()->SetScalars(signedDistances);

    return rgrid;
  }







  // This function creates a vtkRectilinearGrid and fills it with the level set values
  // Then it uses ClipDataSetWithPolyData to create tetra meshes for all wrapped materials
  // Explicit boolean operations are then used to create the final layers in a single mesh
  template<class LevelSetsType>
  void extract_volume(const LevelSetsType& LevelSets, vtkSmartPointer<vtkUnstructuredGrid>& volumeMesh){
    // set debug option
    // constexpr bool debugOutput=true;

    // number of levels required to output
    constexpr int numLayers = 2;

    // const double gridDelta = LevelSets.begin()->grid().grid_delta();

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

      // if(debugOutput) std::cout << "Cutting levelset " << counter << std::endl;

      assert(it->number_of_layers()>=numLayers);  // check if enough information


      // create grid of next LS with slight offset and project into current mesh
      vtkSmartPointer<vtkRectilinearGrid> rgrid = vtkSmartPointer<vtkRectilinearGrid>::New();
      rgrid = LS2RectiLinearGrid(*it, 1, -1e-5);  // number of extra grid points outside and LSOffset

      // if(debugOutput){
      //   vtkSmartPointer<vtkXMLRectilinearGridWriter> gwriter =
      //     vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
      //   gwriter->SetFileName(("grid" + std::to_string(counter) + ".vtr").c_str());
      //   gwriter->SetInputData(rgrid);
      //   gwriter->Write();
      // }

      // now transfer implicit values to mesh points
      vtkSmartPointer<vtkProbeFilter> probeFilter = vtkSmartPointer<vtkProbeFilter>::New();
      probeFilter->SetInputData(*(materialMeshes.rbegin()));  // last element
      probeFilter->SetSourceData(rgrid);
      probeFilter->Update();

      // if(debugOutput){
      //   vtkSmartPointer<vtkXMLUnstructuredGridWriter> owriter =
      //     vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
      //   owriter->SetFileName(("probed_" + std::to_string(counter) + ".vtu").c_str());
      //   owriter->SetInputData(probeFilter->GetOutput());
      //   owriter->Write();
      // }

      // now clip the mesh and save the clipped as the 1st layer and use the inverse for the next layer clipping
      // Use vtkTabelBasedClipDataSet to slice the grid
      vtkSmartPointer<vtkTableBasedClipDataSet> insideClipper =
        vtkSmartPointer<vtkTableBasedClipDataSet>::New();
      insideClipper->SetInputConnection(probeFilter->GetOutputPort());
      insideClipper->GenerateClippedOutputOn();
      // insideClipper->UseValueAsOffsetOff();
      // insideClipper->SetValue(1.5 * gridDelta);
      insideClipper->Update();

      materialMeshes.rbegin()[0] = insideClipper->GetOutput();
      materialMeshes.push_back(insideClipper->GetClippedOutput());

      // if(debugOutput){
      //   vtkSmartPointer<vtkXMLUnstructuredGridWriter> owriter =
      //     vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
      //   owriter->SetFileName(("layer_" + std::to_string(counter) + ".vtu").c_str());
      //   owriter->SetInputData(insideClipper->GetOutput());
      //   owriter->Write();
      //
      //   owriter->SetFileName(("clipped_away_" + std::to_string(counter) + ".vtu").c_str());
      //   owriter->SetInputData(insideClipper->GetClippedOutput());
      //   owriter->Write();
      // }


      // TODO remove; increment counter
      ++counter;
    }

    vtkSmartPointer<vtkAppendFilter> appendFilter =
      vtkSmartPointer<vtkAppendFilter>::New();
    //appendFilter->MergePointsOn();

    for(unsigned i=0; i<materialMeshes.size(); ++i){ // TODO change to use reverse iterator
      // if(debugOutput){
      //   vtkSmartPointer<vtkXMLUnstructuredGridWriter> owriter =
      //     vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
      //   owriter->SetFileName(("material_" + std::to_string(i) + ".vtu").c_str());
      //   owriter->SetInputData(materialMeshes[materialMeshes.size()-1-i]);
      //   owriter->Write();
      // }

      // write material number in mesh
      vtkSmartPointer<vtkIntArray> materialNumberArray = vtkSmartPointer<vtkIntArray>::New();
      materialNumberArray->SetNumberOfComponents(1);
      materialNumberArray->SetName("Material");
      for(unsigned j=0; j<materialMeshes[materialMeshes.size()-1-i]->GetNumberOfCells(); ++j){
        materialNumberArray->InsertNextValue(i+1);
      }
      materialMeshes[materialMeshes.size()-1-i]->GetCellData()->SetScalars(materialNumberArray);

      // delete all cell data, so it is not in ouput
      // TODO this includes signed distance information which could be conserved for debugging
      // also includes wheter a cell was vaild for cutting by the grid
      vtkSmartPointer<vtkPointData> pointData = materialMeshes[i]->GetPointData();
      const int numberOfArrays = pointData->GetNumberOfArrays();
      for(int j=0; j<numberOfArrays; ++j){
        pointData->RemoveArray(0); // remove first array until none are left
      }

      appendFilter->AddInputData(materialMeshes[i]);
    }

    appendFilter->Update();

    // if(debugOutput) std::cout << "Triangulating volume mesh" << std::endl;

    // change all 3D cells into tetras and all 2D cells to triangles
    vtkSmartPointer<vtkDataSetTriangleFilter> triangleFilter =
      vtkSmartPointer<vtkDataSetTriangleFilter>::New();
    triangleFilter->SetInputConnection(appendFilter->GetOutputPort());
    triangleFilter->Update();

    volumeMesh = triangleFilter->GetOutput();
  }
}

#endif
