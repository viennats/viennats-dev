#ifndef LEVELSET2VOLUME_HPP_
#define LEVELSET2VOLUME_HPP_


#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkXMLRectilinearGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkPolyData.h>
#include <vtkGeometryFilter.h>
#include <vtkBooleanOperationPolyDataFilter.h>
#include <vtkTableBasedClipDataSet.h>
#include <vtkDataSetTriangleFilter.h>
#include <vtkAppendFilter.h>
#include "vector.hpp"

// TODO remove
#include <vtkXMLPolyDataWriter.h>


namespace lvlset{


  // This function creates a vtkRectilinearGrid and fills it with the level set values
  // Then it uses ClipDataSetWithPolyData to create tetra meshes for all wrapped materials
  // Explicit boolean operations are then used to create the final layers in a single mesh
  template<class LevelSetsType>
  void extract_volume(const LevelSetsType& LevelSets, vtkSmartPointer<vtkUnstructuredGrid>& volumeMesh){

    typedef typename LevelSetsType::value_type LevelSetType;
    typedef typename LevelSetType::grid_type2 GridType;
    static const int D=GridType::dimensions;

    constexpr bool debugOutput=true;

    // number of levels required to output
    constexpr int numLayers = 2;

    // store volume for each material
    std::vector< vtkSmartPointer<vtkUnstructuredGrid> > materialMeshes;

    unsigned counter = 0; // TODO remove

    // create tetrahedralised volume for each level set
    // start with back of list to make subtraction less memory intensive
    for(typename LevelSetsType::const_reverse_iterator it=LevelSets.rbegin(); it!=LevelSets.rend(); ++it){
      if(it->num_active_pts() == 0) continue; //ignore empty levelsets

      if(debugOutput) std::cout << "Starting levelset " << counter << std::endl;

      assert(it->number_of_layers()>=numLayers);  // check if enough information

      // get extent of rectilinear grid needed
      double gridDelta = it->grid().grid_delta();
      lvlset::vec<int, D> gridMin(0), gridMax(0); // intialise as origin vector

      for(unsigned i=0; i<D; ++i){
        if(it->grid().boundary_conditions(i)==lvlset::INFINITE_BOUNDARY){
          gridMin[i] = it->get_min_runbreak(i);
          gridMax[i] = it->get_max_runbreak(i);
        }else{
          gridMin[i] = it->grid().min_grid_index(i);
          gridMax[i] = it->grid().max_grid_index(i);
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
      typename LevelSetType::const_iterator_runs it_l(*it, gridMin);


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
        vec<typename LevelSetType::index_type, D> indices(it->grid().global_coordinates_2_global_indices(p));

        // if(debugOutput){
        //   std::cout << "Grid point: " << indices << " / " << it_l.end_indices() << " value: " << it_l.value() << std::endl;
        // }

        typename LevelSetType::value_type LSValue = (it_l.is_finished())?(LevelSetType::POS_VALUE):it_l.value();
        if(LSValue == LevelSetType::POS_VALUE){
          LSValue = numLayers;
        }else if(LSValue == LevelSetType::NEG_VALUE){
          LSValue = - numLayers;
        }
        signedDistances->InsertNextValue(LSValue * it->grid().grid_delta());


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
        gwriter->SetFileName(("grid_" + std::to_string(counter) + ".vtr").c_str());
        gwriter->SetInputData(rgrid);
        gwriter->Write();
      }

      // Use vtkClipDataSet to slice the grid
      vtkSmartPointer<vtkTableBasedClipDataSet> clipper = //TODO change to vtkTableBasedClipDataSet
        vtkSmartPointer<vtkTableBasedClipDataSet>::New();
      clipper->SetInputData(rgrid);
      clipper->InsideOutOn();
      clipper->Update();

      // convert generated clip to only include tetras
      vtkSmartPointer<vtkDataSetTriangleFilter> triangleFilter = vtkDataSetTriangleFilter::New();
      triangleFilter->SetInputConnection(clipper->GetOutputPort());
      triangleFilter->Update();

      // transform grid to polydata
      // vtkSmartPointer<vtkGeometryFilter> geometryFilter =
      //   vtkSmartPointer<vtkGeometryFilter>::New();
      // geometryFilter->SetInputConnection(clipper->GetOutputPort());
      // geometryFilter->Update();

      materialMeshes.push_back(triangleFilter->GetOutput());

      // if(it!=LevelSets.rbegin()){
      //   // now subtract second to last material with the current one
      //
      //
      //   materialMeshes.rbegin()[1] = newLayer;
      // }


      // subtract material below for all but first element
      // if(it!=LevelSets.rbegin()){
      //   if(debugOutput) std::cout << "Booling levelset: " << counter << std::endl;
      //   vtkSmartPointer<vtkBooleanOperationPolyDataFilter> booleanOperation =
      //     vtkSmartPointer<vtkBooleanOperationPolyDataFilter>::New();
      //   booleanOperation->SetOperationToDifference();
      //
      //   // subtract material below from this material
      //   booleanOperation->SetInputData( 0, materialMeshes[materialMeshes.size()-1]);
      //   booleanOperation->SetInputData( 1, materialMeshes[materialMeshes.size()-2]);
      //
      //   materialMeshes[materialMeshes.size()-2] = booleanOperation->GetOutput();
      // }

      // TODO remove; increment counter
      ++counter;
    }

    vtkSmartPointer<vtkAppendFilter> appendFilter =
      vtkSmartPointer<vtkAppendFilter>::New();
    //appendFilter->MergePointsOn();

    for(int i=materialMeshes.size()-1; i>=0; --i){
      vtkSmartPointer<vtkXMLUnstructuredGridWriter> owriter =
        vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
      owriter->SetFileName(("material_" + std::to_string(i) + ".vtu").c_str());
      owriter->SetInputData(materialMeshes[i]);
      owriter->Write();

      appendFilter->AddInputData(materialMeshes[i]);
    }

    appendFilter->Update();
    volumeMesh = appendFilter->GetOutput();
  }
}

#endif
