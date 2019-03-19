#ifndef LEVELSET2VOLUME_HPP_
#define LEVELSET2VOLUME_HPP_


#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkTableBasedClipDataSet.h>
#include <vtkDataSetTriangleFilter.h>
#include <vtkAppendFilter.h>
#include <vtkProbeFilter.h>

#include "vector.hpp"



namespace lvlset{

  // This function takes a levelset and converts it to a vtkRectilinearGrid
  // The full domain contains values, which are capped at numLayers * gridDelta
  // gridExtraPoints layers of grid points are added to the domain according to boundary conditions
  template<int gridExtraPoints=0, class LevelSetType>
  vtkSmartPointer<vtkRectilinearGrid> LS2RectiLinearGrid(const LevelSetType& LevelSet, const double LSOffset=0., int infiniteMinimum=std::numeric_limits<int>::max()){
    typedef typename LevelSetType::grid_type2 GridType;
    static const int D=GridType::dimensions;

    // check if enough information
    constexpr int minLayers = 2;
    const int numLayers = LevelSet.number_of_layers();

    assert(numLayers >= minLayers);

    // get extent of rectilinear grid needed
    double gridDelta = LevelSet.grid().grid_delta();

    vtkSmartPointer<vtkFloatArray> coords[3]; // needs to be 3 because vtk only knows 3D
    int gridMin=0, gridMax=0;


    // fill grid with offset depending on orientation
    for(unsigned i=0; i<D; ++i){
      coords[i] = vtkSmartPointer<vtkFloatArray>::New();

      if(LevelSet.grid().boundary_conditions(i)==lvlset::INFINITE_BOUNDARY){
        // add one to gridMin and gridMax for numerical stability
        gridMin = std::min(LevelSet.get_min_runbreak(i), infiniteMinimum)-1;  // choose the smaller number so that for first levelset the overall minimum can be chosen
        gridMax = LevelSet.get_max_runbreak(i)+1;

        for(int x=gridMin; x <= gridMax; ++x){
          coords[i]->InsertNextValue(x * gridDelta);
        }

      }else{
        gridMin = LevelSet.grid().min_grid_index(i)-gridExtraPoints;
        gridMax = LevelSet.grid().max_grid_index(i)+gridExtraPoints;
        // VgridMin[i] = gridMin;

        for(int x=gridMin; x<=gridMax;++x){
          coords[i]->InsertNextValue(x * gridDelta);
        }
      }
    }

    // if we work in 2D, just add 1 grid point at origin
    if(D==2){
      coords[2] = vtkSmartPointer<vtkFloatArray>::New();
      coords[2]->InsertNextValue(0);
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
    bool fixBorderPoints = (gridExtraPoints!=0);

    typename LevelSetType::const_iterator_runs it_l(LevelSet);

    // Make array to store signed distance function
    vtkSmartPointer<vtkFloatArray> signedDistances =
      vtkSmartPointer<vtkFloatArray>::New();
    signedDistances->SetNumberOfComponents(1);
    signedDistances->SetName("SignedDistances");

    //iterate until all grid points have a signed distance value
    while( (pointId < rgrid->GetNumberOfPoints()) ){
      double p[3];
      rgrid->GetPoint(pointId, p);
      // create index vector
      vec<typename LevelSetType::index_type, D> indices(LevelSet.grid().global_coordinates_2_global_indices(p));

      // write the corresponding LSValue
      typename LevelSetType::value_type LSValue;

      // if indices are outside of domain mark point with max value type
      if(LevelSet.grid().is_outside_domain(indices)){
        fixBorderPoints=true;
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
      if(it_l.is_finished()){ // when iterator is done just fill all the remaining points
        ++pointId;
      }
      else{
        // advance iterator until it is at correct point
        while(compare(it_l.end_indices(), indices) < 0){
          it_l.next();
          if(it_l.is_finished()) break;
        }
        // now move iterator with pointId to get to next point
        switch(compare(it_l.end_indices(), indices)) {
            case 0:
                it_l.next();
            default:
                ++pointId;
        }
      }
    }

    // now need to go through again to fix border points, this is done by mapping existing points onto the points outside of the domain according to the correct boundary conditions
    if(fixBorderPoints){
      pointId = 0;
      while( (pointId < rgrid->GetNumberOfPoints()) ){
        if(signedDistances->GetValue(pointId) == signedDistances->GetDataTypeValueMax()){
          double p[3];
          rgrid->GetPoint(pointId, p);

          // create index vector
          vec<typename LevelSetType::index_type, D> indices(LevelSet.grid().global_coordinates_2_global_indices(p));

          // vector for mapped point inside domain
          vec<typename LevelSetType::index_type,D> localIndices = LevelSet.grid().global_indices_2_local_indices(indices);

          // now find Id of point we need to take value from
          int originalPointId=0;
          for(int i=D-1; i>=0; --i){
            originalPointId *= coords[i]->GetNumberOfTuples(); // extent in direction
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







  // This function takes the biggest/top levelset and creates a mesh of it
  // since the top levelset wraps all lower ones, all levelsets below it are used to cut this mesh
  // explicitly. When all the cuts are made, material numbers are assigned and the tetra meshes written
  // to one single file
  template<class LevelSetsType>
  void extract_volume(const LevelSetsType& LevelSets, vtkSmartPointer<vtkUnstructuredGrid>& volumeMesh){

    typedef typename LevelSetsType::value_type::grid_type2 GridType;
    static const int D=GridType::dimensions;

    // number of levels required to output
    constexpr int numLayers = 2;

    // store volume for each material
    std::vector< vtkSmartPointer<vtkUnstructuredGrid> > materialMeshes;


    int totalMinimum = std::numeric_limits<int>::max();
    for(auto it=LevelSets.begin(); it!=LevelSets.end(); ++it){
      for(unsigned i=0; i<D; ++i){
        if(it->grid().boundary_conditions(i)==lvlset::INFINITE_BOUNDARY){
          totalMinimum = std::min(totalMinimum, it->get_min_runbreak(i));
        }
      }
    }


    // create volume mesh for largest LS
    // Use vtkClipDataSet to slice the grid
    vtkSmartPointer<vtkTableBasedClipDataSet> clipper =
      vtkSmartPointer<vtkTableBasedClipDataSet>::New();
    clipper->SetInputData(LS2RectiLinearGrid(*(LevelSets.rbegin()), 0, totalMinimum));  // last LS
    clipper->InsideOutOn();
    clipper->Update();


    materialMeshes.push_back(clipper->GetOutput());

    // now cut large volume mesh with all the smaller ones
    for(typename LevelSetsType::const_reverse_iterator it=++LevelSets.rbegin(); it!=LevelSets.rend(); ++it){
      if(it->num_active_pts() == 0) continue; //ignore empty levelsets

      assert(it->number_of_layers()>=numLayers);  // check if enough information

      // create grid of next LS with slight offset and project into current mesh
      vtkSmartPointer<vtkRectilinearGrid> rgrid = vtkSmartPointer<vtkRectilinearGrid>::New();
      rgrid = LS2RectiLinearGrid<1>(*it, -1e-5);  // number of extra grid points outside and LSOffset

      // now transfer implicit values to mesh points
      vtkSmartPointer<vtkProbeFilter> probeFilter = vtkSmartPointer<vtkProbeFilter>::New();
      probeFilter->SetInputData(*(materialMeshes.rbegin()));  // last element
      probeFilter->SetSourceData(rgrid);
      probeFilter->Update();

      // now clip the mesh and save the clipped as the 1st layer and use the inverse for the next layer clipping
      // Use vtkTabelBasedClipDataSet to slice the grid
      vtkSmartPointer<vtkTableBasedClipDataSet> insideClipper =
        vtkSmartPointer<vtkTableBasedClipDataSet>::New();
      insideClipper->SetInputConnection(probeFilter->GetOutputPort());
      insideClipper->GenerateClippedOutputOn();
      insideClipper->Update();

      materialMeshes.rbegin()[0] = insideClipper->GetOutput();
      materialMeshes.push_back(insideClipper->GetClippedOutput());
    }

    vtkSmartPointer<vtkAppendFilter> appendFilter =
      vtkSmartPointer<vtkAppendFilter>::New();
    appendFilter->MergePointsOn();


    for(unsigned i=0; i<materialMeshes.size(); ++i){

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

    // change all 3D cells into tetras and all 2D cells to triangles
    vtkSmartPointer<vtkDataSetTriangleFilter> triangleFilter =
      vtkSmartPointer<vtkDataSetTriangleFilter>::New();
    triangleFilter->SetInputConnection(appendFilter->GetOutputPort());
    triangleFilter->Update();

    volumeMesh = triangleFilter->GetOutput();
  }
}

#endif
