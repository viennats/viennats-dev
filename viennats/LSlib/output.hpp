#ifndef OUTPUT_HPP_
#define OUTPUT_HPP_

/* =========================================================================
   Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.

                            -----------------
                 ViennaTS - The Vienna Topography Simulator
                            -----------------

   Contact:         viennats@iue.tuwien.ac.at

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


#include <vector>
#include <fstream>
#include <cstdint>

#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>

#include "kernel.hpp"
#include "levelset2surface.hpp"
#include "levelset2volume.hpp"
#include "../message.h"

//Options for levelset output
#define LVST_FILE_VERSION_NUMBER 1

//Limits for 3,5,6,7 byte output; only change if one byte does not have 8 bits
#define UINT24_MAX 16777215 // highest value of 3 unsigned bytes  2^24-1
#define  INT24_MAX  8388607 //  highest value of 3 signed   bytes  2^23-1
#define  INT24_MIN -8388608 //  lowest  value of 3 signed   bytes -2^23

#define UINT40_MAX 1099511627775UL // highest value of 3 unsigned bytes  2^40-1
#define  INT40_MAX  549755813887LL //  highest value of 3 signed   bytes  2^39-1
#define  INT40_MIN -549755813888LL //  lowest  value of 3 signed   bytes -2^39

#define UINT48_MAX 281474976710655ULL  // highest value of 3 unsigned bytes  2^48-1
#define  INT48_MAX  140737488355327LL //  highest value of 3 signed   bytes  2^47-1
#define  INT48_MIN -140737488355328LL //  lowest  value of 3 signed   bytes -2^47

#define UINT56_MAX 72057594037927935ULL  // highest value of 3 unsigned bytes  2^56-1
#define  INT56_MAX  36028797018963967LL //  highest value of 3 signed   bytes  2^55-1
#define  INT56_MIN -36028797018963968LL //  lowest  value of 3 signed   bytes -2^55

namespace lvlset {

    namespace {

        template <int D> class Surface {  //the surface(-inserter) class used by the output functions defined below
        public:
            std::vector<vec<float,D> > Nodes;
            std::vector<vec<unsigned int,D> > Elements;
            std::vector< std::vector<float> > NodeDataValues;
            std::vector< std::string > NodeDataDescription;

            typedef unsigned int node_ref_type;

            template <class V> node_ref_type insert_node(const V& v) {
                Nodes.push_back(vec<float, D>(v));
                return Nodes.size()-1;
            }

            template <class V> node_ref_type insert_node(const V& v, unsigned int loc) {
                Nodes.push_back(vec<float, D>(v));
                std::vector<vec<float,D> > temp;
                for (unsigned int i=0;i<Nodes.size();i++) temp.push_back(vec<float, D>(Nodes[i]));
                Nodes.clear();
                for (unsigned int i=0;i<temp.size()+1;i++) {
                  if (i==loc) {
                    Nodes.push_back(vec<float, D>(v));
                  }
                  Nodes.push_back(vec<float, D>(temp[i]));
                }
                return Nodes.size()-1;
            }

            template <class V> void insert_element(const V& v) {
                Elements.push_back(vec<unsigned int, D>(v));
            }

        };

      template <int D> class SurfaceList {
      public:

        typedef Surface<D> Surfaces;

        int dim;
        int nsurfaces;
        int open_boundary_direction;
        float bottom;
        std::list<Surfaces > surfaces;

        SurfaceList(int nb_surf):dim(D), nsurfaces(nb_surf) {}
      };
    }

    namespace {

        class DefaultDataType {
        public:
            int number_of_series() const {
                return 0;
            }

            template <class PT_ID_TYPE>
            double get_series_data_double(PT_ID_TYPE active_pt_id, int series) const {
                return double();
            }

            template <class PT_ID_TYPE>
            std::string get_series_data(PT_ID_TYPE active_pt_id, int series) const {
                return std::string();
            }

            std::string get_series_label(int series) const {
                return std::string();
            }

            std::string get_series_type(int series) const {
                return std::string();
            }

            bool get_series_output(int series) const {
                return false;
            }


        };

        template <class SizeType, class DataType>
        class GetActivePointType {
        public:
            typedef std::vector<SizeType> result;
        };


        template<class SizeType>
        class GetActivePointType<SizeType, DefaultDataType> {
        public:
            class result {
            public:

                void push_back(const SizeType& v) const {}

                template<class V>
                SizeType operator[](const V& x) const {
                    return 0;
                }
            };
        };

    }


    template <class GridTraitsType, class LevelSetTraitsType, class DataType>
    void write_explicit_surface_opendx(const levelset<GridTraitsType, LevelSetTraitsType>& l, const std::string& filename, const DataType& Data, typename LevelSetTraitsType::value_type eps=0.) {
        //this function etracts the surface using the marching cubes algorithm and writes it to the specified file using OpenDX-file format
        //the parameter eps is forwarded to the surface extraction function

        const int D=GridTraitsType::dimensions;

        Surface<D> s;

        typename GetActivePointType<typename LevelSetTraitsType::size_type, DataType>::result ActivePointList;

        extract(l, s, eps, ActivePointList);

        std::ofstream f(filename.c_str());

        //!print positions
        f<< "object \"positions\" class array type float rank 1 shape " << D << " items "<< s.Nodes.size() <<" data follows" << std::endl;
        for (unsigned int i=0;i<s.Nodes.size();i++) {
            for (int j=0;j<D;j++) f << static_cast<float>(s.Nodes[i][j]) << " ";
            f<< std::endl;
        }

        //! print connections
        f << "object \"connections\" class array type int rank 1 shape " << D << " items "<< s.Elements.size() <<" data follows" << std::endl;
        for (unsigned int i=0;i<s.Elements.size();i++) {
            for (int j=0;j<D;j++) f<< s.Elements[i][j] << " ";
            f << std::endl;
        }

        if (D==2)
            f << "attribute \"element type\" string \"lines\"" << std::endl;
        else if (D==3)
            f << "attribute \"element type\" string \"triangles\"" << std::endl;
        f << "attribute \"ref\" string \"positions\"" << std::endl;


        //output data
        for (int k=0;k<Data.number_of_series();++k) {
            if (Data.get_series_output(k)) {
                f << "object \"" << Data.get_series_label(k) << "_data\" class array type " << Data.get_series_type(k) << " rank 0 items " << s.Nodes.size() << " data follows" << std::endl;
                for (unsigned int i=0;i<s.Nodes.size();i++) {
                    f << Data.get_series_data(ActivePointList[i],k) << std::endl;
                }
                f << "attribute \"dep\" string \"positions\"" << std::endl;
            }
        }

        //! print profile
        f << "object \"profile\" class field" << std::endl;
        f << "  component \"positions\" value \"positions\"" << std::endl;
        f << "  component \"connections\" value \"connections\"" << std::endl;

        for (int k=0;k<Data.number_of_series();++k) {
            if (Data.get_series_output(k)) {
                f << "object \""<< Data.get_series_label(k) << "\" class field" << std::endl;
                f << "  component \"positions\" value \"positions\"" << std::endl;
                f << "  component \"connections\" value \"connections\"" << std::endl;
                f << "  component \"data\" value \"" << Data.get_series_label(k) << "_data\"" << std::endl;
            }
        }

        f << "end" << std::endl;




        f.close();
    }

    template <class GridTraitsType, class LevelSetTraitsType>
    void write_explicit_surface_opendx(const levelset<GridTraitsType, LevelSetTraitsType>& l, const std::string& filename, typename LevelSetTraitsType::value_type eps=0.) {
        write_explicit_surface_opendx(l, filename, DefaultDataType(), eps);
    }

    template <class LevelSetsType, class CounterType, class ParameterType>
    void write_explicit_volume_vtk(const LevelSetsType& LevelSets, const CounterType counter, ParameterType& p, double eps=0.) {
      write_explicit_volume_vtk(LevelSets, counter, p, DefaultDataType(), eps);
    }

    template <class LevelSetsType, class CounterType, class ParameterType, class DataType>
    void write_explicit_volume_vtk(const LevelSetsType& LevelSets, const CounterType counter, ParameterType& p, const DataType& Data, double eps=0.) {
      static const int D=LevelSetsType::value_type::grid_type2::dimensions;

      vtkSmartPointer<vtkUnstructuredGrid> volumeMesh;
      vtkSmartPointer<vtkPolyData> hullMesh;

      // check if any output should be done. if none, we do not need to extract the volume
      if(!(p.print_volume_tetra || p.print_volume_hull)) return;

      // if volumeMesh or hullMesh are not initialised, they will not be created in extract_volume
      if(p.print_volume_tetra) volumeMesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
      if(p.print_volume_hull) hullMesh = vtkSmartPointer<vtkPolyData>::New();

      // depending on the open boundary, volume output needs to be treated differently if the bottom is removed
      if(p.remove_bottom &&
       (p.open_boundary<(D-1))){
        extract_volume<true>(LevelSets, volumeMesh, hullMesh);
      }else{
        extract_volume<false>(LevelSets, volumeMesh, hullMesh);
      }

      if(p.print_volume_tetra){
        std::ostringstream oss;
        oss << p.output_path << "Volume_" << counter << ".vtu";

        vtkSmartPointer<vtkXMLUnstructuredGridWriter> owriter =
          vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        owriter->SetFileName(oss.str().c_str());
        owriter->SetInputData(volumeMesh);
        owriter->Write();
      }


      if(p.print_volume_hull){
        std::ostringstream oss;
        oss << p.output_path << "Hull_" << counter << ".vtp";

        vtkSmartPointer<vtkXMLPolyDataWriter> pwriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        pwriter->SetFileName(oss.str().c_str());
        pwriter->SetInputData(hullMesh);
        pwriter->Write();
      }
    }

    template <class GridTraitsType, class LevelSetTraitsType>
    void write_explicit_levelset(const levelset<GridTraitsType, LevelSetTraitsType>& l, const std::string& filename){
      const int D=GridTraitsType::dimensions;

      vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
      vtkSmartPointer<vtkPoints> polyPoints = vtkSmartPointer<vtkPoints>::New();
      vtkSmartPointer<vtkCellArray> polyCells = vtkSmartPointer<vtkCellArray>::New();
      vtkSmartPointer<vtkFloatArray> polyValues = vtkSmartPointer<vtkFloatArray>::New();

      polyValues->SetNumberOfComponents(1);
      polyValues->SetName("LSValues");

      // start iterator over LS
      typename levelset<GridTraitsType, LevelSetTraitsType>::const_iterator_runs it_l(l);
      double gridDelta = l.grid().grid_delta();

      while(!it_l.is_finished()){
        if(!it_l.is_active()){
          it_l.next();
          continue;
        }
        double p[3];
        for(unsigned i=0; i<D; ++i) p[i] = gridDelta*it_l.start_indices(i);
        vtkIdType pointId = polyPoints->InsertNextPoint(p);
        polyCells->InsertNextCell(1, &pointId); // insert vertex for visualisation
        polyValues->InsertNextValue(it_l.value());
        it_l.next();
      }

      polyData->SetPoints(polyPoints);
      polyData->SetVerts(polyCells);
      polyData->GetPointData()->SetScalars(polyValues);

      vtkSmartPointer<vtkXMLPolyDataWriter> pWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
      pWriter->SetFileName(filename.c_str());
      pWriter->SetInputData(polyData);
      pWriter->Write();
    }

    template <class GridTraitsType, class LevelSetTraitsType, class DataType>
    void write_explicit_surface_vtp(const levelset<GridTraitsType, LevelSetTraitsType>& l, const std::string& filename, const DataType& Data, typename LevelSetTraitsType::value_type eps=0.) {

      const int D=GridTraitsType::dimensions;
      Surface<D> s;
      typename GetActivePointType<typename LevelSetTraitsType::size_type, DataType>::result ActivePointList;

      extract(l, s, eps, ActivePointList);

      // fill pointlist with points from geometry
      vtkSmartPointer<vtkPoints> polyPoints = vtkSmartPointer<vtkPoints>::New();
      for (unsigned int i=0;i<s.Nodes.size();i++) {
          polyPoints->InsertNextPoint(s.Nodes[i][0], s.Nodes[i][1], s.Nodes[i][2]);
      }

      vtkSmartPointer<vtkCellArray> polyCells = vtkSmartPointer<vtkCellArray>::New();
      // fill cellarray
      for(unsigned int i=0;i<s.Elements.size();i++) {
          polyCells->InsertNextCell(D);
          for (int j=0;j<D;j++){
            polyCells->InsertCellPoint(s.Elements[i][j]);
          }
      }

      vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
      polyData->SetPoints(polyPoints);
      polyData->SetPolys(polyCells);

      // for each data series create point data
      for (int k=0;k<Data.number_of_series();++k) {
        if (Data.get_series_output(k)) {
          vtkSmartPointer<vtkFloatArray> pointData = vtkSmartPointer<vtkFloatArray>::New();
          pointData->SetNumberOfComponents(1);
          pointData->SetName(Data.get_series_label(k).c_str());
          for (unsigned int i=0;i<s.Nodes.size();i++) {
              pointData->InsertNextValue(Data.get_series_data_double(ActivePointList[i],k));
          }
          polyData->GetPointData()->AddArray(pointData);
        }
      }

      vtkSmartPointer<vtkXMLPolyDataWriter> pwriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
      pwriter->SetFileName(filename.c_str());
      pwriter->SetInputData(polyData);
      pwriter->Write();
    }


    template <class GridTraitsType, class LevelSetTraitsType>
    void write_explicit_surface_vtp(const levelset<GridTraitsType, LevelSetTraitsType>& l, const std::string& filename, typename LevelSetTraitsType::value_type eps=0.) {
      write_explicit_surface_vtp(l, filename, DefaultDataType(), eps);
    }

    template <class GridTraitsType, class LevelSetTraitsType, class DataType>
    void write_explicit_surface_vtk(const levelset<GridTraitsType, LevelSetTraitsType>& l, const std::string& filename, const DataType& Data, typename LevelSetTraitsType::value_type eps=0.) {
        //this function etracts the surface using the marching cubes algorithm and writes it to the specified file using vtk-file format
        //the parameter eps is forwarded to the surface extraction function

        const int D=GridTraitsType::dimensions;

        Surface<D> s;

        typename GetActivePointType<typename LevelSetTraitsType::size_type, DataType>::result ActivePointList;

        extract(l, s, eps, ActivePointList);

        std::ofstream f(filename.c_str());

        f << "# vtk DataFile Version 2.0" << std::endl;
        f << D << "D Surface" << std::endl;
        f << "ASCII" << std::endl;
        f << "DATASET UNSTRUCTURED_GRID" << std::endl;
        f << "POINTS " << s.Nodes.size() << " float" << std::endl;

        //!print positions
        for (unsigned int i=0;i<s.Nodes.size();i++) {
            for (int j=0;j<D;j++) f << static_cast<float>(s.Nodes[i][j]) << " ";
            if (D==2) f << "0. ";
            f<< std::endl;
        }
        f << "CELLS " << s.Elements.size() << " " << ((D+1)*s.Elements.size()) << std::endl;
        for (unsigned int i=0;i<s.Elements.size();i++) {
            f << D << " ";
            for (int j=0;j<D;j++) f<< s.Elements[i][j] << " ";
            f << std::endl;
        }

        f << "CELL_TYPES " << s.Elements.size() << std::endl;
        for (unsigned int i=0;i<s.Elements.size();i++) {
            f<< ((D==3)?"5":"3") << std::endl;
        }

        //output data
        f << "POINT_DATA " << s.Nodes.size() << std::endl;

        for (int k=0;k<Data.number_of_series();++k) {
            if (Data.get_series_output(k)) {
                f << "SCALARS " << Data.get_series_label(k) << " " << Data.get_series_type(k) << " 1" << std::endl;
                f << "LOOKUP_TABLE default" << std::endl;
                for (unsigned int i=0;i<s.Nodes.size();i++) {
                    f << Data.get_series_data(ActivePointList[i],k) << std::endl;
                }
            }
        }


        f.close();
    }

    template <class GridTraitsType, class LevelSetTraitsType>
    void write_explicit_surface_vtk(const levelset<GridTraitsType, LevelSetTraitsType>& l, const std::string& filename, typename LevelSetTraitsType::value_type eps=0.) {
        write_explicit_surface_vtk(l, filename, DefaultDataType(), eps);
    }

    template <class TriangulationType>
    void write_triangulated_surface_vtk(const TriangulationType& t, const std::string& filename) {
        //this functions writes a given surface triangulation "t" to the vtk-file "filename"
        //the TriangulationType has to have the same properties as for the levelset-initialization
        //function init (see surface2levelset.hpp)

        const int D=TriangulationType::dimension;

        std::map<typename TriangulationType::node_index_type,unsigned int> node_ids;

        unsigned int counter=0;
        for (typename TriangulationType::element_index_type e=0;e<t.number_of_elements();e++) {
            for (int k=0;k<D;k++) {
                if (node_ids.find(t.element_node_id(e,k))==node_ids.end()) node_ids[t.element_node_id(e,k)]=counter++;
            }
        }

        node_ids.clear();

        std::ofstream f(filename.c_str());

        f << "# vtk DataFile Version 2.0" << std::endl;
        f << D << "D Surface" << std::endl;
        f << "ASCII" << std::endl;
        f << "DATASET UNSTRUCTURED_GRID" << std::endl;
        f << "POINTS " << counter << " float" << std::endl;

        //!print positions

        counter=0;
        for (typename TriangulationType::element_index_type e=0;e<t.number_of_elements();e++) {
            for (int k=0;k<D;k++) {
                if (node_ids.find(t.element_node_id(e,k))==node_ids.end()) {
                    for (int j=0;j<D;j++) f << static_cast<float>(t.node_coordinate(t.element_node_id(e,k),j)) << " ";
                    node_ids[t.element_node_id(e,k)]=counter++;
                    if (D==2) f << "0. ";
                    f<< std::endl;
                }
            }
        }

        f << "CELLS " << t.number_of_elements() << " " << ((D+1)*t.number_of_elements()) << std::endl;
        for (unsigned int i=0;i<t.number_of_elements();i++) {
            f << D << " ";
            for (int j=0;j<D;j++) f<< node_ids[t.element_node_id(i,j)] << " ";
            f << std::endl;
        }

        f << "CELL_TYPES " << t.number_of_elements() << std::endl;
        for (unsigned int i=0;i<t.number_of_elements();i++) {
            f<< ((D==3)?"5":"3") << std::endl;
        }

        f.close();
    }



    template <class TriangulationType>
    void write_triangulated_surface_opendx(const TriangulationType& t, const std::string& filename) {
        //this functions writes a given surface triangulation "t" to the OpenDX (Open Data Explorer)-file "filename"
        //the TriangulationType has to have the same properties as for the levelset-initialization
        //function init (see surface2levelset.hpp)


        const int D=TriangulationType::dimensions;

        std::map<typename TriangulationType::node_index_type,unsigned int> node_ids;

        unsigned int counter=0;
        for (typename TriangulationType::element_index_type e=0;e<t.number_of_elements();e++) {
            for (int k=0;k<D;k++) {
                if (node_ids.find(t.element_node_id(e,k))==node_ids.end()) node_ids[t.element_node_id(e,k)]=counter++;
            }
        }

        node_ids.clear();

        std::ofstream f(filename.c_str());

        //!print positions
        f<< "object \"positions\" class array type float rank 1 shape " << D << " items "<< counter <<" data follows" << std::endl;
        counter=0;
        for (typename TriangulationType::element_index_type e=0;e<t.number_of_elements();e++) {
            for (int k=0;k<D;k++) {
                if (node_ids.find(t.element_node_id(e,k))==node_ids.end()) {
                    for (int j=0;j<D;j++) f << static_cast<float>(t.node_coordinate(t.element_node_id(e,k),j)) << " ";
                    node_ids[t.element_node_id(e,k)]=counter++;
                    f<< std::endl;
                }
            }
        }

        //! print connections
        f << "object \"connections\" class array type int rank 1 shape " << D << " items "<< t.number_of_elements() <<" data follows" << std::endl;
        for (unsigned int i=0;i<t.number_of_elements();i++) {
            for (int j=0;j<D;j++) f<< node_ids[t.element_node_id(i,j)] << " ";
            f << std::endl;
        }

        if (D==2)
            f << "attribute \"element type\" string \"lines\"" << std::endl;
        else if (D==3)
            f << "attribute \"element type\" string \"triangles\"" << std::endl;
        f << "attribute \"ref\" string \"positions\"" << std::endl;

        //! print profile
        f << "object \"profile\" class field" << std::endl;
        f << "  component \"positions\" value \"positions\"" << std::endl;
        f << "  component \"connections\" value \"connections\"" << std::endl;
        f << "end" << std::endl;

        f.close();

    }

    template <class GridTraitsType, class LevelSetTraitsType>
    void write_levelset_opendx(const levelset<GridTraitsType, LevelSetTraitsType>& l, std::string FileName,  bool only_defined_grid_points, float limit, bool only_signs) {
        //this functions writes all defined grid points including their level set values
        //      to file using the OpenDX (Open Data Explorer)-file format
        //limit specifies the range of values assigned to the grid points,
        //      if the level set value of a grid point is
        //      out of this range the value is set to +/- limit
        //if only_defined_grid_points is set to false then also the start and end grid points
        //      of undefined runs are written to file, their values are then set to +/-limit
        //if only_signs is set to true then only the signs of the grid points are written to file,
        //      the grid point values are then set either to +1 or -1

        typedef levelset<GridTraitsType, LevelSetTraitsType> LevelSetType;
//        typedef typename LevelSetType::value_type value_type;
        typedef typename LevelSetType::size_type size_type;
        const int D=LevelSetType::dimensions;


        std::ofstream f(FileName.c_str());

        size_type num=0;

        for (typename LevelSetType::const_iterator_runs it(l);
            !it.is_finished();it.next()) {
            if (only_defined_grid_points) if (!it.is_defined()) continue;
            if (!l.grid().is_at_infinity(it.start_indices())) num++;
            if (it.start_indices()!=it.end_indices()) if (!l.grid().is_at_infinity(it.end_indices())) num++;
        }


        f<< "object 1 class array type float rank 1 shape " << D << " items "<< num <<" data follows" << std::endl;

        for (typename LevelSetType::const_iterator_runs it(l);
                !it.is_finished();it.next()) {
            if (only_defined_grid_points) if (!it.is_defined()) continue;
            if (!l.grid().is_at_infinity(it.start_indices())) {
                for (int j=0;j<D;j++) f << (it.start_indices()[j]) << " ";
                f << std::endl;
            }
            if (!l.grid().is_at_infinity(it.end_indices())) {
                if (it.start_indices()!=it.end_indices()) {
                    for (int j=0;j<D;j++) f << (it.end_indices()[j]) << " ";
                    f << std::endl;
                }
            }

        }

        f << "object 2 class array type float rank 0 items "<< num<<" data follows" <<std::endl;
        for (typename LevelSetType::const_iterator_runs it(l);
                !it.is_finished();it.next()) {
            if (only_defined_grid_points) if (!it.is_defined()) continue;

            float dist;
            if (only_signs) {
                if (it.sign()==POS_SIGN) dist=limit; else dist=-limit;
            } else {
                dist=static_cast<float>(it.value());
                dist=std::min(limit,dist);
                dist=std::max(-limit,dist);
            }
            if (!l.grid().is_at_infinity(it.start_indices())) f << dist << std::endl;
            if (it.start_indices()!=it.end_indices()) if (!l.grid().is_at_infinity(it.end_indices())) f << dist << std::endl;
        }

        f << "attribute \"dep\" string \"positions\"" << std::endl;

        f << "object \"data\" class field" << std::endl;
        f << "component \"positions\" value 1" << std::endl;
        f << "component \"data\" value 2" << std::endl;
        f << "end" << std::endl;

        f.close();
    }


    //used to test the systems endianess, which is needed for levelset file io
    union {
      uint16_t shortVar;    // binary  number of length 16 Bits
      uint8_t  charVar[2];  // 2 binary numbers, each 8 Bits
    } test_endianness;

    bool bigEndian(){
      test_endianness.shortVar = 0x8000; // MSB of 16
      return test_endianness.charVar[0] != 0;
    }

    template <class GridTraitsType> GridTraitsType get_grid_traits_from_lvst_file(std::string path){
      //reads the grid information from a levelset file

      std::ifstream fin(path);
      if(!fin.is_open()) msg::print_error("Could not open the file: " + path);
      char buff[9] = {};
      fin.read(buff, 9);
      //Comparing Identification Bytes
      if(std::string(buff).compare(0, 4, "LvSt")) msg::print_error("File is not a levelset file.");
      if(LVST_FILE_VERSION_NUMBER !=  buff[4]-48) msg::print_error("File version does not match!");
      if(bigEndian() != buff[5]-48) msg::print_warning("File was written in a different byte order than it is being read. Results may be incorrect!");
      const int dim = buff[6]-48;
      const int D = GridTraitsType::dimensions;

      double grid_delta;
      long long grid_min[D] = {};
      long long grid_max[D] = {};
      char b_conditions[D] = {};

      int bytes_per_grid_limit = buff[8]-48;
      //since grid limits can be negative we need a mask to fill the remaining bytes with FF if the sign bit is set (Two's complement)
      unsigned long long sign_bit_mask = 0;
      --sign_bit_mask <<= bytes_per_grid_limit * CHAR_BIT;

      //read in the grid properties,which have variable sizes
      for(int i=dim;i--;){
        fin.read((char *)&grid_min[i], bytes_per_grid_limit);
        //if the sign bit is set, fill the remaining bits with 1s(Two's complement)
        if(grid_min[i] >> (bytes_per_grid_limit * CHAR_BIT-1)) grid_min[i] |= sign_bit_mask;
        fin.read((char *)&grid_max[i], bytes_per_grid_limit);
        if(grid_max[i] >> (bytes_per_grid_limit * CHAR_BIT-1)) grid_max[i] |= sign_bit_mask;
        fin.read((char *)&b_conditions[i], 1);
      }
      fin.read((char *)&grid_delta, sizeof(double));

#ifdef VERBOSE
      msg::print_message("\nReading grid properties from levelset file..." + path);
      std::ostringstream oss;
      oss << "Bytes per grid limit: " << (int)bytes_per_grid_limit << std::endl;
      for(int i=dim;i--;){
        oss << "Dimension " << i << ":" << std::endl
            << "    Grid min: " << grid_min[i] << std::endl
            << "    Grid max: " << grid_max[i] << std::endl
            << "    Boundary condition: " << (int)b_conditions[i] << std::endl;
      }
      oss << "Grid delta: " << grid_delta;
      msg::print_message(oss.str());
#endif
      if(!fin.good()) msg::print_error("Could not read grid properties from file: " + path);
      fin.close();
      return GridTraitsType(grid_min, grid_max, b_conditions, grid_delta);
    }

    template <class GridTraitsType, class LevelSetTraitsType>
    void export_levelset_to_file(levelset<GridTraitsType, LevelSetTraitsType>& ls, const std::string& path, const int& bits_per_distance = 8) {
      /*
      **************************************************************************************************************
      ***************************************    THE LEVELSET FILE FORMAT    ***************************************
      **************************************************************************************************************
      *    File Header: 8 Bytes     *                                                                              *
      ********************************                                                                             *
      *    4 Bytes   Identification Bytes (LvSt)                                                                   *
      *    1 Byte    File Version Number                                                                           *
      *    1 Byte    Endianess - Little Endian (0) or Big Endian (1)                                               *
      *    1 Byte    Dimension of the Levelset (2 or 3)                                                            *
      *    1 Byte    Bits Per Distance                                                                             *
      **************************************************************************************************************
      *    Grid: 15+ Bytes    *                                                                                    *
      *************************                                                                                    *
      *    1 Byte    This byte contains the number of bytes used for the grid min and grid max of each dimension.  *
      *    NOTE: The following block is repeated for each dimension.                                               *
      *    x Bytes   Grid minimum                                                                                  *
      *    x Bytes   Grid maximum                                                                                  *
      *    1 Byte    Boundary Condition                                                                            *
      *    NOTE: END                                                                                               *
      *    8 Bytes   GridDelta (sizeof double)                                                                     *
      **************************************************************************************************************
      *    H-RLE Block Header: 14 Bytes    *                                                                       *
      **************************************                                                                       *
      *    1  Byte   This byte contains the number of bytes used for each start index and runtype                  *
      *    1  Byte   This byte contains the number of bytes used for each runbreak.                                *
      *    4  Bytes  Number of saved Start Indices                                                                 *
      *    4  Bytes  Number of saved Runtypes                                                                      *
      *    4  Bytes  Number of saved Runbreaks                                                                     *
      **************************                                                                                   *
      *    H-RLE Block Data    *                                                                                   *
      **************************                                                                                   *
      *    Start Indices               - using adaptive number of bytes(delta encoded)                             *
      *    Runtypes                    - using 2 bits per runtype (-oo, +oo, defined)                              *
      *    Indices of defined runtypes - using adaptive number of bytes(delta encoded)                             *
      *    Runbreaks                   - using adaptive number of bytes                                            *
      **************************************************************************************************************
      *    Distances Header: 4 Bytes    *                                                                          *
      ***********************************                                                                          *
      *    4 Bytes  Number of distances                                                                            *
      ***********************************                                                                          *
      *    Distances Data               *                                                                          *
      ***********************************                                                                          *
      *    Distances - using 8 bits per distance (default)                                                         *
      *    NOTE: Currently up to 64 bits are supported. When bits_per_distance is 0, 8 bits are used.              *
      **************************************************************************************************************
      */
      std::ofstream fout(path);
      std::ostringstream oss;
      if(path.find(".lvst") == std::string::npos) msg::print_warning("File name does not have the correct file ending.");
      if(!fout.is_open()) msg::print_error("Could not open the file: " + path);

      //type & constant redefinitions
      typedef typename levelset<GridTraitsType, LevelSetTraitsType>::size_type size_type;
      typedef typename levelset<GridTraitsType, LevelSetTraitsType>::index_type index_type;
      typedef typename levelset<GridTraitsType, LevelSetTraitsType>::value_type value_type;
      const unsigned int D = ls.dimensions;

      /************************************************   WRITE FILE HEADER   ************************************************/
      fout << "LvSt" << LVST_FILE_VERSION_NUMBER << (bigEndian() ? 1 : 0)  <<  D;
      fout.write((char *)&bits_per_distance, 1); //bits per distance; as binary
      /************************************************ WRITE GRID PROPERTIES ************************************************/
      //get grid properties(minima, maxima and boundary conditions)
      index_type grid_minima[D];
      index_type grid_maxima[D];

      const boundary_type* grid_b_conditions = ls.grid().grid_traits().boundary_conditions();

      for(unsigned i=0; i<D; ++i){
        if(grid_b_conditions[i] == lvlset::INFINITE_BOUNDARY ||
              grid_b_conditions[i] == lvlset::POS_INFINITE_BOUNDARY){
            grid_maxima[i] = ls.get_max_runbreak(i);
        } else{
          grid_maxima[i] = ls.grid().max_grid_index(i);
        }
        if(grid_b_conditions[i] == lvlset::INFINITE_BOUNDARY ||
              grid_b_conditions[i] == lvlset::NEG_INFINITE_BOUNDARY){
            grid_minima[i] = ls.get_min_runbreak(i);
        } else{
          grid_minima[i] = ls.grid().min_grid_index(i);
        }

        //std::cout << grid_minima[i] << ", " << grid_maxima[i] << std::endl;
      }

      //write grid properties, using adaptive size for min/max
      unsigned int bytes_grid_limits;
      //search for smallest grid_min and highest grid_max
      index_type minimum = 0;
      index_type maximum = 0;
      for(unsigned int i=0; i<D; i++){
        if(minimum > grid_minima[i]) minimum = grid_minima[i];
        if(maximum < grid_maxima[i]) maximum = grid_maxima[i];
      }
      //Use bytesize such that all gridMin and gridMax fit into the number of bytes
      if(minimum >= INT8_MIN && maximum <= INT8_MAX) bytes_grid_limits = 1;
      else if(minimum >= INT16_MIN && maximum <= INT16_MAX) bytes_grid_limits = 2;
      else if(minimum >= INT24_MIN && maximum <= INT24_MAX) bytes_grid_limits = 3;
      else if(minimum >= INT32_MIN && maximum <= INT32_MAX) bytes_grid_limits = 4;
      else if(minimum >= INT40_MIN && maximum <= INT40_MAX) bytes_grid_limits = 5;
      else if(minimum >= INT48_MIN && maximum <= INT48_MAX) bytes_grid_limits = 6;
      else if(minimum >= INT56_MIN && maximum <= INT56_MAX) bytes_grid_limits = 7;
      else  bytes_grid_limits = 8;

#ifdef VERBOSE
      oss << "File Version: " << LVST_FILE_VERSION_NUMBER << std::endl
          << (bigEndian() ? "Big Endian" : "Little Endian")
          << std::endl << "Dimensions: " << D << std::endl
          << "Bits per distance: " << bits_per_distance << std::endl
          << "Bytes per grid min/max: " << bytes_grid_limits << std::endl;
#endif
      fout << bytes_grid_limits;
      for (int dim=D;dim--;) {
        fout.write((char *)&grid_minima[dim], bytes_grid_limits);
        fout.write((char *)&grid_maxima[dim], bytes_grid_limits);
        fout.write((char *)&grid_b_conditions[dim], 1);
#ifdef VERBOSE
        oss << "Dimension " << dim << ":" << std::endl
            << "    Grid min: " << grid_minima[dim] << std::endl
            << "    Grid max: " << grid_maxima[dim] << std::endl
            << "    Max Diff: " << grid_maxima[dim] - grid_minima[dim] << std::endl
            << "    Boundary Condition: " << grid_b_conditions[dim] << std::endl;
#endif
      }
      double delta = ls.grid().grid_traits().grid_position(0, 1);//1st argument is not used. Returns 1 * GridDelta.
      fout.write((char *)&delta, sizeof(double)); //double, aka binary64, has 64 bits(8 bytes) on all platforms as of the IEEE-754 standard(1985; newest revision 2008)
#ifdef VERBOSE
      oss << "Grid delta: " << delta << std::endl
          << "Offset H-RLE block: " << fout.tellp() << std::endl;
      msg::print_message_2(oss.str());
#endif

      int count;
      unsigned char byte = 0;
      // makes parallelized sub_levelsets into one levelset
      ls.serialize();
      /************************************************   WRITE H-RLE BLOCKS  ************************************************/
      for (int dim=D-1;dim>=0;--dim) {
        //get start indices, runbreaks and runtypes
        const std::vector<size_type> & start_indices = ls.start_indices(dim);
        const std::vector<size_type> & runtypes = ls.runtypes(dim);
        const std::vector<index_type> & runbreaks = ls.runbreaks(dim);

        uint8_t bytes_per_index = 1;

        //Delta Encoding. Save only the difference to the next index.
        //The max difference is: gridMax - gridMin
        const unsigned long tmp = grid_maxima[dim] - grid_minima[dim];

        if(tmp <= UINT8_MAX) bytes_per_index = 1;
        else if(tmp <= UINT16_MAX) bytes_per_index = 2;
        else if(tmp <= UINT24_MAX) bytes_per_index = 3;
        else if(tmp <= UINT32_MAX) bytes_per_index = 4;
        else if(tmp <= UINT40_MAX) bytes_per_index = 5;
        else if(tmp <= UINT48_MAX) bytes_per_index = 6;
        else if(tmp <= UINT56_MAX) bytes_per_index = 7;
        else bytes_per_index = 8;

        //bytesize for runbreaks
        //The smallest runbreak is >= grid_min and the largest runbreak is <= grid_max.
        uint8_t bytes_per_runbreak = 1;
        if(grid_minima[dim] >= INT8_MIN && grid_maxima[dim] <= INT8_MAX) bytes_per_runbreak = 1;
        else if(grid_minima[dim] >= INT16_MIN && grid_maxima[dim] <= INT16_MAX) bytes_per_runbreak = 2;
        else if(grid_minima[dim] >= INT24_MIN && grid_maxima[dim] <= INT24_MAX) bytes_per_runbreak = 3;
        else if(grid_minima[dim] >= INT32_MIN && grid_maxima[dim] <= INT32_MAX) bytes_per_runbreak = 4;
        else if(grid_minima[dim] >= INT40_MIN && grid_maxima[dim] <= INT40_MAX) bytes_per_runbreak = 5;
        else if(grid_minima[dim] >= INT48_MIN && grid_maxima[dim] <= INT48_MAX) bytes_per_runbreak = 6;
        else if(grid_minima[dim] >= INT56_MIN && grid_maxima[dim] <= INT56_MAX) bytes_per_runbreak = 7;
        else bytes_per_runbreak = 8;

        //write 14 byte H-RLE block header
        fout.write((char *)&bytes_per_index, 1);
        fout.write((char *)&bytes_per_runbreak, 1);
        uint32_t num = start_indices.size();
        fout.write((char *)&num, 4);
        num = runtypes.size();
        fout.write((char *)&num, 4);
        num = runbreaks.size();
        fout.write((char *)&num, 4);

        uint32_t values_written = 0;
        //Write start indices; only save the difference to the next start index (delta encoding)
        unsigned long diff;
        //First index is always 0, no need to write explicitly
        for(unsigned int i=0; i<start_indices.size()-1; i++){
          diff = start_indices[i+1] - start_indices[i];
          fout.write((char *)&diff, bytes_per_index);
          values_written++;
        }

#ifdef VERBOSE
        oss.str("");
        oss << "Dimension " << dim << ":" << std::endl
            << "    " << (int)bytes_per_index << " byte(s) per start index or runtype." << std::endl
            << "    " << (int)bytes_per_runbreak << " byte(s) per runbreak." << std::endl
            << "    " << (values_written+1) << " of " << start_indices.size() << " start indices written.";
        msg::print_message_2(oss.str());
#endif

        //write all runtypes to the file, skipping all segments and indices (using 2 bits per runtype)
        count = CHAR_BIT/2 - 1;
        byte = 0;
        values_written = 0;
        std::vector<size_type> def_run_indices = {}; //store all indices of defined runtypes
        for (typename std::vector<size_type>::const_iterator it=runtypes.begin();it!=runtypes.end();++it) {
          if(*it == ls.POS_PT) //01 - positive undefined runtype
            byte |= 1 << count * 2;
          else if(*it == ls.NEG_PT) //11 - negative undefined runtype
            byte |= 3 << count * 2;
          else if(*it == ls.UNDEF_PT) //10 - uninitialized runtype
            byte |= 2 << count * 2;
          else if(!ls.is_defined(*it)){ // skip Segments --> there shouldn't be any because of serialize()
            msg::print_warning("Segment detected during levelset file export.");
            continue;
          }
          else { //00 - defined runtype
              byte |= 0 << count * 2;
              def_run_indices.push_back(*it);
          }
          values_written++;
          count--;
          if(count < 0){ //if the byte contains 4 runtypes, write it
            fout << byte;
            count = CHAR_BIT/2 - 1;
            byte = 0;
          }
        }
        if(count >= 0 && count < CHAR_BIT/2 -1)//if the last byte contains less than 4 runtyes, write it
            fout << byte;

        //Write indices of defined runtypes; only save the difference to the next defined runtype
        //write the first runtype(always 0) explicitly; makes reading easier
        fout.write((char *)&def_run_indices[0], bytes_per_index);
        for(unsigned int i=0; i<def_run_indices.size()-1; i++){
          diff = def_run_indices[i+1] - def_run_indices[i];
          fout.write((char *)&diff, bytes_per_index);
        }

#ifdef VERBOSE
      oss.str("");
      oss << "    " << values_written << " of " << runtypes.size() << " runtypes written. Defined runtypes: " << def_run_indices.size();
      msg::print_message_2(oss.str());
#endif
        //Write runbreaks
        values_written = 0;
        for (typename std::vector<index_type>::const_iterator it=runbreaks.begin();it!=runbreaks.end();++it) {
          fout.write((char *)&(*it), bytes_per_runbreak);
          values_written++;
        }
#ifdef VERBOSE
      oss.str("");
      oss << "    " << values_written << " of " << runbreaks.size() << " runbreaks written.";
      msg::print_message_2(oss.str());
#endif
      }
      /************************************************    WRITE  DISTANCES    ************************************************/
      const std::vector<value_type> & distances = ls.distances();
      uint32_t num = distances.size(), values_written = 0;
      fout.write((char *)&num, 4);

      count = std::ceil((double)CHAR_BIT/bits_per_distance) -1;
      const int num_overflow_bytes = bits_per_distance/CHAR_BIT -1; //-1 because bpd is mapped to the next higher power of 2; 16/8 = 2 but overflow bytes is only 1
      const int bits_per_byte = num_overflow_bytes > 0 ? CHAR_BIT : bits_per_distance;
      const long double value = std::pow(2.0L, (long double) bits_per_distance-1)-0.5L;

#ifdef VERBOSE
      oss.str("");
      oss << "Offset Distances: " << fout.tellp() << std::endl
          << "Value: " << value << std::endl
          << "Distances per byte: " << count << std::endl
          << "Number of overflow bytes: " << num_overflow_bytes << std::endl;
      msg::print_message_2(oss.str());
#endif
      /*
      *************************************************************************************
      *  NOTE: In the current implementation x = 1                                        *
      *  Levelset values range from -x .... +x, where x = levelset.num_layers/2           *
      *  With n bits we can represent values from 0 .... (2^n -1)                         *
      *         -x .... +x            |* (2^n-1)/(2x)                                     *
      * -(2^n-1)/2 .... +(2^n-1)/2    |+ (2^n-1)/2                                        *
      *          0 .... +(2^n-1)                                                          *
      *************************************************************************************
      */
      unsigned long long discrete_distance, overflow;
      long double tmp;
      byte = 0;
      for (typename std::vector<value_type>::const_iterator it=distances.begin();it!=distances.end();++it) {
        tmp = *it * value + value + 0.5; //+0.5 for rounding, std::llround returns a signed long long and not an unsigned one
        discrete_distance = (unsigned long long)tmp;
        overflow = discrete_distance >> bits_per_byte;
        discrete_distance  <<= count * bits_per_byte;
        byte |= discrete_distance;
        count--;
        if(count < 0){
          fout << byte;
          count = CHAR_BIT/bits_per_distance-1;
          if(count < 0) count = 0;
          //write overflow, even if it is 0, so each distance has the same amount of bytes
          for(int i=0; i<num_overflow_bytes; i++){
            byte = 0;
            byte |= overflow;
            fout << byte;
            overflow >>= bits_per_byte;
          }
          byte=0;
        }
        values_written++;
      }
      if(count >= 0 && count < CHAR_BIT/bits_per_distance-1 )
          fout << byte;

#ifdef VERBOSE
      oss.str("");
      oss << values_written << " of " << distances.size() << " distances written.";
      msg::print_message_2(oss.str());
#endif

      if(!fout.good()) msg::print_error("Writing to file " + path + " failed.");
      fout.close();

      ls.finalize(2);//sets the segmentation points
      //ls.prune() is called before writing the levelset, therefore we only need to set up the segmentation
      ls.segment();//parallelize the levelset
    }

    template <class GridTraitsType, class LevelSetTraitsType>
    void import_levelset_from_file(levelset<GridTraitsType, LevelSetTraitsType>& ls, const std::string& path){
      //this function takes an empty levelset and reads it from a file
      std::ifstream fin(path);
      std::ostringstream oss;
      if(path.find(".lvst") == std::string::npos)
        msg::print_warning("File name does not have the correct file ending.");
      if(!fin.is_open()) msg::print_error("Could not open the file: " + path);

      typedef typename levelset<GridTraitsType, LevelSetTraitsType>::size_type size_type;
      typedef typename levelset<GridTraitsType, LevelSetTraitsType>::index_type index_type;
      typedef typename levelset<GridTraitsType, LevelSetTraitsType>::value_type value_type;

      /***********************************************  READ FILE HEADER  ***********************************************/
      char buff[9] = {};
      unsigned char byte;
      fin.read(buff, 9);
      //Comparing Identification Bytes
      if(std::string(buff).compare(0, 4, "LvSt")) msg::print_error("File is not a levelset file.");
      if(LVST_FILE_VERSION_NUMBER !=  buff[4]-48) msg::print_warning("File version does not match!");
      if(bigEndian() != buff[5]-48) msg::print_warning("File was written in a different byte order than it is being read. Results may be incorrect!");

      const int dim = buff[6]-48;
      int bits_per_distance = buff[7];
      const int bytes_per_grid_limits = buff[8]-48;
      /*
      **************************** SKIP GRID PROPERTIES(they have to be read before the levelset) ******************************
      *  Current position:  fin.tellg()                                                                                        *
      *  For each dimension add: 2 * bytes_per_grid_limits + 1 byte (per boundary condition)                                   *
      *  Add: the additional number of bytes for the grid delta; sizeof(double)                                                *
      **************************************************************************************************************************
      */
      fin.seekg(int(fin.tellg()) + (bytes_per_grid_limits*2+1)*dim + sizeof(double));
#ifdef VERBOSE
      oss << std::endl << "Dimensions: " << dim << std::endl
          << "Bits per distance: " << bits_per_distance << std::endl
          << "Bytes per grid min/max: " << bytes_per_grid_limits << std::endl
          << "Offset H-RLE block: " << fin.tellg();
      msg::print_message_2(oss.str());
#endif

      //initialize the empty levelset to sub_levelsets size 1 (no segmentation)
      ls.initialize();
      /*********************************************** READ H-RLE BLOCKS ***********************************************/
      for(int i=dim;i--;){
        //get the start indices, runtypes and runbreaks vectors
        std::vector<size_type>& start_indices = ls.start_indices(i);
        std::vector<size_type>& runtypes = ls.runtypes(i);
        std::vector<index_type>& runbreaks = ls.runbreaks(i);
        uint32_t num_start_indices, num_runtypes, num_runbreaks;
        uint8_t bytes_per_index, bytes_per_runbreak;
        //reading in the 14 byte H-RLE header
        fin.read((char *)&bytes_per_index, 1);
        fin.read((char *)&bytes_per_runbreak, 1);
        fin.read((char *)&num_start_indices, 4);
        fin.read((char *)&num_runtypes, 4);
        fin.read((char *)&num_runbreaks, 4);

        uint32_t values_read = 0;
        if(start_indices.size() > 0) start_indices.clear();
        //read differences between start indices
        unsigned long long sum = 0;
        unsigned long long current = 0;
        start_indices.push_back(current);//push the 0, it was not written to the file
        for(unsigned int i=0; i<num_start_indices-1; i++){//-1 cause we didnt write the 0
          current = 0;
          fin.read((char *)&current, bytes_per_index);
          sum += current;
          start_indices.push_back(sum);
          values_read++;
        }

#ifdef VERBOSE
        oss.str("");
        oss << "Dimension " << i << ":";
        msg::print_message_2(oss.str());
        oss.str("");
        oss << "    " << (int)bytes_per_index << " byte(s) per start index and runtype." << std::endl
            << "    " << (int)bytes_per_runbreak << " byte(s) per runbreak." << std::endl
            << "    " << (values_read+1) << " of " << num_start_indices << " start indices read.";//first index (0) is not written to file
        msg::print_message_2(oss.str());
#endif
        uint32_t count = 0;
        values_read = 0;
        //reading runtypes
        if(runtypes.size() > 0) runtypes.clear();
        uint32_t bytes_to_read = std::ceil(num_runtypes/4.0);
        uint32_t runtypes_per_byte = CHAR_BIT/2;
        unsigned long long uInt;
        std::ifstream tmp_fin(path);//temporary stream to read defined runtypes from
        if(!tmp_fin.is_open()) msg::print_error("Could not open the tmp file: " + path);
        tmp_fin.seekg(fin.tellg() + (long)bytes_to_read);//set position to the defined indices for runtypes
        sum = 0;
        for(unsigned int y = 0; y < bytes_to_read; y++){
          fin.read((char *)&byte, 1);
          for(unsigned int z = runtypes_per_byte; z--;){
            if(values_read == num_runtypes) break;
            uInt = byte >> z*2 & 0x3;
            if(uInt == 1) runtypes.push_back(ls.POS_PT);
            else if(uInt == 3) runtypes.push_back(ls.NEG_PT);
            else if(uInt == 2) runtypes.push_back(ls.UNDEF_PT);
            else if(uInt == 0) {
                current = 0;
                tmp_fin.read((char *)&current, bytes_per_index);//read index of defined runtype
                sum += current;
                runtypes.push_back(sum);
                count++;
            }
            values_read++;
          }
        }
        fin.seekg(tmp_fin.tellg());//set position of stream to after the defined runtypes
        tmp_fin.close();
#ifdef VERBOSE
        oss.str("");
        oss << "    " << values_read << " of " << num_runtypes << " runtypes read. Defined runtypes: " << count;
        msg::print_message_2(oss.str());
#endif

        values_read = 0;
        //reading runbreaks
        long long sInt;
        /*
        *********************************************************************************************************************************
        * Example: We have 2 bytes and the runbreak is -7, which corresponds to FFF9 (Two's complement).                                *
        * Now if we read it into a 4 byte integer we will have 0000 FFF9, which is not interpreted as -7, but as a positive number.     *
        * -7 as a 4 byte integer is: FFFF FFF9.                                                                                         *
        * We need a mask to set the higher bytes to FF if the sign bit is set                                                           *
        *********************************************************************************************************************************
        */
        unsigned long long sign_bit_mask = 0;
        --sign_bit_mask <<= bytes_per_runbreak * CHAR_BIT;

        if(runbreaks.size() > 0) runbreaks.clear();
        for(unsigned int z = 0; z < num_runbreaks; z++){
          sInt = 0;
          fin.read((char *) &sInt, bytes_per_runbreak);
          //if the sign bit is set, fill up the higher bytes
          if(sInt >> (bytes_per_runbreak * CHAR_BIT -1)) sInt |= sign_bit_mask;
          runbreaks.push_back(sInt);
          values_read++;
        }
#ifdef VERBOSE
        oss.str("");
        oss << "    " << values_read << " of " << num_runbreaks << " runbreaks." << std::endl;
        msg::print_message_2(oss.str());
#endif
      }
      /*********************************************** READ DISTANCES ***********************************************/
      uint32_t num_distances = 0, values_read = 0;
      //reading the number of distances to read
      fin.read((char *)&num_distances, 4);
#ifdef VERBOSE
      oss.str("");
      oss << "Offset distances: " << fin.tellg();
      msg::print_message_2(oss.str());
#endif
      std::vector<value_type> & distances = ls.distances();
      const long double value = 1.0L / (std::pow(2.0L, (long double) bits_per_distance-1)-0.5L);//value ^ -1 to avoid division
      const int count = std::ceil((double)CHAR_BIT/bits_per_distance);
      const int num_overflow_bytes = bits_per_distance/CHAR_BIT -1; //-1 because bpd is mapped to the next higher power of 2
      const int bits_per_byte = num_overflow_bytes > 0 ? CHAR_BIT : bits_per_distance;
      const int num_reads = std::ceil((double)num_distances/count);
      const unsigned char mask = 0xFF >> (bits_per_distance < CHAR_BIT ? CHAR_BIT - bits_per_distance : 0); //for bits per distance > 8 the mask will be 0xFF

#ifdef VERBOSE
      oss.str("");
      oss << "Value: " << value << std::endl;
      oss << "Distance(s) per byte: " << count << std::endl;
      oss << "Number of overflow bytes: " << num_overflow_bytes << std::endl;
      msg::print_message_2(oss.str());
#endif

      //reading distances
      if(distances.size() > 0) distances.clear();
      unsigned long long discrete_distance = 0;
      unsigned long long tmp_distance = 0;
      for(int i = 0; i < num_reads; i++){
        fin.read((char *)&byte, 1);
        for(int z = count; z--;){
          if(values_read == num_distances) break; //if distances are odd, skip padding bits
          discrete_distance = 0;
          discrete_distance |= byte >> z*bits_per_byte & mask;
          //read in overflow
          for(int j=0; j<num_overflow_bytes; j++){
            fin.read((char *)&byte, 1);
            tmp_distance = byte;
            //shift the byte to the correct position
            discrete_distance |= tmp_distance << (j+1) * bits_per_byte;
          }
          distances.push_back(discrete_distance * value - 1.0);

          values_read++;
        }
      }
#ifdef VERBOSE
      oss.str("");
      oss << values_read << " of " << num_distances << " distances read.";
      msg::print_message_2(oss.str());
#endif
      if(!fin.good()) msg::print_error("Reading from " + path + " failed.");
      fin.close();

      ls.finalize(2);
      //ls.prune() is called before writing the levelset file, therefore we only need to set up segmentation.
      ls.segment();
      ls.set_levelset_id();
    }

} //NAMESPACE LVLSET END

#endif /*OUTPUT_HPP_*/
