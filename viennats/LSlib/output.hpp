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
#include <float.h>
#include <cstdint>
#include "kernel.hpp"
#include "levelset2surface.hpp"
#include "../message.h"

//Options for levelset output
#define LVST_FILE_VERSION_NUMBER 1

//Limits for 3,5,6,7 byte output; only change if one byte does not have 8 bits
#define UINT24_MAX 16777215L // highest value of 3 unsigned bytes  2^24-1
#define  INT24_MAX  8388607L //  highest value of 3 signed   bytes  2^23-1
#define  INT24_MIN -8388608L //  lowest  value of 3 signed   bytes -2^23

#define UINT40_MAX 1099511627775L // highest value of 3 unsigned bytes  2^40-1
#define  INT40_MAX  549755813887L //  highest value of 3 signed   bytes  2^39-1
#define  INT40_MIN -549755813888L //  lowest  value of 3 signed   bytes -2^39

#define UINT48_MAX 281474976710655L  // highest value of 3 unsigned bytes  2^48-1
#define  INT48_MAX  140737488355327L //  highest value of 3 signed   bytes  2^47-1
#define  INT48_MIN -140737488355328L //  lowest  value of 3 signed   bytes -2^47

#define UINT56_MAX 72057594037927935L  // highest value of 3 unsigned bytes  2^56-1
#define  INT56_MAX  36028797018963967L //  highest value of 3 signed   bytes  2^55-1
#define  INT56_MIN -36028797018963968L //  lowest  value of 3 signed   bytes -2^55

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

    template <class GridTraitsType> GridTraitsType getGridFromLVSTFile(std::string path){
      //reads the grid information from a levelset file

      std::ifstream fin(path);
      if(!fin.is_open()) msg::print_error("ERROR: Couldn't open the file: " + path);
      char buff[10] = {};
      fin.read(buff, 10);
      const int dim = buff[4]-48;
      const int D = GridTraitsType::dimensions;
      //Comparing Identification Bytes
      if(std::string(buff).compare(0, 4, "LvSt")) msg::print_error("File is not a levelset file.");
      if(LVST_FILE_VERSION_NUMBER !=  buff[5]-48) msg::print_warning("File version does not match!");
      if(bigEndian() != buff[6]-48) msg::print_warning("File was written in a different byte order than it is being read. Results may be incorrect!");

      double gridDelta;
      int32_t grid_min[D] = {};
      int32_t grid_max[D] = {};
      char bconditions[D] = {};


      uint8_t bytes_per_grid_limit = buff[9]>>4;

      uint32_t sign_bit_mask = 0;
      --sign_bit_mask <<= bytes_per_grid_limit * CHAR_BIT;

      //read in the grid properties,which have variable sizes
      for(int i=dim;i--;){
        fin.read((char *)&grid_min[i], bytes_per_grid_limit);
        if(grid_min[i] >> (bytes_per_grid_limit * CHAR_BIT-1) & 0x1 ) grid_min[i] |= sign_bit_mask;
        fin.read((char *)&grid_max[i], bytes_per_grid_limit);
        if(grid_max[i] >> (bytes_per_grid_limit * CHAR_BIT-1) & 0x1 ) grid_max[i] |= sign_bit_mask;
        fin.read((char *)&bconditions[i], 1);
      }
      fin.read((char *)&gridDelta, buff[9] & 0xF);

#ifdef VERBOSE
  msg::print_message("Reading grid properties from levelset file..." + path);
  std::ostringstream oss;
  oss << "Bytes per grid limit: " << (char)bytes_per_grid_limit << std::endl;
  msg::print_message_2(oss.str());
  for(int i=dim;i--;){
    oss.str("");
    oss << "Dim " << i << ":";
    msg::print_message_2(oss.str());
    oss.str("");
    oss << "\tGrid min: " << grid_min[i];
    msg::print_message_2(oss.str());
    oss.str("");
    oss << "\tGrid max: " << grid_max[i];
    msg::print_message_2(oss.str());
    oss.str("");
    oss << "\tBoundary condition: " << (int)bconditions[i];
    msg::print_message_2(oss.str());
  }
  oss.str("");
  oss << "Grid delta: " << gridDelta;
  msg::print_message(oss.str());
#endif
      if(fin.fail()) msg::print_error("Couldn't read grid properties from file: " + path);
      fin.close();
      return GridTraitsType(grid_min, grid_max, bconditions, gridDelta);
    }

    template <class GridTraitsType, class LevelSetTraitsType>
    void exportLevelsetToFile(levelset<GridTraitsType, LevelSetTraitsType>& ls, const std::string& path, const int& bits_per_distance = 4) {
        /*************************************************************************************************************
        ***************************************    THE LEVELSET FILE FORMAT    ***************************************
        **************************************************************************************************************
        *    File Header: 33 - 46 Bytes    *                                                                         *
        ************************************                                                                         *
        *     4 Bytes   Identification Bytes (LvSt)                                                                 *
        *     1 Byte    Dimension of the Levelset (2 or 3)                                                           *
        *     1 Byte    File Version Number                                                                          *
        *     1 Byte    Endianess - Little Endian (0) or Big Endian (1)                                              *
        *     1 Byte    Bits Per Distance                                                                            *
        *     1 Byte    Bits per Byte (char)                                                                         *
        *     1 Byte    This byte contains the number of bytes used for the grid min and grid max of each dimension. *
        *               In addition it also contains the number of bytes used for the grid delta.                    *
        *               The maximum for each is 16, since we use 4 bits(0 ... 15) we have to subtract 1.             *
        *  For example: 0x37 means that each grid min/max takes 4 bytes each and the grid delta 8 bytes.             *
        *   NOTE: The following block is repeated for each dimension.                                                *
        *     x Bytes   Grid minimum                                                                                 *
        *     x Bytes   Grid maximum                                                                                 *
        *     1 Byte    Boundary Condition                                                                           *
        *   NOTE: END                                                                                                *
        *   4-8 Bytes   GridDelta (sizeof double)                                                                    *
        **************************************************************************************************************
        *    H - RLE Block Header: 15 Bytes    *                                                                     *
        ****************************************                                                                     *
        *    1  Byte   This byte contains the number of bytes used for each start index.                             *
        *    1  Byte   This byte contains the number of bytes used for each runtype.                                 *
        *    1  Byte   This byte contains the number of bytes used for each runbreak.                                *
        *    4  Bytes  Number of saved Start Indices                                                                 *
        *    4  Bytes  Number of saved Runtypes                                                                      *
        *    4  Bytes  Number of saved Runbreaks                                                                     *
        **************************************************************************************************************
        *    H - RLE Block Data    *                                                                                 *
        ****************************                                                                                 *
        *    Start Indices               - using adaptive number of bytes (default)                                  *
        *    Runtypes                    - using 2 bits per runtype(-oo, +oo, defined)                                                  *
        *    Indices of defined runtypes - using adaptive number of bytes (default)                                  *
        *    Runbreaks                   - using adaptive number of bytes (default)                                  *
        **************************************************************************************************************
        *    Distances Header: 4 Bytes    *                                                                          *
        ***********************************                                                                          *
        *    4 Bytes  Number of distances                                                                            *
        **************************************************************************************************************
        *    Distances - using 4 bits per distance (default)                                                         *
        *    //NOTE: Currently up to 8 bits are supported.                                                           *
        **************************************************************************************************************/
        std::ofstream fout(path);
        std::ostringstream oss;
        if(!fout.is_open()) {msg::print_error("Couldn't open the file: " + path); return;}

        //type & constant redefinitions
        typedef typename levelset<GridTraitsType, LevelSetTraitsType>::size_type size_type;
        typedef typename levelset<GridTraitsType, LevelSetTraitsType>::index_type index_type;
        typedef typename levelset<GridTraitsType, LevelSetTraitsType>::value_type value_type;
        const unsigned int D = ls.dimensions;

        /************************************************ WRITE THE FILE HEADER ************************************************/
        fout << "LvSt" <<  D << LVST_FILE_VERSION_NUMBER << (bigEndian() ? 1 : 0);
        fout.write((char *)&bits_per_distance, 1); //bits per distance
        fout << CHAR_BIT; //bits per char(byte)
        //get grid properties
        const index_type* gridMinima = ls.grid().grid_traits().getMinima();
        const index_type* gridMaxima = ls.grid().grid_traits().getMaxima();
        const boundary_type* gridBConditions = ls.grid().grid_traits().getBoundaryConditions();
        unsigned char byte = 0;
        //write grid properties, using adaptive size for min/max
        unsigned int bytesGridLimits;
        //search for smallest gridMin and highest gridMax
        index_type minGridMin = 0;
        index_type maxGridMax = 0;
        for(unsigned int i=0; i<D; i++){
          if(minGridMin > gridMinima[i]) minGridMin = gridMinima[i];
          if(maxGridMax < gridMaxima[i]) maxGridMax = gridMaxima[i];
        }
        //Use bytesize such that all gridMin and gridMax fit into the number of bytes
        if(minGridMin >= INT8_MIN && maxGridMax <= INT8_MAX) bytesGridLimits = 1;
        else if(minGridMin >= INT16_MIN && maxGridMax <= INT16_MAX) bytesGridLimits = 2;
        else if(minGridMin >= INT24_MIN && maxGridMax <= INT24_MAX) bytesGridLimits = 3;
        else if(minGridMin >= INT32_MIN && maxGridMax <= INT32_MAX) bytesGridLimits = 4;
        else if(minGridMin >= INT40_MIN && maxGridMax <= INT40_MAX) bytesGridLimits = 5;
        else if(minGridMin >= INT48_MIN && maxGridMax <= INT48_MAX) bytesGridLimits = 6;
        else if(minGridMin >= INT56_MIN && maxGridMax <= INT56_MAX) bytesGridLimits = 7;
        else  bytesGridLimits = 8;
#ifdef VERBOSE
        oss << std::endl << "Dimensions: " << D << std::endl
            << "Bits per distance: " << bits_per_distance << std::endl
            << "Char bit: " << CHAR_BIT << std::endl
            << "Bytes per grid min/max: " << bytesGridLimits << std::endl
            << "Bytes per grid delta: " << sizeof(double);
        msg::print_message_2(oss.str());
#endif
        byte |= bytesGridLimits << 4;//bytes used to store grid min and grid max; boundary condition is always 1 byte
        byte |= sizeof(double);//bytes used to store grid delta
        fout << byte;
        for (int dim=D;dim--;) {
            fout.write((char *)&gridMinima[dim], bytesGridLimits);
            fout.write((char *)&gridMaxima[dim], bytesGridLimits);
            fout.write((char *)&gridBConditions[dim], 1);
        }
        double delta = ls.grid().grid_traits().grid_position(0, 1);//1st argument is not used. Returns 1 * GridDelta.
        fout.write((char *)&delta, sizeof(double));
        /************************************************ WRITE FILE HEADER END ************************************************/
#ifdef VERBOSE
        oss.str("");
        oss << "Offset FILE HEADER: " << fout.tellp() << std::endl;
        msg::print_message_2(oss.str());
#endif
        int count;
        // makes parallelized sub_levelsets into one levelset
        ls.serialize();
        /************************************************ WRITE HRLE BLOCKS ************************************************/
        for (int dim=D-1;dim>=0;--dim) {
          //get start indices, runbreaks and runtypes
          //NOTE: the functions return non-const references
          const std::vector<size_type> & startIndices = ls.startIndices(dim);
          const std::vector<size_type> & runTypes = ls.runTypes(dim);
          const std::vector<index_type> & runBreaks = ls.runBreaks(dim);

          //Getting the largest possible values for a start index, runtype and runbreak
          //This is done to use the minimum amount of bytes per start index, runtype and runbreak per dimension.

          //bytesize for start indices
          //The largest start index is <= the size of runtypes. (start index points to a runtype)
          uint8_t bytesPerStIndex = 0;
          if(runTypes.size() <= UINT8_MAX) bytesPerStIndex = 1;
          else if(runTypes.size() <= UINT16_MAX) bytesPerStIndex = 2;
          else if(runTypes.size() <= UINT24_MAX) bytesPerStIndex = 3;
          else if(runTypes.size() <= UINT32_MAX) bytesPerStIndex = 4;
          else if(runTypes.size() <= UINT40_MAX) bytesPerStIndex = 5;
          else if(runTypes.size() <= UINT48_MAX) bytesPerStIndex = 6;
          else if(runTypes.size() <= UINT56_MAX) bytesPerStIndex = 7;
          else bytesPerStIndex = 8;

          //bytesize for runtypes
          //The largest runtype is <= the size of start indices of the next dimension.
          //If dim == 0 the largest runtype <= the size of distances. (defined runtype points to a start index/distance)
          uint8_t bytesPerRnType = 0;
          //runtypes point to start indices for higher dimensions, and to distances for the lowest dimension (0)
          unsigned long tmp = (dim == 0) ? ls.distances().size() : ls.startIndices(dim-1).size();
          if(tmp <= UINT8_MAX) bytesPerRnType = 1;
          else if(tmp <= UINT16_MAX) bytesPerRnType = 2;
          else if(tmp <= UINT24_MAX) bytesPerRnType = 3;
          else if(tmp <= UINT32_MAX) bytesPerRnType = 4;
          else if(tmp <= UINT40_MAX) bytesPerRnType = 5;
          else if(tmp <= UINT48_MAX) bytesPerRnType = 6;
          else if(tmp <= UINT56_MAX) bytesPerRnType = 7;
          else bytesPerRnType = 8;

          //bytesize for runbreaks
          //The smallest runbreak is >= grid_min.
          //The largest runbreak is <= grid_max.
          uint8_t bytesPerRnBreak = 0;
          if(gridMinima[dim] >= INT8_MIN && gridMaxima[dim] <= INT8_MAX) bytesPerRnBreak = 1;
          else if(gridMinima[dim] >= INT16_MIN && gridMaxima[dim] <= INT16_MAX) bytesPerRnBreak = 2;
          else if(gridMinima[dim] >= INT24_MIN && gridMaxima[dim] <= INT24_MAX) bytesPerRnBreak = 3;
          else if(gridMinima[dim] >= INT32_MIN && gridMaxima[dim] <= INT32_MAX) bytesPerRnBreak = 4;
          else if(gridMinima[dim] >= INT40_MIN && gridMaxima[dim] <= INT40_MAX) bytesPerRnBreak = 5;
          else if(gridMinima[dim] >= INT48_MIN && gridMaxima[dim] <= INT48_MAX) bytesPerRnBreak = 6;
          else if(gridMinima[dim] >= INT56_MIN && gridMaxima[dim] <= INT56_MAX) bytesPerRnBreak = 7;
          else bytesPerRnBreak = 8;

          if(!bytesPerStIndex || !bytesPerRnType || !bytesPerRnBreak) msg::print_warning("Bytesize for levelset file is 0!");
#ifdef VERBOSE
          oss.str("");
          oss << int(bytesPerStIndex) << " byte(s) per start index." << std::endl
              << int(bytesPerRnType) << " byte(s) per runtype." << std::endl
              << int(bytesPerRnBreak) << " byte(s) per runbreak.";
          msg::print_message_2(oss.str());
#endif
          //write 15 byte H-RLE block header
          fout.write((char *)&bytesPerStIndex, 1);
          fout.write((char *)&bytesPerRnType, 1);
          fout.write((char *)&bytesPerRnBreak, 1);

          uint32_t num = startIndices.size();
          fout.write((char *)&num, 4);
          num = runTypes.size();
          fout.write((char *)&num, 4);
          num = runBreaks.size();
          fout.write((char *)&num, 4);

          uint32_t values_written = 0;
          //write the start indices to the file
          for (typename std::vector<size_type>::const_iterator it=startIndices.begin();it!=startIndices.end();++it) {
              fout.write((char *)&(*it), bytesPerStIndex);
              values_written++;
          }

#ifdef VERBOSE
          oss.str("");
          oss << "Dimension " << dim;
          msg::print_message_2(oss.str());
          oss.str("");
          oss << "\t" << values_written << " of " << startIndices.size() << " start indices written.";
          msg::print_message_2(oss.str());
#endif

          //write all runtypes to the file, skipping all segments and indices (using 2 bits per runtype)
          count = CHAR_BIT/2 - 1;
          byte=0;
          values_written = 0;
          std::vector<size_type> def_run_indices = {}; //store all indices for defined runtypes
          for (typename std::vector<size_type>::const_iterator it=runTypes.begin();it!=runTypes.end();++it) {
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
            if(count < 0){ //if 4 runtypes are written into the byte, write it to the file
              fout << byte;
              count = CHAR_BIT/2 - 1;
              byte=0;
            }
          }
          if(count >= 0 && count < CHAR_BIT/2 -1)//if byte is not yet full, but the loop finished
              fout << byte;

#ifdef VERBOSE
        oss.str("");
        oss << "\t" << values_written << " of " << runTypes.size() << " runtypes written." << " Defined runtypes: " << def_run_indices.size();
        msg::print_message_2(oss.str());
#endif

          //write defined runtypes
          for (typename std::vector<size_type>::const_iterator it=def_run_indices.begin();it!=def_run_indices.end();++it) {
            fout.write((char *)&(*it), bytesPerRnType);
          }

          //Write runbreaks
          values_written = 0;
          for (typename std::vector<index_type>::const_iterator it=runBreaks.begin();it!=runBreaks.end();++it) {
            fout.write((char *)&(*it), bytesPerRnBreak);
            values_written++;
          }
#ifdef VERBOSE
        oss.str("");
        oss << "\t" << values_written << " of " << runBreaks.size() << " runbreaks written.";
        msg::print_message_2(oss.str());
#endif
        }
        /************************************************ WRITE HRLE BLOCKS END ************************************************/
#ifdef VERBOSE
        oss.str("");
        oss << "Offset Distances: " << fout.tellp() << std::endl;
        msg::print_message_2(oss.str());
#endif
        /************************************************ WRITE DISTANCES ************************************************/
        const std::vector<value_type> & distances = ls.distances();
        uint32_t num = distances.size(), values_written = 0;
        fout.write((char *)&num, 4);
        count = CHAR_BIT/bits_per_distance -1;
        if(count < 0) count = 0;
        int overflow_num_bytes = bits_per_distance/CHAR_BIT -1; //-1 because bpd is mapped to the next higher power of 2; 16/8 = 2 but overflow bytes is only 1
        int bits_per_byte = overflow_num_bytes > 0 ? CHAR_BIT : bits_per_distance; //for 8 bits 1 ... 8
        byte = 0;
        long double value = std::pow((long double) 2, (long double) bits_per_distance-1)-0.5L;
        /*switch (bits_per_distance) {
          case 16:
            value = UINT16_MAX/2.0L;
            break;
          case 24:
            value = UINT24_MAX/2.0L;
            break;
          case 32:
            value = UINT32_MAX/2.0L;
            break;
          case 40:
            value = UINT40_MAX/2.0L;
            break;
          case 48:
            value = UINT48_MAX/2.0L;
            break;
          case 56:
            value = UINT56_MAX/2.0L;
            break;
          case 64:
            value = UINT64_MAX/2.0L;
            break;
          default:
            value = std::pow((long double) 2, (long double) bits_per_distance - 1) - 0.5L;
            break;
        }*/
#ifdef VERBOSE
        oss.str("");
        oss << "Value: " << value << std::endl;
        oss << "Number of overflow bytes: : " << overflow_num_bytes << std::endl;
        msg::print_message_2(oss.str());
#endif
        /************************************************************************************
        *  NOTE: In the current implementation x = 1                                        *
        *  Levelset values range from -x .... +x, where x = levelset.num_layers/2           *
        *  With n bits we can represent values from 0 .... (2^n -1)                         *
        *         -x .... +x            |* (2^n-1)/(2x)                                     *
        * -(2^n-1)/2 .... +(2^n-1)/2    |+ (2^n-1)/2                                        *
        *          0 .... +(2^n-1)                                                          *
        *************************************************************************************/
        for (typename std::vector<value_type>::const_iterator it=distances.begin();it!=distances.end();++it) {
          unsigned long long discrete_distance = std::llround(*it * value + value);
          unsigned long long overflow = discrete_distance >> bits_per_byte;
          byte |= discrete_distance  << count * bits_per_byte;
          count--;
          if(count < 0){
            fout << byte;
            count = CHAR_BIT/bits_per_distance-1;
            if(count < 0) count = 0;
            //write overflow, even if it is 0, so each distance has the same amount of bytes
            for(int i=0; i<overflow_num_bytes; i++){
              byte = 0;
              byte |= overflow;
              fout << byte;
              overflow = overflow >> bits_per_byte;
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
        if(!fout.good()) {
          msg::print_message_2("Error while writing to file " + path + ".");
          if(fout.bad()){
            fout.close();
            msg::print_error("Write error on output operation.");
          }
          else if(fout.fail()){
            fout.close();
            msg::print_error("Logical error on output operation.");
          }
          else{
            fout.close();
            msg::print_error("Undetermined error on output operation.");
          }
        }
        fout.close();

        ls.finalize(2);
        //ls.prune();
        //prune() is called before writing the levelset, therefore we only need to set up the segmentation
        ls.print();
        ls.segment();
    }

    template <class GridTraitsType, class LevelSetTraitsType>
    void importLevelsetFromFile(levelset<GridTraitsType, LevelSetTraitsType>& ls, const std::string& path){
      //this function takes an empty levelset and reads from a file
      std::ifstream fin(path);
      std::ostringstream oss;
      if(!fin.is_open()) {msg::print_error("Couldn't open the file: " + path);return;}

      typedef typename levelset<GridTraitsType, LevelSetTraitsType>::size_type size_type;
      typedef typename levelset<GridTraitsType, LevelSetTraitsType>::index_type index_type;
      typedef typename levelset<GridTraitsType, LevelSetTraitsType>::value_type value_type;

      char buff[10] = {};
      unsigned char byte;
      uint32_t uInt;
      fin.read(buff, 10);
      //Comparing Identification Bytes
      if(std::string(buff).compare(0, 4, "LvSt")) {msg::print_error("File is not a levelset file."); return;}

      const int dim = buff[4]-48;
      if(LVST_FILE_VERSION_NUMBER !=  buff[5]-48) msg::print_warning("File version does not match!");
      if(bigEndian() != buff[6]-48) msg::print_warning("File was written in a different byte order than it is being read. Results may be incorrect!");
      int bits_per_distance = buff[7];
      int char_bit = buff[8]-48;

      if(char_bit != CHAR_BIT) msg::print_error("Bits per byte do not match this operating system.");
#ifdef VERBOSE
      oss << std::endl << "Dimensions: " << dim << std::endl
          << "Bits per distance: " << bits_per_distance << std::endl
          << "Bits per char: " << char_bit << std::endl
          << "Bytes per grid min/max: " << (buff[9] >> 4) << std::endl
          << "Bytes per grid delta: " << (buff[9] & 0xF) << std::endl
          << "Offset file header: " << fin.tellg();
      msg::print_message_2(oss.str());
#endif
      /*
      Skip grid properties because they were read before the levelset:
      From our current position  - fin.tellg()
      We go 2 * Bytes per grid min/max + 1 Byte per boundary condition, per dimension  (  (buff[10]>>4) *2  +1)  *dim
      plus the additional number of bytes for the grid delta  (buff[10] & 0xF)
      */
      fin.seekg(int(fin.tellg()) + ((buff[9]>>4)*2+1)*dim + (buff[9] & 0xF));

      //initialize the empty levelset to sublevelsets size 1
      ls.initialize();

      for(int i=dim;i--;){
        //get the start indices, runtypes and runbreaks vectors
        std::vector<size_type>& startIndices = ls.startIndices(i);
        std::vector<size_type>& runTypes = ls.runTypes(i);
        std::vector<index_type>& runBreaks = ls.runBreaks(i);
        uint32_t num_st_indices, num_run_types, num_run_breaks;
        uint8_t bytesPerStIndex, bytesPerRnType, bytesPerRnBreak;
        //reading in the HRLE header
        fin.read((char *)&bytesPerStIndex, 1);
        fin.read((char *)&bytesPerRnType, 1);
        fin.read((char *)&bytesPerRnBreak, 1);
        fin.read((char *)&num_st_indices, 4);
        fin.read((char *)&num_run_types, 4);
        fin.read((char *)&num_run_breaks, 4);

#ifdef VERBOSE
        oss.str("");
        oss << int(bytesPerStIndex) << " byte(s) per start index." << std::endl
            << int(bytesPerRnType) << " byte(s) per runtype." << std::endl
            << int(bytesPerRnBreak) << " byte(s) per runbreak.";
        msg::print_message_2(oss.str());
#endif
        uint32_t values_read = 0;
        //reading start indices
        if(startIndices.size() > 0) startIndices.clear();
        for(unsigned int x = 0; x < num_st_indices; x++){
          uInt = 0;
          fin.read((char*) &uInt, bytesPerStIndex);
          startIndices.push_back(uInt);
          values_read++;
        }
#ifdef VERBOSE
        oss.str("");
        oss << "Dimension " << i;
        msg::print_message_2(oss.str());
        oss.str("");
        oss << "\t" << values_read << " of " << num_st_indices << " start indices read.";
        msg::print_message_2(oss.str());
#endif
        uint32_t count = 0;
        values_read = 0;
        //reading runtypes
        if(runTypes.size() > 0) runTypes.clear();
        uint32_t bytesToRead = std::ceil(num_run_types/4.0);
        uint8_t rnTypesPerByte = CHAR_BIT/2;
        for(unsigned int y = 0; y < bytesToRead; y++){
          fin.read((char *)&byte, 1);
          for(unsigned int z = rnTypesPerByte; z--;){
            if(values_read == num_run_types) break;
            uInt = byte >> z*2 & 0x3;
            if(uInt == 1) runTypes.push_back(ls.POS_PT);
            else if(uInt == 3) runTypes.push_back(ls.NEG_PT);
            else if(uInt == 2) runTypes.push_back(ls.UNDEF_PT);
            else if(uInt == 0) runTypes.push_back(0), count++;
            values_read++;
          }
        }

#ifdef VERBOSE
        oss.str("");
        oss << "\t" << values_read << " of " << num_run_types << " runtypes read." << " Defined runtypes: " << count;
        msg::print_message_2(oss.str());
#endif
        //reading defined runtypes
        for(unsigned int j=0; j<runTypes.size(); j++){
          if(runTypes[j] == 0){
            uInt = 0;
            fin.read((char *) &uInt, bytesPerRnType);
            runTypes[j] = uInt;
          }
        }

        values_read = 0;

        //Example: We have 2 bytes and the runbreak is -7, which corresponds to FFF9 (Two's complement).
        //Now if we read it into a 4 byte integer we will have 0000 FFF9, which is not interpreted as -7, but as a positive number.
        //-7 as a 4 byte integer is: FFFF FFF9.

        //reading runbreaks
        long long sInt;
        uint64_t sign_bit_mask = 0;
        --sign_bit_mask <<= bytesPerRnBreak * CHAR_BIT;

        if(runBreaks.size() > 0) runBreaks.clear();
        for(unsigned int z = 0; z < num_run_breaks; z++){
          sInt = 0;
          //fin.read((char *) sInt, bytesPerRnBreak);
          fin.read((char *) &sInt, bytesPerRnBreak);
          //if the sign bit is set, fill up the upper bits
          if(sInt >> (bytesPerRnBreak * CHAR_BIT-1) & 0x1 ) sInt |= sign_bit_mask;
          runBreaks.push_back(sInt);
          values_read++;
        }
        #ifdef VERBOSE
              oss.str("");
              oss << "\t" << values_read << " of " << num_run_breaks << " runbreaks." << std::endl;
              msg::print_message_2(oss.str());
        #endif
      }

      uint32_t num_distances = 0, values_read = 0;
      //reading the number of distances to read
      fin.read((char *)&num_distances, 4);
#ifdef VERBOSE
      oss.str("");
      oss << "Offset distances: " << fin.tellg();
      msg::print_message_2(oss.str());
#endif
      std::vector<value_type> & distances = ls.distances();
      long double value = (std::pow((long double) 2, (long double) bits_per_distance)-1.0L)/2.0L;
      int count = char_bit/bits_per_distance;
      if(count < 1) count  = 1;
      const int overflow_num_bytes = bits_per_distance/char_bit -1; //-1 because bpd is mapped to the next higher power of 2
      int bits_per_byte = overflow_num_bytes > 0 ? char_bit : bits_per_distance; //for 8 bits 1 ... 8
      int num_reads = std::ceil((double)num_distances/count);
      unsigned char mask = 0xFF >> bits_per_distance%char_bit; //for bits per distance > 8 the mask will be 0xFF
      unsigned long long m = 0;
      --m >>= sizeof(long long)*CHAR_BIT-bits_per_distance;
      byte = 0;
      //reading distances
      if(distances.size() > 0) distances.clear();
      unsigned long long discrete_distance = 0;
      unsigned long long tmp_distance = 0;
      for(int i = 0; i < num_reads; i++){
        //byte = 0;
        fin.read((char *)&byte, 1);
        for(int z = count; z--;){
          if(values_read == num_distances) break; //if distances are odd, skip padding bits
          discrete_distance = 0;
          //tmp_distance = 0;
          discrete_distance |= byte >> z*bits_per_byte & mask;
          //read in overflow
          for(int j=0; j<overflow_num_bytes; j++){
            //byte = 0;
            fin.read((char *)&byte, 1);
            tmp_distance = byte & mask;
            tmp_distance <<= (j+1) * bits_per_byte;//shift the byte to the correct position
            discrete_distance |= tmp_distance;
            //if you dont use a temporary instead this will fail with bits > 32, because temporary values have 32 bits only!!!!
          }
          discrete_distance &= m; //for some reason the program adds FF as padding. remove those ....
          distances.push_back( (discrete_distance - value) / value); //- 1.0L);

          values_read++;
        }
        std::cout << std::dec;
      }
#ifdef VERBOSE
      oss.str("");
      oss << values_read << " of " << num_distances << " distances read.";
      msg::print_message_2(oss.str());
#endif
      if(!fin.good()) {
        msg::print_message_2("Error while reading from " + path + ".");
        if(fin.eof()) {
          oss.str("");
          long tmp = fin.tellg();
          fin.close();
          oss << "EOF reached on input operation.(" << tmp << ")";
          msg::print_error(oss.str());
        }
        else if(fin.bad()) {
          fin.close();
          msg::print_error("Read error on input operation.");
        }
        else if(fin.fail()) {
          fin.close();
          msg::print_error("Logical error on input operation.");
        }
        else {
          fin.close();
          msg::print_error("Undetermined error on input operation.");
        }
      }
      fin.close();

      ls.finalize(2);
#ifdef VERBOSE
      msg::print_start("NOTE: If the program crashes during prune()/segment(), it is likely that the grid was not read properly or the read/write funtion changed.\nExecuting prune()/segment()...");
#endif
      //print to debug
      ls.print();
      //ls.prune();
      //prune() is called before writing the levelset file, therefore we only need to set up segmentation.
      ls.segment();
#ifdef VERBOSE
      msg::print_done();
#endif
    }

} //NAMESPACE LVLSET END

#endif /*OUTPUT_HPP_*/
