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
#include "kernel.hpp"
#include "levelset2surface.hpp"
#include "../message.h"

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
      char buff[11] = {};
      fin.read(buff, 11);
      const int dim = buff[5]-48;
      const int D = GridTraitsType::dimensions;
      //Comparing Identification Bytes
      if(std::string(buff).compare(0, 5, "LVSTx")) msg::print_error("File is not a levelset file.");
      if(LVST_FILE_VERSION_NUMBER !=  buff[6]-48) msg::print_warning("File version does not match!");
      if(bigEndian() != buff[7]-48) msg::print_warning("File was written in a different byte order than it is being read. Results may be incorrect!");

      double gridDelta;
      int32_t grid_min[D];
      int32_t grid_max[D];
      char bconditions[D];

      //read in the grid properties
      for(int i=dim;i--;){
        fin.read((char *)&grid_min[i], buff[10]>>4);
        fin.read((char *)&grid_max[i], buff[10]>>4);
        fin.read((char *)&bconditions[i], 1);
      }
      fin.read((char *)&gridDelta, buff[10] & 0xF);
      if(fin.fail()) msg::print_error("Couldn't read grid properties from file: " + path);
      fin.close();
      return GridTraitsType(grid_min, grid_max, bconditions, gridDelta);
    }

    template <class GridTraitsType, class LevelSetTraitsType>
    void exportLevelsetToFile(levelset<GridTraitsType, LevelSetTraitsType>& ls, const std::string& path) {
        /*************************************************************************************************************
        ***************************************    THE LEVELSET FILE FORMAT    ***************************************
        **************************************************************************************************************
        *    File Header: 33 - 46 Bytes    *                                                                         *
        ************************************                                                                         *
        *     5 Bytes   Identification Bytes (LVSTx)                                                                 *
        *     1 Byte    Dimension of the Levelset (2 or 3)                                                           *
        *     1 Byte    File Version Number                                                                          *
        *     1 Byte    Endianess - Little Endian (0) or Big Endian (1)                                              *
        *     1 Byte    Bits Per Runtype                                                                             *
        *     1 Byte    Bits Per Distance                                                                            *
        *     1 Byte    This byte contains the number of bytes used for the grid min and grid max of each dimension. *
        *               In addition it also contains the number of bytes used for the grid delta.                    *
        *  For example: 0x48 means that each grid min/max takes 4 bytes each and the grid delta 8 bytes.             *
        *   NOTE: The following block is repeated for each dimension.                                                *
        *     4 Bytes   Grid minimum                                                                                 *
        *     4 Bytes   Grid maximum                                                                                 *
        *     1 Byte    Boundary Condition                                                                           *
        *   NOTE: END                                                                                                *
        *   4-8 Bytes   GridDelta (sizeof double)                                                                    *
        **************************************************************************************************************
        *    H - RLE Block Header: 13 Bytes    *                                                                     *
        ****************************************                                                                     *
        *    1  Byte   This byte contains the number of bytes used for each start index, runtype and runbreak.       *
        *              The bytesize goes from [1, 4], but since only 2 bit are used it is shifted(-1) to [0, 3].     *
        *              //TODO: If bytesize is larger than 4 we need more than one byte to store the information.     *
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
        if(!fout.is_open()) {msg::print_error("Couldn't open the file: " + path); return;}

        //type & constant redefinitions
        typedef typename levelset<GridTraitsType, LevelSetTraitsType>::size_type size_type;
        typedef typename levelset<GridTraitsType, LevelSetTraitsType>::index_type index_type;
        typedef typename levelset<GridTraitsType, LevelSetTraitsType>::value_type value_type;
        const unsigned int D = ls.dimensions;

        //write the file header
        fout << "LVSTx" <<  D << LVST_FILE_VERSION_NUMBER << (bigEndian() ? 1 : 0)  << BITS_PER_RUNTYPE << BITS_PER_DISTANCE;
        const index_type* gridMinima = ls.grid().grid_traits().getMinima();
        const index_type* gridMaxima = ls.grid().grid_traits().getMaxima();
        const boundary_type* gridBConditions = ls.grid().grid_traits().getBoundaryConditions();
        char byte = 0;
        //write grid properties
        byte |= BYTES_GRID_LIMITS << 4;//bytes used to store grid min and grid max; boundary condition is always 1 byte
        byte |= sizeof(double);//bytes used to store grid delta
        fout << byte;
        for (int dim=D;dim--;) {
            fout.write((char *)&gridMinima[dim], BYTES_GRID_LIMITS);
            fout.write((char *)&gridMaxima[dim], BYTES_GRID_LIMITS);
            fout.write((char *)&gridBConditions[dim], 1);
        }
        double delta = ls.grid().grid_traits().grid_position(0, 1);//1st argument is not used. Returns 1 * GridDelta.
        fout.write((char *)&delta, sizeof(double));

        int count;

        // makes parallelized sublevelset into one levelset
        ls.serialize();

        for (int dim=D-1;dim>=0;--dim) {
          //get start indices, runbreaks and runtypes
          //NOTE: the functions return non-const references
          const std::vector<size_type> & startIndices = ls.startIndices(dim);
          const std::vector<size_type> & runTypes = ls.runTypes(dim);
          const std::vector<index_type> & runBreaks = ls.runBreaks(dim);

          //Getting the largest possible values for a start index, runtype and runbreak
          //This is done to use the minimum amount of bytes per start index, runtype and runbreak per dimension.

          //The largest start index is <= the size of runtypes. (start index points to a runtype)

          //The largest runtype is <= the size of start indices of the next dimension.
          //If dim == 0 the largest runtype <= the size of distances. (defined runtype points to a start index/distance)

          //The smallest runbreak is >= grid_min.
          //The largest runbreak is <= grid_max.

    //FORCE_BYTESIZE lets you manually override the bytesize if you ever need to (you shouldn't)
    #undef FORCE_BYTESIZE
    #ifdef FORCE_BYTESIZE
          //if bytesize has to be forced manually
          uint32_t bytesPerStIndex = BYTES_START_INDEX;
          uint32_t bytesPerRnType = BYTES_RUNTYPE;
          int32_t bytesPerRnBreak = BYTES_RUNBREAK;
    #else
          uint32_t bytesPerStIndex;// = *(startIndices.end()-1);
          uint32_t bytesPerRnType = (dim == 0) ? ls.distances().size() : ls.startIndices(dim-1).size();
          int32_t bytesPerRnBreak;

          //bytesize for start indices
          if(runTypes.size() <= UINT8_MAX) bytesPerStIndex = 1;
          else if(runTypes.size() <= UINT16_MAX) bytesPerStIndex = 2;
              else if(runTypes.size() <= UINT24_MAX) bytesPerStIndex = 3;
          else bytesPerStIndex = 4;
          //bytesize for runtypes
          if(bytesPerRnType <= UINT8_MAX) bytesPerRnType = 1;
          else if(bytesPerRnType <= UINT16_MAX) bytesPerRnType = 2;
              else if(bytesPerRnType <= UINT24_MAX) bytesPerRnType = 3;
          else bytesPerRnType = 4;
          //bytesize for runbreaks
          if(gridMinima[dim] >= INT8_MIN && gridMaxima[dim] <= INT8_MAX) bytesPerRnBreak = 1;
          else if(gridMinima[dim] >= INT16_MIN && gridMaxima[dim] <= INT16_MAX) bytesPerRnBreak = 2;
          else if(gridMinima[dim] >= INT24_MIN && gridMaxima[dim] <= INT24_MAX) bytesPerRnBreak = 3;
          else bytesPerRnBreak = 4;
    #endif
          byte = 0;
          //use 2 bits to store the number of bytes needed
          //The number of bytes can be [1, 4], but with 2 bits you can store only [0, 3]
          byte |= bytesPerStIndex-1;
          byte |= (bytesPerRnType-1) << 2;
          byte |= (bytesPerRnBreak-1) << 4;

          //write 13 byte H-RLE block header
          fout.write(&byte, 1);
          uint32_t num = startIndices.size();
          fout.write((char *)&num, 4);
          num = runTypes.size();
          fout.write((char *)&num, 4);
          num = runBreaks.size();
          fout.write((char *)&num, 4);

          //write the start indices to the file
          for (typename std::vector<size_type>::const_iterator it=startIndices.begin();it!=startIndices.end();++it) {
              fout.write((char *)&(*it), bytesPerStIndex);
          }

          //write all runtypes to the file, skipping all segments and indices (using 2 bits per runtype)
          count = 3;byte=0;
          std::vector<size_type> def_run_indices = {}; //store all indices for defined runtypes
          for (typename std::vector<size_type>::const_iterator it=runTypes.begin();it!=runTypes.end();++it) {
            if(*it == ls.POS_PT){// 01 - positive undefined runtype
              byte |= 1 << count * BITS_PER_RUNTYPE;
              count--;
            }
            else if(*it == ls.NEG_PT){// 11 - negative undefined runtype
              byte |= 3 << count * BITS_PER_RUNTYPE;
              count--;
            }
            else if(*it == ls.UNDEF_PT){// 10 - uninitialized runtype
              byte |= 2 << count * BITS_PER_RUNTYPE;
              count--;
            }
            else if(!ls.is_defined(*it)){// skip Segments --> there shouldn't be any because of serialize()
              msg::print_warning("Segment detected during levelset file export.");
              continue;
            }
            else {// 00 - defined runtype
                byte |= 0 << count * BITS_PER_RUNTYPE;
                def_run_indices.push_back(*it);
                count--;
            }
            if(count < 0){ //if 4 runtypes are written into the byte, write it to the file
              fout << byte;
              count = 3;byte=0;
            }
          }
          if(count >= 0 && count < 3) //number of runtypes % 4 > 0
              fout << byte;

          //write defined runtypes
          for (typename std::vector<size_type>::const_iterator it=def_run_indices.begin();it!=def_run_indices.end();++it) {
            fout.write((char *)&(*it), bytesPerRnType);
          }

          //Write runbreaks
          for (typename std::vector<index_type>::const_iterator it=runBreaks.begin();it!=runBreaks.end();++it) {
            fout.write((char *)&(*it), bytesPerRnBreak);
          }
        }
        //get distances
        const std::vector<value_type> & distances = ls.distances();
        long num = distances.size();
        fout.write((char *)&num, 4);
        count = CHAR_BIT/BITS_PER_DISTANCE-1;
        byte = 0;
        double value = (std::pow(2, BITS_PER_DISTANCE)-1)/2;

        for (typename std::vector<value_type>::const_iterator it=distances.begin();it!=distances.end();++it) {
          //  Levelset values range from -1 .... +1
          //  With n bits we can represent values from 0 .... (2^n -1)
          //         -1 .... +1          |*(2^n-1)/2
          // -(2^n-1)/2 .... +(2^n-1)/2  |+(2^n-1)/2
          //          0 .... +(2^n-1)
          byte |= std::lround(*it * value + value) << count * BITS_PER_DISTANCE;
          count--;
          if(count < 0){
            fout << byte;
            count = CHAR_BIT/BITS_PER_DISTANCE-1;byte=0;
          }
        }
        if(count == 0)
            fout << byte;

        if(fout.fail()) msg::print_error("ERROR: Couldn't write to file: " + path);
        fout.close();
        ls.finalize(2);
        ls.thin_out();// sets up segmentation for one levelset, calls finalize
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

      //std::cout << "Reading in LevelSet from " << path << std::endl;
      char buff[11] = {};
      char byte;
      uint32_t uInt;
      fin.read(buff, 11);
      //Comparing Identification Bytes
      if(std::string(buff).compare(0, 5, "LVSTx")) {msg::print_error("File is not a levelset file."); return;}

      const int dim = buff[5]-48;
      if(LVST_FILE_VERSION_NUMBER !=  buff[6]-48) msg::print_warning("File version does not match!");
      if(bigEndian() != buff[7]-48) msg::print_warning("File was written in a different byte order than it is being read. Results may be incorrect!");
      int bits_per_runtype = buff[8]-48;
      int bits_per_distance = buff[9]-48;
#ifdef VERBOSE
      oss << std::endl << "Dimensions: " << dim << std::endl << "Bits per runtype:" << bits_per_runtype << std::endl
          << "Bits per distance:" << bits_per_distance << std::endl << "Bytes per grid min/max: " << (buff[10] >> 4) << std::endl
          << "Bytes per grid delta: " << (buff[10] & 0xF);
      msg::print_message_2(oss.str());
#endif
      /*
      Skip grid properties because they was read before we read the levelset:
      From our current position  fin.tellg()
      We go 2 * Bytes per grid min/max + 1 Byte per boundary condition, per dimension  (  (buff[10]>>4) *2  +1)  *dim
      plus the additional number of bytes for the grid delta  (buff[10] & 0xF)
      */
      fin.seekg(int(fin.tellg()) + ((buff[10]>>4)*2+1)*dim + (buff[10] & 0xF));

      //initialize the empty levelset to sublevelsets size 1
      ls.initialize();

      for(int i=dim;i--;){
        //get the start indices, runtypes and runbreaks vectors
        std::vector<size_type>& startIndices = ls.startIndices(i);
        std::vector<size_type>& runTypes = ls.runTypes(i);
        std::vector<index_type>& runBreaks = ls.runBreaks(i);
        uint32_t num_st_indices, num_run_types, num_run_breaks;
        int32_t bytesPerStIndex, bytesPerRnType, bytesPerRnBreak;
        //reading in the HRLE header
        fin.read(&byte, 1);
        fin.read((char *)&num_st_indices, 4);
        fin.read((char *)&num_run_types, 4);
        fin.read((char *)&num_run_breaks, 4);

        bytesPerStIndex = (byte & 0x3) +1;
        bytesPerRnType = (byte >> 2 & 0x3) +1;
        bytesPerRnBreak = (byte >> 4 & 0x3) +1;
#ifdef VERBOSE
        oss.str("");
        oss << bytesPerStIndex << " byte(s) per start index." << std::endl
            << bytesPerRnType << " byte(s) per runtype." << std::endl
            << bytesPerRnBreak << " byte(s) per runbreak.";
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
        for(int y = 0; y < std::ceil(num_run_types/4.0); y++){
          fin.read(&byte, 1);
          for(int z = 4; z--;){
            if(values_read == num_run_types) break;
            uInt = byte >> z*bits_per_runtype & 0x3;
            if(uInt == 1) runTypes.push_back(ls.POS_PT);
            else if(uInt == 3) runTypes.push_back(ls.NEG_PT);
            else if(uInt == 2) runTypes.push_back(ls.UNDEF_PT);
            else if(uInt == 0) runTypes.push_back(101), count++;
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
          if(runTypes[j] == 101){
            uInt = 0;
            fin.read((char *) &uInt, bytesPerRnType);
            runTypes[j] = uInt;
          }
        }

        values_read = 0;

        //Allocate the exact amount of bytes you need for a runbreak
        //Example: We have 2 bytes and the runbreak is -7, which corresponds to FFF9 (Two's complement).
        //Now if we read it into a 4 byte integer we will have 0000 FFF9, which is not interpreted as -7, but as a positive number.
        //-7 as a 4 byte integer is: FFFF FFF9.

        //int8_t *sInt = new int8_t[bytesPerRnBreak];
        int8_t *sInt = (int8_t *) malloc(bytesPerRnBreak*sizeof(int8_t));
        //reading runbreaks
        if(runBreaks.size() > 0) runBreaks.clear();
        for(unsigned int z = 0; z < num_run_breaks; z++){
          fin.read((char *) sInt, bytesPerRnBreak);
          runBreaks.push_back(*sInt);
          values_read++;
        }
        //delete[] sInt;
        free(sInt);
#ifdef VERBOSE
        oss.str("");
        oss << "\t" << values_read << " of " << num_run_breaks << " runbreaks.";
        msg::print_message_2(oss.str());
#endif
      }

      uint32_t num_distances = 0, values_read = 0;
      //reading the number of distances to read
      fin.read((char *)&num_distances, 4);
      std::vector<value_type> & distances = ls.distances();

      double value = (std::pow(2, bits_per_distance)-1)/2;
      int count = CHAR_BIT/bits_per_distance;
      int num_reads = std::ceil(num_distances/count);
      char mask = 0xFF >> (CHAR_BIT-bits_per_distance);
      //reading distances
      if(distances.size() > 0) distances.clear();
      uint8_t uInt8;
      for(int i = 0; i < num_reads; i++){
        fin.read(&byte, 1);
        for(unsigned int z = count; z--;){
          if(values_read == num_distances) break; //if distances are odd, skip padding bits
          uInt8 = byte >> z*bits_per_distance & mask;
          distances.push_back(uInt8 / value - 1.0);
          values_read++;
        }
      }
#ifdef VERBOSE
      oss.str("");
      oss << values_read << " of " << num_distances << " distances read.";
      msg::print_message_2(oss.str());
#endif
      if(fin.fail()) msg::print_error("Couldn't read from file: " + path);
      fin.close();

      ls.finalize(2);
      ls.thin_out();
    }

} //NAMESPACE LVLSET END

#endif /*OUTPUT_HPP_*/
