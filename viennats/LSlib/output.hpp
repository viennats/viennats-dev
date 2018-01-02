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

#include "kernel.hpp"
#include "levelset2surface.hpp"
#include <float.h>

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
    void write_explicit_hollow_surface_vtk(const levelset<GridTraitsType, LevelSetTraitsType>& l, const std::string& filename, typename LevelSetTraitsType::value_type eps=0.){
        write_explicit_hollow_surface_vtk(l, filename, DefaultDataType(), eps);
    }

    template <class GridTraitsType, class LevelSetTraitsType, class DataType>
    void write_explicit_hollow_surface_vtk(const levelset<GridTraitsType, LevelSetTraitsType>& l, const std::string& filename, const DataType& Data, typename LevelSetTraitsType::value_type eps=0.){

        const int D=levelset<GridTraitsType, LevelSetTraitsType>::dimensions;

        typename GetActivePointType<typename LevelSetTraitsType::size_type, DataType>::result ActivePointList;
        Surface<D> s;

        extract(l, s, eps, ActivePointList, true);

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
}

#endif /*OUTPUT_HPP_*/
