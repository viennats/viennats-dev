#ifndef GEOMETRY_H_
#define GEOMETRY_H_

/* =========================================================================
   Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.

                            -----------------
                 ViennaTS - The Vienna Topography Simulator
                            -----------------

   Contact:         viennats@iue.tuwien.ac.at

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


#include <vector>
#include <string>
#include <fstream>
#include <bitset>
#include <set>
#include "vector.hpp"
#include "message.h"
#include <stdexcept>

namespace geometry {

	template <int D> class geometry {
	public:

		std::vector<lvlset::vec<double, D> > Nodes;
		std::vector<lvlset::vec<unsigned int, D+1> > Elements;
		std::vector<int> Materials;

		lvlset::vec<double, D> Min,Max; //for bounding box

		geometry()  {}

    // Silvaco Structure Reader
                void ReadSTR(std::string   const& FileName,
              double               scale,//=1.,
              std::vector<int>   & InputTransformationDirections,//=std::vector<int>(),
              std::vector<bool>  & InputTransformationSigns,//=std::vector<bool>(),
              bool                change_input_parity,//=false,
              std::vector<double>& shift,//=std::vector<double>()
              std::vector<int>   & ignore_materials
             )
    {
      std::cout << "Reading STR file!";


      std::ifstream reader(FileName.c_str(), ios::binary);

      if (!reader)
        throw std::runtime_error("Cannot open file: "+FileName+" - Aborting!");

      std::string line;
      std::string dump;

      int dim_geometry;
      int dim_topology;
      int number_of_points;
      int number_of_triangles;
      int number_of_tetrahedrons;


      while(1)
      {
        std::getline(reader, line);

        if(!line.empty() && line.at(0) == 'k')
        {
          // traverse the line via stringstream
          std::stringstream ss(line);
          ss >> dump;
          ss >> dim_geometry;
          ss >> number_of_points;
          ss >> dim_topology;
          ss >> number_of_triangles;
          ss >> number_of_tetrahedrons;
          break;
        }
      }

      std::cout << "dim geometry: " << dim_geometry << std::endl;
      std::cout << "dim topology: " << dim_topology << std::endl;
      std::cout << "number of points: " << number_of_points << std::endl;
      std::cout << "number of triangles: " << number_of_triangles << std::endl;
      std::cout << "number of tetrahedrons: " << number_of_tetrahedrons << std::endl;

      if(D != dim_geometry)
        throw std::runtime_error("Geometry dimension of reader does not fit input file - Aborting!");

      // read one more dummy line, won't need this one, this is just to advance the file pointer
      std::getline(reader, line);

      // -------------------------------
      // import points
      //
      Nodes.resize(number_of_points);
      for (int i=0;i<number_of_points;i++)
      {
        double coords[D];

        // c - ignored
        reader >> dump;
        // point-id - ignored
        reader >> dump;

        // x y z
        reader >> coords[0];
        reader >> coords[1];
        if(D == 3) reader >> coords[2];

        for (int j=0;j<D;j++)
        {
          Nodes[i][j]=coords[InputTransformationDirections[j]];
          int shift_size=shift.size();
          if (shift_size>j) Nodes[i][j]+=shift[j];
          if (InputTransformationSigns[j]) Nodes[i][j]=-Nodes[i][j];
          Nodes[i][j]*=scale;
        }
      }

      // -------------------------------
      // import elements (+materials), i.e., triangles (2D) or tetrahedrons (3D)
      //
      int region_id;
//      int vertex_id;
//      int cell_indices[D+1];
      lvlset::vec<unsigned int, D+1>  cell_indices;
      bool ignore;
      if(D == 3)
      {
        // in the 3D case, we need to skip the triangles which are also
        // provided in the input file. Therefore we have to advance the file pointer
        // by the number of triangles - the '+1' is required to move the file poiniter
        // into the first 'tetrahedron' line
        for (int i=0;i<number_of_triangles+1;i++)
          std::getline(reader, line);

//        Elements.resize(number_of_tetrahedrons);
//        Materials.resize(number_of_tetrahedrons);

        for (int i=0;i<number_of_tetrahedrons;i++)
        {
          // T - ignored
          reader >> dump;
//          std::cout << "Type: " << dump << std::endl;
          // tetrahedron-id - ignored
          reader >> dump;
//          std::cout << "ID: " << dump << std::endl;
          // region
//          reader >> Materials[i];
          reader >> region_id;
//          std::cout << "Region: " << Materials[i] << std::endl;

          if(std::find(ignore_materials.begin(), ignore_materials.end(), region_id) != ignore_materials.end())
          {
            ignore = true;
          }
          else ignore = false;

          if(!ignore)
            Materials.push_back(region_id);

          // tetrahedron vertex ids
          for (int j=0;j<D+1;j++)
          {
            reader >> cell_indices[j];

//            reader >> Elements[i][j];
            // Silvaco uses 1-based indices, we use 0-based,
            // so we need to manually decrease it by one to map the index spaces correctly
//            Elements[i][j]--;
            cell_indices[j]--;
          }

          if(!ignore)
          {
            if (change_input_parity) std::swap(cell_indices[0], cell_indices[1]);
            Elements.push_back(cell_indices);
          }

          // advance file pointer, omitting the additional meta data for each tetrahedron
          for (int j=0;j<D+1;j++) reader >> dump;
        }
      }
      else
        throw std::runtime_error("2D files are currently not supported - Aborting!");

      // correct material IDs
      // note: material IDs need to start with 1 and have a 1-increment ID
      // if we ignored some materials, having a small ID, e.g., 1, we need to decrement all other
      // material ids accordingly
      // we start from the largest region to be ignored and process the material set one after another
      for(std::vector<int>::reverse_iterator rit = ignore_materials.rbegin(); rit != ignore_materials.rend(); ++rit)
      {
        for(std::vector<int>::iterator it = Materials.begin(); it != Materials.end(); ++it)
        {
          if(*it > *rit) --(*it);
        }
      }

      reader.close();
    }

    // DX Reader
                void ReadDX(std::string   const& FileName,
                double               scale,//=1.,
                std::vector<int>   & InputTransformationDirections,//=std::vector<int>(),
                std::vector<bool>  & InputTransformationSigns,//=std::vector<bool>(),
                bool                change_input_parity,//=false,
                std::vector<double>& shift//=std::vector<double>()
               )
    {
      std::cout << "Reading DX file!";
                        std::ifstream f(FileName.c_str());

                        if (!f) msg::print_error("Failed reading geometry file!");

                        std::string c;

			//read nodes
			std::getline(f,c);

			int num_nodes=atoi(&c[63]);

			Nodes.resize(num_nodes);
			for (int i=0;i<num_nodes;i++) {
			    double coords[D];

				for (int j=0;j<D;j++) f>>coords[j];

				for (int j=0;j<D;j++) {
				    Nodes[i][j]=coords[InputTransformationDirections[j]];

				    int shift_size=shift.size();
				    if (shift_size>j) Nodes[i][j]+=shift[j];
		    if (InputTransformationSigns[j]) Nodes[i][j]=-Nodes[i][j];
		    Nodes[i][j]*=scale;
				}
			}

			std::getline(f,c);

			//read elements
			std::getline(f,c);
			int num_elems=atoi(&c[63]);

			Elements.resize(num_elems);
			for (int i=0;i<num_elems;i++) {
				for (int j=0;j<D+1;j++) f>>Elements[i][j];
				if (change_input_parity) std::swap(Elements[i][0],Elements[i][1]);
			}


			std::getline(f,c);
			std::getline(f,c);
			std::getline(f,c);
			std::getline(f,c);

			//read materials
			Materials.resize(num_elems);
			for (int i=0;i<num_elems;i++) {
				f>>Materials[i];
			}

                        f.close();
    }

                void Read(std::string  const& FileName,
              double              scale,//=1.,
              std::vector<int>    InputTransformationDirections,//=std::vector<int>(),
              std::vector<bool>   InputTransformationSigns,//=std::vector<bool>(),
              bool                change_input_parity,//=false,
              std::vector<int>    MapMaterials,//=std::vector<int>(),
              std::vector<double> shift,//=std::vector<double>()
              std::vector<int>    ignore_materials
             )
    {


      // make sure that the material IDs to be ignored are sorted
      std::sort(ignore_materials.begin(), ignore_materials.end());

      //input transformation
      while(InputTransformationDirections.size()<D) InputTransformationDirections.push_back(InputTransformationDirections.size());        //TODO test if directions are unique
      while(InputTransformationSigns.size()<D) InputTransformationSigns.push_back(false);

      if ((InputTransformationDirections[0]+1)%D!=InputTransformationDirections[1]) change_input_parity=!change_input_parity;

      for (int i=0;i<D;++i)
        if (InputTransformationSigns[i]) change_input_parity=!change_input_parity;
      // determine whether we need the DX or the STR reader - call appropriate subroutine accordingly
      if(FileName.substr(FileName.find_last_of(".") + 1) == "dx")
        ReadDX(FileName, scale, InputTransformationDirections, InputTransformationSigns, change_input_parity, shift);
      else if(FileName.substr(FileName.find_last_of(".") + 1) == "str")
        ReadSTR(FileName, scale, InputTransformationDirections, InputTransformationSigns, change_input_parity, shift, ignore_materials);
      else throw std::runtime_error("Input geometry file type \""+FileName.substr(FileName.find_last_of(".") + 1)+"\" not supported - Aborting!");


      //map materials - this is the same regardless of the input file format
      if (!MapMaterials.empty())
      {
        std::vector<lvlset::vec<unsigned int, D+1> > oldElements;
        std::vector<int> oldMaterials;

        std::swap(oldElements, Elements);
        std::swap(oldMaterials, Materials);

        for (int a=0;a<oldElements.size();++a)
        {
          bool add=true;
          int mat=oldMaterials[a];
          if (mat<static_cast<int>(MapMaterials.size()))
          {
            mat=MapMaterials[mat];
            if (mat<=0) add=false;
          }
          else
          {
            assert(0);  //Material mapping failed
            mat=false;
          }
          if (add)
          {
            Elements.push_back(oldElements[a]);
            Materials.push_back(mat);
          }
        }
      }

      //determine min and max of geometry (bounding box)
      CalculateExtensions();

                }

		void Write(const std::string& FileName) const {
			std::ofstream f(FileName.c_str());

			f << "# vtk DataFile Version 2.0" << std::endl;
			f << D << "D Volume" << std::endl;
			f << "ASCII" << std::endl;
			f << "DATASET UNSTRUCTURED_GRID" << std::endl;
			f << "POINTS " << Nodes.size() << " float" << std::endl;

			//!print positions
			for (unsigned int i=0;i<Nodes.size();i++) {
				for (int j=0;j<D;j++) f << static_cast<float>(Nodes[i][j]) << " ";
				if (D==2) f << "0. ";
				f<< std::endl;
			}
			f << "CELLS " << Elements.size() << " " << ((D+2)*Elements.size()) << std::endl;
			for (unsigned int i=0;i<Elements.size();i++) {
				f << (D+1) << " ";
				for (int j=0;j<D+1;j++) f<< Elements[i][j] << " ";
				f << std::endl;
			}

			f << "CELL_TYPES " << Elements.size() << std::endl;
			for (unsigned int i=0;i<Elements.size();i++) {
				f<< ((D==3)?"10":"5") << std::endl;
			}
			f << "CELL_DATA " << Elements.size() << std::endl;
			f << "SCALARS Material int 1" << std::endl;
			f << "LOOKUP_TABLE default" << std::endl;
			for (unsigned int i=0;i<Elements.size();i++) {
				f<< Materials[i] << std::endl;
			}
			f.close();


		}

		void CalculateExtensions() {
			assert(!Elements.empty());
			Min=Nodes[Elements[0][0]];
			Max=Nodes[Elements[0][0]];
			/*for (unsigned int i=1;i<Nodes.size();i++) {
				Min=lvlset::Min(Min,Nodes[i]);
				Max=lvlset::Max(Max,Nodes[i]);
			}*/
			for (unsigned int a=0;a<Elements.size();++a) {
				for (unsigned int i=0;i<D+1;i++) {
					Min=lvlset::Min(Min,Nodes[Elements[a][i]]);
					Max=lvlset::Max(Max,Nodes[Elements[a][i]]);
				}
			}
		}
	};



	template <int D> class surface {


	public:

		typedef unsigned int node_index_type;
		typedef unsigned int element_index_type;

		std::vector<lvlset::vec<double, D> > Nodes;
		std::vector<lvlset::vec<unsigned int, D> > Elements;

		unsigned int number_of_nodes() const {
			return Nodes.size();
		}

		unsigned int number_of_elements() const {
			return Elements.size();
		}

		unsigned int element_node_id(unsigned ElementIndex, int NodeNumber) const {
			return Elements[ElementIndex][NodeNumber];
		}

		double node_coordinate(unsigned int NodeIndex, int Dimension) const {
			return Nodes[NodeIndex][Dimension];
		}

		void Write(std::string FileName) const {
			std::ofstream f(FileName.c_str());

			//!print positions
			f<< "object \"positions\" class array type float rank 1 shape " << D << " items "<< Nodes.size() <<" data follows" << std::endl;
			for (unsigned int i=0;i<Nodes.size();i++) {
				for (int j=0;j<D;j++) f << static_cast<float>(Nodes[i][j]) << " ";
				f<< std::endl;
			}

			//! print connections
			f << "object \"connections\" class array type int rank 1 shape " << D << " items "<< Elements.size() <<" data follows" << std::endl;
			for (unsigned int i=0;i<Elements.size();i++) {
				for (int j=0;j<D;j++) f<< Elements[i][j] << " ";
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

		void WriteVTK(std::string FileName) const {
			std::ofstream f(FileName.c_str());

			f << "# vtk DataFile Version 2.0" << std::endl;
			f << D << "D Surface" << std::endl;
			f << "ASCII" << std::endl;
			f << "DATASET UNSTRUCTURED_GRID" << std::endl;
			f << "POINTS " << Nodes.size() << " float" << std::endl;

			//!print positions
			for (unsigned int i=0;i<Nodes.size();i++) {
				for (int j=0;j<D;j++) f << static_cast<float>(Nodes[i][j]) << " ";
				if (D==2) f << "0. ";
				f<< std::endl;
			}
			f << "CELLS " << Elements.size() << " " << ((D+1)*Elements.size()) << std::endl;
			for (unsigned int i=0;i<Elements.size();i++) {
				f << D << " ";
				for (int j=0;j<D;j++) f<< Elements[i][j] << " ";
				f << std::endl;
			}

			f << "CELL_TYPES " << Elements.size() << std::endl;
			for (unsigned int i=0;i<Elements.size();i++) {
				f<< ((D==3)?"5":"3") << std::endl;
			}

			f.close();
		}

		void ReadVTK(std::string FileName)  {

			std::ifstream f(FileName.c_str());

			std::string c;

			std::getline(f,c);

			std::getline(f,c);
			assert(atoi(&c[0])==D);

			std::getline(f,c);

			std::getline(f,c);

			std::getline(f,c);
			int num_nodes=atoi(&c[7]);

			std::cout << num_nodes << std::endl;

			Nodes.resize(num_nodes);

			for (int i=0;i<num_nodes;i++) {
				for (int j=0;j<D;j++) f>>Nodes[i][j];
			}

			std::getline(f,c);
			std::getline(f,c);
			int num_elems=atoi(&c[6]);

			std::cout << num_elems << std::endl;

			Elements.resize(num_elems);
			int tmp;

			for (unsigned int i=0;i<Elements.size();i++) {
				f >> tmp;
				for (int j=0;j<D;j++) f>> Elements[i][j];
			}
			f.close();
		}


	};


	template <int D, class SurfacesType> void TransformGeometryToSurfaces(
		const geometry<D>& Geometry,
		SurfacesType &Surfaces,
		std::bitset<2*D> remove_flags,
		double eps) {


		//determine maximum number of materials
		unsigned int max_mat= *std::max_element(Geometry.Materials.begin(),Geometry.Materials.end());
		assert(max_mat>=1);


		typedef std::map<lvlset::vec<unsigned int,D>, std::pair<unsigned int, unsigned int> > triangle_map;
		triangle_map Triangles;

		for (unsigned int i=0;i<Geometry.Elements.size();++i) {

			lvlset::vec<unsigned int,D> tmp;

			for (int j=0;j<D+1;j++) {

				for (int k=0;k<D;k++) tmp[k]=Geometry.Elements[i][(j+k)%(D+1)];

				std::bitset<2*D> flags;
				flags.set();

				//if triangle at border skip

				for (int k=0;k<D;k++) {
					for (int l=0;l<D;l++) {
						if (Geometry.Nodes[tmp[k]][l]<Geometry.Max[l]-eps) {
							flags.reset(l+D);
						}
						if (Geometry.Nodes[tmp[k]][l]>Geometry.Min[l]+eps) {
			    flags.reset(l);
			}
					}
				}

				flags &=remove_flags;

				//if (is_open_boundary_negative) flags.reset(open_boundary_direction); else flags.reset(open_boundary_direction+D);

				if (flags.any()) continue;

				tmp.sort();

				lvlset::vec<double,D> pts[D+1];
				for (int k=0;k<D;k++) pts[k]=Geometry.Nodes[tmp[k]];
				pts[D]=Geometry.Nodes[Geometry.Elements[i][(j+D)%(D+1)]];

				typename triangle_map::iterator it=Triangles.lower_bound(tmp);
				if ((it!=Triangles.end()) && (it->first==tmp)) {
					if (lvlset::Orientation(pts)) {
						assert(it->second.second==max_mat+1);
						it->second.second=Geometry.Materials[i];
					} else {
						assert(it->second.first==max_mat+1);
						it->second.first=Geometry.Materials[i];
					}

					if (it->second.first==it->second.second) Triangles.erase(it);

				} else {
					if (lvlset::Orientation(pts)) {
						Triangles.insert(it,std::make_pair(tmp,std::make_pair(max_mat+1,Geometry.Materials[i])));
					} else {
						Triangles.insert(it,std::make_pair(tmp,std::make_pair(Geometry.Materials[i],max_mat+1)));
					}
				}
			}
		}


		Surfaces.resize(max_mat);

		//for all materials/for each surface
		typename SurfacesType::iterator srf_it=Surfaces.begin();
		for (unsigned int m=0;m<max_mat;++m) {

			for (typename triangle_map::iterator it=Triangles.begin();it!=Triangles.end();++it) {
				if ((m>=it->second.first-1) && (m<it->second.second-1)) {
					srf_it->Elements.push_back(it->first);
				} else if ((m>=it->second.second-1) && (m<it->second.first-1)) {
					srf_it->Elements.push_back(it->first);
					std::swap(srf_it->Elements.back()[0],srf_it->Elements.back()[1]);
				}
			}

			//replace Nodes of Geometry by Nodes of individual surface

			const unsigned int undefined_node=std::numeric_limits<unsigned int>::max();
			std::vector<unsigned int> NodeReplacements(Geometry.Nodes.size(),undefined_node);
			unsigned int NodeCounter=0;

			for (unsigned int k=0;k<srf_it->Elements.size();++k) {

				for (int h=0;h<D;h++) {
					unsigned int origin_node=srf_it->Elements[k][h];
					if (NodeReplacements[origin_node]==undefined_node) {
						NodeReplacements[origin_node]=NodeCounter++;
						srf_it->Nodes.push_back(Geometry.Nodes[origin_node]);

					}
					srf_it->Elements[k][h]=NodeReplacements[origin_node];
				}
			}
			++srf_it;
		}
	}

	template <int D, class SurfacesType> void TransformGeometryToSurfaces2(
		const geometry<D>& Geometry,
		SurfacesType &Surfaces,
		std::bitset<2*D> remove_flags,
		double eps) {


		//determine maximum number of materials
		unsigned int max_mat= *std::max_element(Geometry.Materials.begin(),Geometry.Materials.end());
		assert(max_mat>=1);


		typedef std::map<lvlset::vec<unsigned int,D>, std::pair<unsigned int, unsigned int> > triangle_map;
		triangle_map Triangles;

		for (unsigned int i=0;i<Geometry.Elements.size();++i) {

			lvlset::vec<unsigned int,D> tmp;

			for (int j=0;j<D+1;j++) {

				for (int k=0;k<D;k++) tmp[k]=Geometry.Elements[i][(j+k)%(D+1)];

				std::bitset<2*D> flags;
				flags.set();

				//if triangle at border skip

				for (int k=0;k<D;k++) {
					for (int l=0;l<D;l++) {
						if (Geometry.Nodes[tmp[k]][l]<Geometry.Max[l]-eps) {
							flags.reset(l+D);
						}
						if (Geometry.Nodes[tmp[k]][l]>Geometry.Min[l]+eps) {
			    flags.reset(l);
			}
					}
				}

				flags &=remove_flags;

				//if (is_open_boundary_negative) flags.reset(open_boundary_direction); else flags.reset(open_boundary_direction+D);

				if (flags.any()) continue;

				tmp.sort();

				lvlset::vec<double,D> pts[D+1];
				for (int k=0;k<D;k++) pts[k]=Geometry.Nodes[tmp[k]];
				pts[D]=Geometry.Nodes[Geometry.Elements[i][(j+D)%(D+1)]];

				typename triangle_map::iterator it=Triangles.lower_bound(tmp);
				if ((it!=Triangles.end()) && (it->first==tmp)) {
					if (lvlset::Orientation(pts)) {
						assert(it->second.second==max_mat+1);
						it->second.second=Geometry.Materials[i];
					} else {
						assert(it->second.first==max_mat+1);
						it->second.first=Geometry.Materials[i];
					}

					if (it->second.first==it->second.second) Triangles.erase(it);

				} else {
					if (lvlset::Orientation(pts)) {
						Triangles.insert(it,std::make_pair(tmp,std::make_pair(max_mat+1,Geometry.Materials[i])));
					} else {
						Triangles.insert(it,std::make_pair(tmp,std::make_pair(Geometry.Materials[i],max_mat+1)));
					}
				}
			}
		}


		Surfaces.resize(max_mat);

		//for all materials/for each surface
		typename SurfacesType::iterator srf_it=Surfaces.begin();
		for (unsigned int m=0;m<max_mat;++m) {
//			std::cout << "m: " << m << std::endl;
			for (typename triangle_map::iterator it=Triangles.begin();it!=Triangles.end();++it) {
//				std::cout << "second.first-1: " << it->second.first-1 << std::endl;
//				std::cout << "second.second-1: " << it->second.second-1 << std::endl;
				if ((m==it->second.first-1) && (m<it->second.second-1)) {
//					std::cout << "FIRST!\n";
					srf_it->Elements.push_back(it->first);
				} else if ((m==it->second.second-1) && (m<it->second.first-1)) {
//					std::cout << "SECOND!\n";
					srf_it->Elements.push_back(it->first);
					std::swap(srf_it->Elements.back()[0],srf_it->Elements.back()[1]);
				}
/*				if ((m>=it->second.first-1) && (m<it->second.second-1)) {
					std::cout << "FIRST!\n";
					srf_it->Elements.push_back(it->first);
				} else if ((m>=it->second.second-1) && (m<it->second.first-1)) {
					std::cout << "SECOND!\n";
					srf_it->Elements.push_back(it->first);
					std::swap(srf_it->Elements.back()[0],srf_it->Elements.back()[1]);
				} */
			}

			//replace Nodes of Geometry by Nodes of individual surface

			const unsigned int undefined_node=std::numeric_limits<unsigned int>::max();
			std::vector<unsigned int> NodeReplacements(Geometry.Nodes.size(),undefined_node);
			unsigned int NodeCounter=0;

			for (unsigned int k=0;k<srf_it->Elements.size();++k) {

				for (int h=0;h<D;h++) {
					unsigned int origin_node=srf_it->Elements[k][h];
					if (NodeReplacements[origin_node]==undefined_node) {
						NodeReplacements[origin_node]=NodeCounter++;
						srf_it->Nodes.push_back(Geometry.Nodes[origin_node]);

					}
					srf_it->Elements[k][h]=NodeReplacements[origin_node];
				}
			}
			++srf_it;
		}
	}

}


#endif /*GEOMETRY_H_*/
