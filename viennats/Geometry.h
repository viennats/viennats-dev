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
#include "message.h"
#include <stdexcept>

#include <vtkIdList.h>
#include <vtkSmartPointer.h>
#include <vtkDataArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLPolyDataReader.h>

#ifdef USE_HDF5
#include "HDF.h"
#endif


///  Includes all Geometry In- and Output-related functions.
namespace geometry {

  // keeps how many nodes to expect for each vtk cell_type
  // 0 means all numbers are allowed
  // first element is a dummy
  unsigned vtk_nodes_for_cell_type[15] = {0, 1, 0, 2, 0, 3, 0, 0, 4, 4, 4, 8, 8, 6, 5};

  template <int D> class geometry {
  public:

    std::vector<lvlset::vec<double, D> > Nodes;
    std::vector<lvlset::vec<unsigned int, D+1> > Elements;
    std::vector<int> Materials;

    lvlset::vec<double, D> Min,Max; //for bounding box

    geometry()  {};

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
      //std::cout << "Reading STR file!" << std::endl;


      std::ifstream reader(FileName.c_str(), std::ios::binary);

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

        for (int i=0;i<number_of_tetrahedrons;i++)
        {
          reader >> dump;
          reader >> dump;
          reader >> region_id;

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

    #ifdef USE_HDF5

    void ReadTDR(  const std::string& FileName,
      double scale,//=1.,
      std::vector<int> InputTransformationDirections,//=std::vector<int>(),
      std::vector<bool> InputTransformationSigns,//=std::vector<bool>(),
      bool change_input_parity,//=false,
      std::vector<double> shift//=std::vector<double>()
    ) {

      H5File* file = new H5File(FileName.c_str(), H5F_ACC_RDWR);

      tdr_geometry geometry;
      geometry.read_collection(file->openGroup("collection"));

      if (D!=geometry.dim) msg::print_error("Dimension in parameters file does not match geometry!");
      delete file;

      Nodes.resize(geometry.nvertices);
      for (unsigned int i=0;i<geometry.nvertices;i++){
        double coords[D];

        for (int j=0;j<D;j++) coords[j]=geometry.vertex[geometry.dim*i+j];

        for (int j=0;j<D;j++) {
          Nodes[i][j]=coords[InputTransformationDirections[j]];

          int shift_size=shift.size();
          if (shift_size>j) Nodes[i][j]+=shift[j];
          if (InputTransformationSigns[j]) Nodes[i][j]=-Nodes[i][j];
          Nodes[i][j]*=scale;
        }
      }

      int num_elems=0;
      for (map<string,region_t>::iterator S=geometry.region.begin(); S!=geometry.region.end(); S++) num_elems+=S->second.nelements;

      Elements.resize(num_elems);
      Materials.resize(num_elems);

      int i=0;
      for (map<string,region_t>::iterator S=geometry.region.begin(); S!=geometry.region.end(); S++) {

        vector<vector<int> > &elements = S->second.elements;
        for (vector<vector<int> >::iterator E=elements.begin(); E!=elements.end(); E++) {
          std::vector<int> e=*E;
          for (int j=0;j<D+1;j++) Elements[i][j]=e[j];
          Materials[i]=2-S->second.regnr;
          i++;
        }
      }
    }
    #endif

    void ReadGRD(  const std::string& FileName,
      double scale,//=1.,
      std::vector<int> InputTransformationDirections,//=std::vector<int>(),
      std::vector<bool> InputTransformationSigns,//=std::vector<bool>(),
      bool change_input_parity,//=false,
      std::vector<double> shift//=std::vector<double>()
    ) {
      std::ifstream f(FileName.c_str());

      if (!f) msg::print_error("Failed reading geometry file!");

      std::string c;

      //read nodes
      std::getline(f,c);

      int num_nodes;
      int num_elems;
      int num_mater;

      std::vector<lvlset::vec<unsigned int, 2> > Edges;
      std::vector<lvlset::vec<unsigned int, 3> > Faces;

      unsigned int num_edges;
      unsigned int num_faces;
      while (c.find("nb_vertices")>=c.npos) std::getline(f,c);
      num_nodes=atoi(&c[14+c.find("nb_vertices")]);
      //      std::cout << "num_nodes = " << num_nodes << "\n";
      while (c.find("nb_edges")>=c.npos) std::getline(f,c);
      num_edges=atoi(&c[14+c.find("nb_edges")]);
      //      std::cout << "num_edges = " << num_edges << "\n";
      while (c.find("nb_faces")>=c.npos) std::getline(f,c);
      num_faces=atoi(&c[14+c.find("nb_faces")]);
      //      std::cout << "num_faces = " << num_faces << "\n";
      while (c.find("nb_elements")>=c.npos) std::getline(f,c);
      num_elems=atoi(&c[14+c.find("nb_elements")]);
      //      std::cout << "num_elems = " << num_elems << "\n";
      Elements.resize(num_elems);
      Materials.resize(num_elems);
      while (c.find("nb_regions")>=c.npos) std::getline(f,c);
      num_mater=atoi(&c[14+c.find("nb_regions")]);
      //      std::cout << "num_mater = " << num_mater << "\n";

      while (c.find("Vertices")>=c.npos) std::getline(f,c);
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

      while (c.find("Edges")>=c.npos) std::getline(f,c);
      Edges.resize(num_edges);
      for (unsigned int i=0;i<num_edges;i++) {
        for (int j=0;j<2;j++) {
          f>>Edges[i][j];
        }
      }

      if (num_faces>0) while (c.find("Faces")>=c.npos) std::getline(f,c);
      Faces.resize(num_faces);
      unsigned int faces_tmp;
      for (unsigned int i=0;i<num_faces;i++) {
        for (int j=-1;j<D;j++) {
          if (j==-1) f>>faces_tmp;
          else {
            f>>faces_tmp;
            if (faces_tmp>num_edges) {
              faces_tmp=UINT_MAX-faces_tmp;
              Faces[i][j]=Edges[faces_tmp][1];
            } else {
              Faces[i][j]=Edges[faces_tmp][0];
            }
          }
        }
      }
      while (c.find("Elements")>=c.npos) std::getline(f,c);

      unsigned int compare=(D==2)?num_edges:num_faces;
      unsigned int elem_tmp;

      for (int i=0;i<num_elems;i++) {
        for (int j=-1;j<D+1;j++) {
          f>>elem_tmp;
          int write=0;
          if (elem_tmp>compare) {
            elem_tmp=UINT_MAX-elem_tmp;
            write=1;
          }
          //---------- 2D ----------
          if ((D==2)&&(j!=-1)) {
            Elements[i][j]=Edges[elem_tmp][write];
          }
          //---------- 3D ----------
          if ((D==3)&&(j==0)) {
            for (int k=0;k<D;k++) Elements[i][j+k]=Faces[elem_tmp][k];
          } else if ((D==3)&&(j==1)) {
            int k=0;
            do {
              Elements[i][3]=Faces[elem_tmp][k++];
            } while( (Elements[i][3]==Elements[i][0]) || (Elements[i][3]==Elements[i][1]) || (Elements[i][3]==Elements[i][2]) );
          }

        }
      }

      unsigned int elem_in_mat;
      int current;
      for (int i=num_mater;i>0;i--){
        while (c.find("Region")>=c.npos) std::getline(f,c);
        while (c.find("Elements")>=c.npos) std::getline(f,c);
        elem_in_mat=atoi(&c[10+c.find("Elements")]);
        for (unsigned int j=0;j<elem_in_mat;j++) {
          f>>current;
          Materials[current]=i;
        }
      }
      f.close();
    }

    void ReadDX(  const std::string& FileName,
      double scale,//=1.,
      std::vector<int> InputTransformationDirections,//=std::vector<int>(),
      std::vector<bool> InputTransformationSigns,//=std::vector<bool>(),
      bool change_input_parity,//=false,
      std::vector<double> shift//=std::vector<double>()
    ) {
      std::ifstream f(FileName.c_str());

      if (!f) msg::print_error("Failed reading geometry file!");

      std::string c;

      //read nodes
      std::getline(f,c);

      int num_nodes;
      int num_elems;

      //std::cout << "DX\n";
      num_nodes=atoi(&c[63]);

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
      num_elems=atoi(&c[63]);

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

    void ReadVTK(  const std::string& FileName,
      double scale,//=1.,
      std::vector<int> InputTransformationDirections,//=std::vector<int>(),
      std::vector<bool> InputTransformationSigns,//=std::vector<bool>(),
      bool change_input_parity,//=false,
      std::vector<double> shift//=std::vector<double>()
    ) {
      // open geometry file
      std::ifstream f(FileName.c_str());
      if (!f) msg::print_error("Could not open geometry file!");
      std::string c;

      // Check if geometry is an unstructured grid as required
      while(std::getline(f,c)){
        if(c.find("DATASET") != std::string::npos) break;
      }
      if(c.find("UNSTRUCTURED_GRID") == std::string::npos){
        msg::print_error("DATASET is not an UNSTRUCTURED_GRID!");
      }

      // Find POINTS in file to know number of nodes to read in
      while(std::getline(f,c)){
        if(c.find("POINTS") != std::string::npos) break;
      }
      int num_nodes=atoi(&c[c.find(" ")+1]);

      Nodes.resize(num_nodes);

      for (int i=0;i<num_nodes;i++) {
        double coords[3];

        for (int j=0;j<3;j++) f>>coords[j];
        for (int j=0;j<D;j++) {
          Nodes[i][j]=coords[InputTransformationDirections[j]];
          int shift_size=shift.size();
          if (shift_size>j) Nodes[i][j]+=shift[j];//Assign desired shift
          if (InputTransformationSigns[j]) Nodes[i][j]=-Nodes[i][j];//Assign sign transformation, if needed
          Nodes[i][j]*=scale;//Scale the geometry according to parameters file
        }
      }

      while(std::getline(f,c)){
        if(c.find("CELLS") == 0) break;
      }


      int num_elems=atoi(&c[c.find(" ")+1]);

      std::ifstream f_ct(FileName.c_str());   // stream to read cell CELL_TYPES
      std::ifstream f_m(FileName.c_str());  //stream for material numbers if they exist

      // advance to cell types and check if there are the right number
      while(std::getline(f_ct, c)){
        if(c.find("CELL_TYPES") == 0) break;
      }
      int num_cell_types = atoi(&c[c.find(" ")+1]);
      // need a cell_type for each cell
      if(num_elems != num_cell_types){
        msg::print_error("Corrupt input geometry! Number of CELLS and CELL_TYPES is different!");
      }

      bool is_material = true;
      // advance to material if it is specified
      while(std::getline(f_m,c)){
        if(c.find("CELL_DATA") != std::string::npos){
          std::getline(f_m,c);
          if((c.find("SCALARS material") != std::string::npos)||(c.find("SCALARS Material") != std::string::npos)){
            std::getline(f_m,c);
            break;
          }
        }
      }
      if(f_m.eof()){
        is_material = false;
      }

      Elements.clear();
      Elements.reserve(num_elems);

      Materials.clear();

      unsigned elems_fake;
      unsigned cell_type;
      unsigned cell_material;
      for (int i=0;i<num_elems;i++) {
        f >> elems_fake;
        f_ct >> cell_type;
        if(is_material) f_m >> cell_material;
        else cell_material = 1; // if there are no materials specified make all the same

        lvlset::vec<unsigned int, D+1> elem = lvlset::vec<unsigned int, D+1>();

        // check if the correct number of nodes for cell_type is given
        unsigned number_nodes = vtk_nodes_for_cell_type[cell_type];
        if(number_nodes == elems_fake || number_nodes == 0){
          // check for different types to subdivide them into supported types
          switch (cell_type) {
            case 5: //triangle for 2D
            case 10: //tetra for 3D
              for(unsigned j=0; j<number_nodes; ++j){
                f >> elem[j];
              }
              Elements.push_back(elem);
              Materials.push_back(cell_material);
              break;

            case 9:       //this is a quad, so just plit it into two triangles
              for(unsigned j=0; j<3; ++j){
                f >> elem[j];
              }
              Elements.push_back(elem);   //push the first three nodes as a triangle
              Materials.push_back(cell_material);

              f >> elem[1]; //replace middle element to create other triangle
              Elements.push_back(elem);
              Materials.push_back(cell_material);
              break;

            default:
              std::ostringstream oss;
              oss << "VTK Cell type " << cell_type << " is not supported. Cell ignored..." << std::endl;
              msg::print_warning(oss.str());
          }
        }else{
          msg::print_wrong_cell_type(D+1, number_nodes);
          //ignore rest of lines
          f.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
      }

      f_ct.close();
      f_m.close();
      f.close();

    }

    void ReadVTU(const std::string& FileName,
      double scale,//=1.,
      std::vector<int> InputTransformationDirections,//=std::vector<int>(),
      std::vector<bool> InputTransformationSigns,//=std::vector<bool>(),
      bool change_input_parity,//=false,
      std::vector<double> shift//=std::vector<double>()
    ){

      vtkSmartPointer<vtkXMLUnstructuredGridReader> greader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
      greader->SetFileName(FileName.c_str());
      greader->Update();

      vtkSmartPointer<vtkUnstructuredGrid> ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
      ugrid = greader->GetOutput();

      // get all points
      Nodes.resize(ugrid->GetNumberOfPoints());
      for(unsigned i=0; i<Nodes.size(); ++i){
        double p[3];
        ugrid->GetPoint(i, p);
        for(unsigned j=0; j<D; ++j){
          Nodes[i][j]=p[InputTransformationDirections[j]];
          if(shift.size()>j) Nodes[i][j]+=shift[j];//Assign desired shift
          if(InputTransformationSigns[j]) Nodes[i][j]=-Nodes[i][j];//Assign sign transformation, if needed
          Nodes[i][j]*=scale;//Scale the geometry according to parameters file
        }
      }

      // get all cells
      Elements.resize(ugrid->GetNumberOfCells());
      for(unsigned i=0; i<Elements.size(); ++i){
        vtkIdList* pointList = vtkIdList::New();
        ugrid->GetCellPoints(i, pointList);
        lvlset::vec<unsigned int, D+1> elem;

        for(unsigned j=0; j<D+1; ++j){
          elem[j]=pointList->GetId(j);
        }
        Elements[i] = elem;
      }

      // get materials
      vtkSmartPointer<vtkCellData> cellData = vtkSmartPointer<vtkCellData>::New();
      cellData = ugrid->GetCellData();

      int arrayIndex;
      vtkDataArray* matArray = cellData->GetArray("Material", arrayIndex);
      Materials.reserve(Elements.size()); // must be the same number
      if(arrayIndex>=0){  // if array exists
        for(unsigned i=0; i<matArray->GetNumberOfTuples(); ++i){
          Materials.push_back(matArray->GetTuple1(i));
        }
      }else{  // if no material specified, use the same for all elements
        for(unsigned i=0; i<Elements.size(); ++i){
          Materials.push_back(1);
        }
      }
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

      //Check if the file is of format .tdr, .grd, .dx, or .vtk
      std::string GeometryFile=FileName.c_str();
      #ifdef USE_HDF5
      if (GeometryFile.find(".tdr") == (GeometryFile.size()-4)) {
        ReadTDR(GeometryFile, scale, InputTransformationDirections, InputTransformationSigns, change_input_parity, shift);
      } else
      #endif
      if (FileName.substr(FileName.find_last_of(".") + 1) == "str") {
        ReadSTR(FileName, scale, InputTransformationDirections, InputTransformationSigns, change_input_parity, shift, ignore_materials);
      } else if (GeometryFile.find(".grd") == (GeometryFile.size()-4)) {
        ReadGRD(GeometryFile, scale, InputTransformationDirections, InputTransformationSigns, change_input_parity, shift);
      } else if (GeometryFile.find(".dx") == (GeometryFile.size()-3)) {
        ReadDX(GeometryFile, scale, InputTransformationDirections, InputTransformationSigns, change_input_parity, shift);
      } else if (GeometryFile.find(".vtk") == (GeometryFile.size()-4)) {
        ReadVTK(GeometryFile, scale, InputTransformationDirections, InputTransformationSigns, change_input_parity, shift);
      } else if(GeometryFile.find(".vtu") == (GeometryFile.size()-4)){
        ReadVTU(GeometryFile, scale, InputTransformationDirections, InputTransformationSigns, change_input_parity, shift);
      } else {
        #ifdef USE_HDF5
        msg::print_error("This software accepts only STR, TDR, GRD, DX and VTK geometry files!");
        #else
        msg::print_error("This software accepts only STR, GRD, DX and VTK geometry files!");
        #endif
      }

      //map materials - this is the same regardless of the input file format
      if (!MapMaterials.empty()){
        std::vector<lvlset::vec<unsigned int, D+1> > oldElements;
        std::vector<int> oldMaterials;

        std::swap(oldElements, Elements);
        std::swap(oldMaterials, Materials);

        for (unsigned int a=0;a<oldElements.size();++a)
        {
          bool add=true;
          int mat=oldMaterials[a];
          if (mat<static_cast<int>(MapMaterials.size()+1))
          {
            mat=MapMaterials[mat-1];
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

    void WriteVTK(const std::string& FileName) const {
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

    void Write(const std::string& FileName) const {
      WriteVTK(FileName);
    }

    void CalculateExtensions() {
      assert(!Elements.empty());
      Min=Nodes[Elements[0][0]];
      Max=Nodes[Elements[0][0]];

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
    std::vector<int> Materials;

    static constexpr int dimension=D;

    lvlset::vec<double, D> Min,Max; //for bounding box

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

    void ReadVTK(  std::string FileName,
      double scale,
      std::vector<int> InputTransformationDirections,
      std::vector<bool> InputTransformationSigns,
      bool change_input_parity,
      std::vector<double> shift
    ) {
      //Assign desired transformation:
      //-------------------------------------------------------------------------------------------------------------------------
      while(InputTransformationDirections.size()<D) InputTransformationDirections.push_back(InputTransformationDirections.size());        //TODO test if directions are unique
      while(InputTransformationSigns.size()<D) InputTransformationSigns.push_back(false);

      if ((InputTransformationDirections[0]+1)%D!=InputTransformationDirections[1]) change_input_parity=!change_input_parity;

      for (int i=0;i<D;++i) {
        if (InputTransformationSigns[i]) change_input_parity=!change_input_parity;
      }
      //-------------------------------------------------------------------------------------------------------------------------
      std::ifstream f(FileName.c_str());
      if (!f) msg::print_error("Could not open geometry file!");
      std::string c;

      // Check if geometry is an unstructured grid as required
      while(std::getline(f,c)){
        if(c.find("DATASET") != std::string::npos) break;
      }
      if(c.find("UNSTRUCTURED_GRID") == std::string::npos){
        msg::print_error("DATASET is not an UNSTRUCTURED_GRID!");
      }

      // Find POINTS
      while(std::getline(f,c)){
        if(c.find("POINTS") != std::string::npos) break;
      }
      int num_nodes=atoi(&c[c.find(" ")+1]);

      Nodes.resize(num_nodes);

      for (int i=0;i<num_nodes;i++) {
        double coords[3];

        for (int j=0;j<3;j++) f>>coords[j];
        for (int j=0;j<D;j++) {
          Nodes[i][j]=coords[InputTransformationDirections[j]];
          int shift_size=shift.size();
          if (shift_size>j) Nodes[i][j]+=shift[j];//Assign desired shift
          if (InputTransformationSigns[j]) Nodes[i][j]=-Nodes[i][j];//Assign sign transformation, if needed
          Nodes[i][j]*=scale;//Scale the geometry according to parameters file
        }
      }

      while(std::getline(f,c)){
        if(c.find("CELLS") == 0) break;
      }

      int num_elems=atoi(&c[c.find(" ")+1]);

      Elements.clear();
      Elements.reserve(num_elems);

      double elems_fake;
      for (int i=0;i<num_elems;i++) {

        lvlset::vec<unsigned int, D> elem = lvlset::vec<unsigned int, D>();

        // TODO: change the following to read in different types of cells as well
        f >> elems_fake;
        if(elems_fake != D){
          msg::print_wrong_cell_type(D, elems_fake);
          //ignore rest of lines
          f.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        } else{
          for(int j=0; j<D; ++j){
            f >> elem[j];
          }
          Elements.push_back(elem);
        }
      }

      CalculateExtensions();
      f.close();
    }

    void ReadVTP(  std::string FileName,
      double scale,
      std::vector<int> InputTransformationDirections,
      std::vector<bool> InputTransformationSigns,
      bool change_input_parity,
      std::vector<double> shift
    ) {
      // Prepare transformation properties
      while(InputTransformationDirections.size()<D) InputTransformationDirections.push_back(InputTransformationDirections.size());
      while(InputTransformationSigns.size()<D) InputTransformationSigns.push_back(false);

      if ((InputTransformationDirections[0]+1)%D!=InputTransformationDirections[1]) change_input_parity=!change_input_parity;
      for (int i=0;i<D;++i) {
        if (InputTransformationSigns[i]) change_input_parity=!change_input_parity;
      }

      // READ FROM VTP FILE
      vtkSmartPointer<vtkXMLPolyDataReader> pReader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
      pReader->SetFileName(FileName.c_str());
      pReader->Update();

      vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
      polyData = pReader->GetOutput();


      // read points(Nodes)
      Nodes.resize(polyData->GetNumberOfPoints());
      for (unsigned i=0;i<Nodes.size();i++) {
        double coords[3];
        polyData->GetPoint(i, coords);

        for (int j=0;j<D;j++) {
          Nodes[i][j]=coords[InputTransformationDirections[j]];
          int shift_size=shift.size();
          if (shift_size>j) Nodes[i][j]+=shift[j];//Assign desired shift
          if (InputTransformationSigns[j]) Nodes[i][j]=-Nodes[i][j];//Assign sign transformation, if needed
          Nodes[i][j]*=scale;//Scale the geometry according to parameters file
        }
      }

      //read cells (Elements)
      vtkCellArray* cellArray;
      if(D==3){
        Elements.reserve(polyData->GetNumberOfPolys());
        cellArray = polyData->GetPolys();
      }else{
        Elements.reserve(polyData->GetNumberOfLines());
        cellArray = polyData->GetLines();
      }

      cellArray->InitTraversal();
      vtkIdList* pointList = vtkIdList::New();
      while(cellArray->GetNextCell(pointList)){
        lvlset::vec<unsigned int, D> elem;
        for(unsigned j=0; j<D; ++j){
          elem[j]=pointList->GetId(j);
        }
        Elements.push_back(elem);
      }

      // get materials
      vtkSmartPointer<vtkCellData> cellData = vtkSmartPointer<vtkCellData>::New();
      cellData = polyData->GetCellData();

      int arrayIndex;
      vtkDataArray* matArray = cellData->GetArray("Material", arrayIndex);
      Materials.reserve(Elements.size()); // must be the same number
      if(arrayIndex>=0){  // if array exists
        for(unsigned i=0; i<matArray->GetNumberOfTuples(); ++i){
          Materials.push_back(matArray->GetTuple1(i));
        }
      }else{  // if no material specified, use the same for all elements
        for(unsigned i=0; i<Elements.size(); ++i){
          Materials.push_back(1);
        }
      }


      CalculateExtensions();
    }

    void CalculateExtensions() {
      assert(!Elements.empty());
      Min=Nodes[Elements[0][0]];
      Max=Nodes[Elements[0][0]];

      for (unsigned int a=0;a<Elements.size();++a) {
        for (unsigned int i=0;i<D;i++) {
          Min=lvlset::Min(Min,Nodes[Elements[a][i]]);
          Max=lvlset::Max(Max,Nodes[Elements[a][i]]);
        }
      }
    }
  };


  template<class MaterialsType> void GetMaterialNumbers(
      const MaterialsType& materials,
      MaterialsType& discreteMaterials){
    for(auto it=materials.begin(); it!=materials.end(); ++it){
      // if material is not yet in discreteMaterials
      if(std::find(discreteMaterials.begin(), discreteMaterials.end(), *it)==discreteMaterials.end()){
        discreteMaterials.push_back(*it);
      }
    }
    std::sort(discreteMaterials.begin(), discreteMaterials.end());
  }


  template <int D, class SurfacesType> void TransformGeometryToSurfaces(
    const geometry<D>& Geometry,
    SurfacesType &Surfaces,
    std::bitset<2*D> remove_flags,
    double eps,
    bool report_import_errors) {
      // get the unique material numbers for explicit booling
      std::vector<int> materialInts;
      GetMaterialNumbers(Geometry.Materials, materialInts);
      std::vector<unsigned> materialNumbers(materialInts.begin(), materialInts.end());

      //determine maximum number of materials
      unsigned max_mat = materialNumbers.back();


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

          //if (is_open_boundary_negative) flags.reset(open_boundary_direction);
          //else flags.reset(open_boundary_direction+D);

          if (flags.any()) continue;

          tmp.sort();

          lvlset::vec<double,D> pts[D+1];
          for (int k=0;k<D;k++) pts[k]=Geometry.Nodes[tmp[k]];

          pts[D]=Geometry.Nodes[Geometry.Elements[i][(j+D)%(D+1)]];

          typename triangle_map::iterator it=Triangles.lower_bound(tmp);
          if ((it!=Triangles.end()) && (it->first==tmp)) {
            if (lvlset::Orientation(pts)) {
              if (report_import_errors && it->second.second!=max_mat+1){
                std::ostringstream oss;
                oss << "Coinciding triangles with same orientation in Element: " << i << std::endl;
                msg::print_error(oss.str());
              }
              it->second.second=Geometry.Materials[i];
            } else {
              if (report_import_errors && it->second.first!=max_mat+1){
                std::ostringstream oss;
                oss << "Coinciding triangles with same orientation in Element: " << i << std::endl;
                msg::print_error(oss.str());
              }
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


      Surfaces.resize(materialNumbers.size());

      //for all materials/for each surface
      typename SurfacesType::iterator srf_it=Surfaces.begin();
      for (auto matIt=materialNumbers.begin(); matIt!=materialNumbers.end(); ++matIt) {
        for (typename triangle_map::iterator it=Triangles.begin();it!=Triangles.end();++it) {
          if (((*matIt)>=it->second.first) && ((*matIt)<it->second.second)) {
            srf_it->Elements.push_back(it->first);
          } else if (((*matIt)>=it->second.second) && ((*matIt)<it->second.first)) {
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


    template<int D, class SurfacesType> void TransformHullToSurfaces(
        const surface<D>& Geometry,
        SurfacesType &Surfaces,
        std::bitset<2*D> remove_flags,
        double eps,
        bool report_import_errors){

      // get the unique material numbers for explicit booling
      std::vector<int> materialNumbers;
      GetMaterialNumbers(Geometry.Materials, materialNumbers);

      // sort all elements based on materials, so they are in ascending order
      typedef std::pair< lvlset::vec<unsigned int, D>, int> materialElementType;
      std::vector<materialElementType> materialElements;

      materialElements.reserve(Geometry.Elements.size());
      for(unsigned i=0; i<Geometry.Elements.size(); ++i){
        materialElements.push_back(std::make_pair(Geometry.Elements[i], Geometry.Materials[i]));
      }

      // elements will now be sorted based on material number
      std::sort(materialElements.begin(), materialElements.end(), [](const materialElementType& i, const materialElementType& j) { return i.second < j.second; });


      // cycle through all available materials and use all triangles of lower materials as well
      for(auto materialIt=materialNumbers.begin(); materialIt!=materialNumbers.end(); ++materialIt){
        // for each material we need a new surface
        Surfaces.push_back(typename SurfacesType::value_type());
        typename SurfacesType::value_type searchSurface;
        // get position in element vector based on material
        // upper_bound returns iterator to first element greater than current material
        // therefore elemIterator is essentially our .end() element
        // same goes for elemStart as .begin() element
        auto elemStart = std::lower_bound(materialElements.begin(), materialElements.end(), *materialIt, [](const materialElementType& i, const int& j) { return i.second < j; });
        auto elemEnd = std::upper_bound(materialElements.begin(), materialElements.end(), *materialIt, [](const int& i, const materialElementType& j) { return i < j.second; });

        // go through all elements and remove duplicates and boundary elements
        for(auto elemIt=elemStart; elemIt!=elemEnd; ++elemIt){
          lvlset::vec<unsigned, D> element;
          for(unsigned i=0; i<D; ++i) element[i]=(elemIt->first)[i];

          std::bitset<2*D> flags;
          flags.set();
          //if triangle at border skip
          for (int k=0;k<D;k++) {
            for (int l=0;l<D;l++) {
              if (Geometry.Nodes[element[k]][l]<Geometry.Max[l]-eps) {
                flags.reset(l+D);
              }
              if (Geometry.Nodes[element[k]][l]>Geometry.Min[l]+eps) {
                flags.reset(l);
              }
            }
          }
          flags &=remove_flags;
          if (flags.any()) continue;

          element.sort();
          // if triangle already exists, check if there already is a degenerate one
          auto degenerateIt = std::find(searchSurface.Elements.begin(), searchSurface.Elements.end(), element);
          if(degenerateIt != searchSurface.Elements.end()){ //degenerate found
            searchSurface.Elements.erase(degenerateIt); // remove element
            Surfaces.back().Elements.erase(Surfaces.back().Elements.begin()+std::distance(searchSurface.Elements.begin(), degenerateIt)); // remove actual element
            continue; // ignore current element as well
          }else{
            searchSurface.Elements.push_back(element); // push element if it does not exist yet
            Surfaces.back().Elements.push_back(elemIt->first); // push element in actual surfaces
          }
        }

        // Take only the required nodes for the new surface
        const unsigned undefinedNode=std::numeric_limits<unsigned int>::max();
        std::vector<unsigned> NodeReplacements(Geometry.Nodes.size(),undefinedNode);
        unsigned NodeCounter=0;
        // go over all elements and replace their node ids with the new ones
        for (unsigned k=0;k<Surfaces.back().Elements.size();++k) {
          for (int h=0;h<D;h++) {
            unsigned origin_node=Surfaces.back().Elements[k][h];
            if (NodeReplacements[origin_node]==undefinedNode) {
              NodeReplacements[origin_node]=NodeCounter++;
              Surfaces.back().Nodes.push_back(Geometry.Nodes[origin_node]);
            }
            Surfaces.back().Elements[k][h]=NodeReplacements[origin_node];
          }
        }
      }
    }

    //Reads in surfaces and transforms it to a levelset
    template<int D, class GridTraitsType, class ParameterType, class LevelSetType>
    void import_levelsets_from_surface(GridTraitsType& GridProperties, lvlset::grid_type<GridTraitsType>& grid,
                                      ParameterType& p, std::list<LevelSetType>& LevelSets)
    {
      int grid_min[D]={ };
      int grid_max[D]={ };
      std::list< surface<D> > surfaces;


      std::cout << "The geometry consists of " << p.geometry_files.size() <<" input surfaces. \n";

      //!surface.ReadVTK(...) reads surface file/s and modifies it/them according to the user-set parameters
      for(unsigned i=0; i<p.geometry_files.size(); ++i){
        surfaces.push_back(surface<D>());
        msg::print_start("Read surface input file " + p.geometry_files[i] + "...");
        if(p.geometry_files[i].find(".vtk") == (p.geometry_files[i].size()-4)){
          surfaces.back().ReadVTK(p.geometry_files[i], p.input_scale, p.input_transformation,
                    p.input_transformation_signs, p.change_input_parity, p.input_shift);
        }else if(p.geometry_files[i].find(".vtp") == (p.geometry_files[i].size()-4)){
          surfaces.back().ReadVTP(p.geometry_files[i], p.input_scale, p.input_transformation,
                    p.input_transformation_signs, p.change_input_parity, p.input_shift);
        }else{
          msg::print_error("Unknown filetype of " + p.geometry_files[i]);
        }


        for (int h = 0; h < D; ++h) {
          grid_min[h] = std::min(grid_min[h],int(std::ceil(surfaces.back().Min[h] / p.grid_delta - p.snap_to_boundary_eps)));
          grid_max[h] = std::max(grid_max[h],int(std::floor(surfaces.back().Max[h] / p.grid_delta + p.snap_to_boundary_eps)));
        }
        msg::print_done();
#ifdef VERBOSE
        std::cout << "min = " << (surfaces.back().Min) << "   " << "max = " << (surfaces.back().Max) << std::endl;
        std::cout << "min = " << (surfaces.back().Min / p.grid_delta) << "   " << "max = " << (surfaces.back().Max / p.grid_delta) << std::endl;
#endif
      }

      //!Determine boundary conditions for level set domain
      lvlset::boundary_type bnc[D];
      for (int hh = 0; hh < D; ++hh) {
        if (p.boundary_conditions[hh].min == bnc::PERIODIC_BOUNDARY &&
            p.boundary_conditions[hh].max == bnc::PERIODIC_BOUNDARY) {
              bnc[hh] = lvlset::PERIODIC_BOUNDARY;
        } else if(p.boundary_conditions[hh].min == bnc::INFINITE_BOUNDARY &&
                  p.boundary_conditions[hh].max == bnc::INFINITE_BOUNDARY) {
              bnc[hh] = lvlset::INFINITE_BOUNDARY;
        } else if (p.boundary_conditions[hh].min == bnc::INFINITE_BOUNDARY) {
              bnc[hh] = lvlset::NEG_INFINITE_BOUNDARY;
        } else if (p.boundary_conditions[hh].max == bnc::INFINITE_BOUNDARY) {
              bnc[hh] = lvlset::POS_INFINITE_BOUNDARY;
        } else {
              bnc[hh] = lvlset::SYMMETRIC_BOUNDARY;
        }
      }

      //Set the level set GridProperties
      GridProperties = GridTraitsType(grid_min, grid_max, bnc, p.grid_delta);
      //Generate the grid with the GridProperties
      grid = lvlset::grid_type<GridTraitsType>(GridProperties);

      msg::print_start("Distance transformation...");

      //!Initialize the level set with "lvlset::init(...)"
      for(typename std::list< surface<D> >::const_iterator it= surfaces.begin(); it!=surfaces.end(); ++it){
        LevelSets.push_back(LevelSetType(grid));
        lvlset::init(LevelSets.back(), *it, p.report_import_errors);
        LevelSets.back().set_levelset_id();
      }

      msg::print_done();
    }









    template<int D, class GridTraitsType, class ParameterType, class LevelSetType>
    void import_levelsets_from_volume(GridTraitsType& GridProperties, lvlset::grid_type<GridTraitsType>& grid,
                                      ParameterType& p, std::list<LevelSetType>& LevelSets)
    {
      int grid_min[D]={ };
      int grid_max[D]={ };
      //!Read Geometry and populate geometry class
      geometry<D> g;
      // g.Read reads a geometry file and modifies it according to the user-set parameters
      std::cout << "Read geometry input file " << p.geometry_files[0];
      msg::print_start("...");
      g.Read(p.geometry_files[0], p.input_scale, p.input_transformation,
             p.input_transformation_signs, p.change_input_parity, p.material_mapping,
             p.input_shift, p.ignore_materials);
      {
        // output to a .vtk file the modified initial geometry
        std::ostringstream oss;
        oss << p.output_path << "Initial_Volume_Mesh.vtk";
        g.Write(oss.str());
      }
      for (int h = 0; h < D; ++h) {
        grid_min[h] = std::ceil(g.Min[h] / p.grid_delta - p.snap_to_boundary_eps);
        grid_max[h] = std::floor(g.Max[h] / p.grid_delta + p.snap_to_boundary_eps);
      }
    #ifdef VERBOSE
      std::cout << "min = " << (g.Min) << "   " << "max = " << (g.Max) << std::endl;
      std::cout << "min = " << (g.Min / p.grid_delta) << "   " << "max = " << (g.Max / p.grid_delta) << std::endl;
    #endif
      msg::print_done();

      //!Determine boundary conditions for level set domain
      lvlset::boundary_type bnc[D];
      for (int hh = 0; hh < D; ++hh) {
        if (p.boundary_conditions[hh].min == bnc::PERIODIC_BOUNDARY &&
            p.boundary_conditions[hh].max == bnc::PERIODIC_BOUNDARY) {
              bnc[hh] = lvlset::PERIODIC_BOUNDARY;
        } else if(p.boundary_conditions[hh].min == bnc::INFINITE_BOUNDARY &&
                  p.boundary_conditions[hh].max == bnc::INFINITE_BOUNDARY) {
              bnc[hh] = lvlset::INFINITE_BOUNDARY;
        } else if (p.boundary_conditions[hh].min == bnc::INFINITE_BOUNDARY) {
              bnc[hh] = lvlset::NEG_INFINITE_BOUNDARY;
        } else if (p.boundary_conditions[hh].max == bnc::INFINITE_BOUNDARY) {
              bnc[hh] = lvlset::POS_INFINITE_BOUNDARY;
        } else {
              bnc[hh] = lvlset::SYMMETRIC_BOUNDARY;
        }
      }

      //Set the level set GridProperties
      GridProperties = GridTraitsType(grid_min, grid_max, bnc, p.grid_delta);
      //Generate the grid with the GridProperties
      grid = lvlset::grid_type<GridTraitsType>(GridProperties);

      //!Transform the input volume geometry to surfaces and interfaces "TransformGeometryToSurfaces(...)"
      msg::print_start("Extract surface and interfaces...");
      std::bitset<2 * D> remove_flags;

      for (int i = 0; i < D; ++i) {
        if (p.boundary_conditions[i].min == bnc::PERIODIC_BOUNDARY ||
            p.boundary_conditions[i].min == bnc::REFLECTIVE_BOUNDARY ||
            p.boundary_conditions[i].min == bnc::EXTENDED_BOUNDARY) {
              remove_flags.set(i);
        } else if (i == p.open_boundary && !p.open_boundary_negative && p.remove_bottom) {
              remove_flags.set(i);
        }
        if (p.boundary_conditions[i].max == bnc::PERIODIC_BOUNDARY ||
            p.boundary_conditions[i].max == bnc::REFLECTIVE_BOUNDARY ||
            p.boundary_conditions[i].max == bnc::EXTENDED_BOUNDARY) {
              remove_flags.set(i + D);
        } else if (i == p.open_boundary && p.open_boundary_negative && p.remove_bottom) {
              remove_flags.set(i + D);
        }
      }

      typedef std::list<surface<D> > SurfacesType;
      SurfacesType Surfaces;

      //std::cout << "transform to surface\n";
      TransformGeometryToSurfaces(g, Surfaces, remove_flags, p.grid_delta * p.snap_to_boundary_eps, p.report_import_errors);
      msg::print_done();

      msg::print_start("Distance transformation...");

      //!Initialize each level set with "lvlset::init(...)"
      for (typename SurfacesType::const_iterator it = Surfaces.begin(); it != Surfaces.end(); ++it) {
        LevelSets.push_back(LevelSetType(grid));
        lvlset::init(LevelSets.back(), *it, p.report_import_errors);
        LevelSets.back().set_levelset_id();
      }

      msg::print_done();
    }





    /// This function takes a hull mesh, as output by ViennaTS and separates it into different
    /// surfaces for each material, applying explicit wrapping, thus conserving thin layers.
    /// The respective layers are then converted to level sets
    template<int D, class GridTraitsType, class ParameterType, class LevelSetType>
    void import_levelsets_from_hull(GridTraitsType& GridProperties, lvlset::grid_type<GridTraitsType>& grid,
                                      ParameterType& p, std::list<LevelSetType>& LevelSets)
    {

      // Read hull surfaces from geometry file
      surface<D> geometry;
      std::list< surface<D> > surfaces;

      std::cout << "Read geometry input file " << p.geometry_files[0];
      msg::print_start("...");

      // read surface
      geometry.ReadVTP(p.geometry_files[0], p.input_scale, p.input_transformation,
                p.input_transformation_signs, p.change_input_parity, p.input_shift);


      // get information on which triangles should be removed
      std::bitset<2 * D> remove_flags;
      for (int i = 0; i < D; ++i) {
        if (p.boundary_conditions[i].min == bnc::PERIODIC_BOUNDARY ||
            p.boundary_conditions[i].min == bnc::REFLECTIVE_BOUNDARY ||
            p.boundary_conditions[i].min == bnc::EXTENDED_BOUNDARY) {
              remove_flags.set(i);
        } else if (i == p.open_boundary && !p.open_boundary_negative && p.remove_bottom) {
              remove_flags.set(i);
        }
        if (p.boundary_conditions[i].max == bnc::PERIODIC_BOUNDARY ||
            p.boundary_conditions[i].max == bnc::REFLECTIVE_BOUNDARY ||
            p.boundary_conditions[i].max == bnc::EXTENDED_BOUNDARY) {
              remove_flags.set(i + D);
        } else if (i == p.open_boundary && p.open_boundary_negative && p.remove_bottom) {
              remove_flags.set(i + D);
        }
      }

      // now get surfaces from hull mesh
      TransformHullToSurfaces(geometry, surfaces, remove_flags, p.grid_delta * p.snap_to_boundary_eps, p.report_import_errors);


      // EXTRACT GRID AND READ SURFACES INTO LEVELSETS
      int grid_min[D]={ };
      int grid_max[D]={ };

      for (int h = 0; h < D; ++h) {
        grid_min[h] = std::ceil(geometry.Min[h] / p.grid_delta - p.snap_to_boundary_eps);
        grid_max[h] = std::floor(geometry.Max[h] / p.grid_delta + p.snap_to_boundary_eps);
      }
    #ifdef VERBOSE
      std::cout << "min = " << (geometry.Min) << "   " << "max = " << (geometry.Max) << std::endl;
      std::cout << "min = " << (geometry.Min / p.grid_delta) << "   " << "max = " << (geometry.Max / p.grid_delta) << std::endl;
    #endif
      msg::print_done();


      // Determine boundary conditions for level set domain
      lvlset::boundary_type bnc[D];
      for (int hh = 0; hh < D; ++hh) {
        if (p.boundary_conditions[hh].min == bnc::PERIODIC_BOUNDARY &&
            p.boundary_conditions[hh].max == bnc::PERIODIC_BOUNDARY) {
              bnc[hh] = lvlset::PERIODIC_BOUNDARY;
        } else if(p.boundary_conditions[hh].min == bnc::INFINITE_BOUNDARY &&
                  p.boundary_conditions[hh].max == bnc::INFINITE_BOUNDARY) {
              bnc[hh] = lvlset::INFINITE_BOUNDARY;
        } else if (p.boundary_conditions[hh].min == bnc::INFINITE_BOUNDARY) {
              bnc[hh] = lvlset::NEG_INFINITE_BOUNDARY;
        } else if (p.boundary_conditions[hh].max == bnc::INFINITE_BOUNDARY) {
              bnc[hh] = lvlset::POS_INFINITE_BOUNDARY;
        } else {
              bnc[hh] = lvlset::SYMMETRIC_BOUNDARY;
        }
      }

      //Set the level set GridProperties
      GridProperties = GridTraitsType(grid_min, grid_max, bnc, p.grid_delta);
      //Generate the grid with the GridProperties
      grid = lvlset::grid_type<GridTraitsType>(GridProperties);

      msg::print_start("Distance transformation...");

      // Initialize each level set with "lvlset::init(...) and wrap with material below"
      for (typename std::list< surface<D> >::const_iterator it = surfaces.begin(); it != surfaces.end(); ++it) {
        LevelSets.push_back(LevelSetType(grid));
        lvlset::init(LevelSets.back(), *it, p.report_import_errors);
        LevelSets.back().set_levelset_id();
        // need to fix layer wrapping, so wrap with levelset below
        if(it!=surfaces.begin()){
          LevelSets.back().min(*(--(--LevelSets.end())));
        }
      }

      msg::print_done();
    }

  }//Namespace geometry END

  #endif /*GEOMETRY_H_*/
