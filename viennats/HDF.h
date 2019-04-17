#include <map>
#include <typeinfo>
using namespace std;

#include <iostream>
#include <H5Cpp.h>

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

#define mythrow(a) { cerr << a << endl; throw; }
#define mythrow1(a) { cerr << a << endl; throw; }

int read_int(const H5Object &g, const string name)
{
  int i;
  Attribute a=g.openAttribute(name);
  if (a.getTypeClass()!=H5T_INTEGER)
    mythrow("Wrong class in atrribute");
  a.read(a.getDataType(),&i);
  return i;
}
double read_double(const H5Object &g, const string name)
{
  double i;
  Attribute a=g.openAttribute(name);
  if (a.getTypeClass()!=H5T_FLOAT)
    mythrow("Wrong class " << typeid(a.getTypeClass()).name() << " in atrribute");
  a.read(a.getDataType(),&i);
  return i;
}
string read_string(const H5Object &g, const string name)
{
  Attribute a=g.openAttribute(name);
  if (a.getTypeClass()!=H5T_STRING)
    mythrow("Wrong class in atrribute");

  char *buf = new char[a.getDataType().getSize()+1];
//  char buf[a.getDataType().getSize()+1];
  a.read(a.getDataType(),buf);
  buf[a.getDataType().getSize()]='\0';

  return string(buf);
//  delete[] buf;
}

struct dataset_t
{
  string name, quantity, unit;
  unsigned int nvalues;
  double conversion_factor;
  vector<double> values;
};

struct region_t
{
  int regnr;
  string name,material;
  int nelements,npointidx;
  vector<vector<int> > elements;
  map<string,dataset_t> dataset;
};

struct tdr_geometry
{
//  int dim,nvertices,nregions,ndatasets;
  int dim,nregions,ndatasets;
  hsize_t nvertices;
  vector<double> vertex;
  map<string,region_t> region;
  double trans_matrix[9],trans_move[3];

  void read_transformation(const Group &trans)
  {
    const DataSet &A=trans.openDataSet("A");
    const DataSet &b=trans.openDataSet("b");
    A.read( trans_matrix, PredType::NATIVE_DOUBLE);
    b.read( trans_move, PredType::NATIVE_DOUBLE);
  }

  typedef struct coord2_t {
           double x[3];
      } coord2_t;


  void read_vertex(const DataSet &vert)
  {
    const DataSpace &dataspace = vert.getSpace();
//    int rank = dataspace.getSimpleExtentNdims();
    dataspace.getSimpleExtentNdims();
    hsize_t dims[10];
//    int ndims = dataspace.getSimpleExtentDims( dims, NULL);
    dataspace.getSimpleExtentDims( dims, NULL);
    if (nvertices!=dims[0])
      mythrow("nvertices not equal vertices.dim");

    CompType mtype2( sizeof(coord2_t) );
    const H5std_string MEMBER1( "x" );
    const H5std_string MEMBER2( "y" );
    const H5std_string MEMBER3( "z" );
           mtype2.insertMember( "x", HOFFSET(coord2_t, x[0]), PredType::NATIVE_DOUBLE);
    if (dim>1)
    mtype2.insertMember( "y", HOFFSET(coord2_t, x[1]), PredType::NATIVE_DOUBLE);
    if (dim>2)
    mtype2.insertMember( "z", HOFFSET(coord2_t, x[2]), PredType::NATIVE_DOUBLE);
           /*
           * Read two fields x and y from s2 dataset. Fields in the file
           * are found by their names "x_name" and "y_name".
           */
//           coord2_t s2[dims[0]];
           coord2_t *s2 = new coord2_t[dims[0]];
           vert.read( s2, mtype2 );

    for (unsigned int i=0; i<dims[0]; i++)
    {
      vertex.push_back(s2[i].x[0]);
      if (dim>1)
      vertex.push_back(s2[i].x[1]);
      if (dim>2)
      vertex.push_back(s2[i].x[2]);
    }
    delete[] s2;
  }

  void read_elements(region_t &region, const DataSet &elem)
  {
//    cerr << "nelements: " << region.nelements << endl;

    const DataSpace &dataspace = elem.getSpace();
    int rank = dataspace.getSimpleExtentNdims();
    hsize_t dims[10];
    int ndims = dataspace.getSimpleExtentDims( dims, NULL);

    if (rank!=1)
      mythrow("rank of elements in region " << region.name << " is not one");
    if (ndims!=1)
      mythrow("ndims of elements in region " << region.name << " is not one");

//    int el[dims[0]];
    int *el = new int [dims[0]];
    elem.read( el, PredType::NATIVE_INT);

//    int elct=0,eldim;
    unsigned int elct=0;
    while (elct<dims[0])
    {
      region.elements.push_back(vector<int>());
      switch (el[elct++])
      {
        case 1:  for (int i=0; i<2; i++)
            region.elements.back().push_back(el[elct++]);
          break;
        case 2:  for (int i=0; i<3; i++)
            {
              region.elements.back().push_back(el[elct++]);
            }
          break;
        case 5:  for (int i=0; i<4; i++)
            region.elements.back().push_back(el[elct++]);
          break;
        default:mythrow("Element type " << el[elct-1] << " in region " << region.name << " not known");
      }
    }
    delete[] el;
  }

  void read_region(const int regnr, const Group &reg)
  {
    string name0=read_string(reg,"name");
    string name;
    std::size_t i;
    while ((i=name0.find_first_of("_."))!=std::string::npos)
    {
      name=name+name0.substr(0,i);
      name0=name0.substr(i+1);
    }
    name=name+name0;
    string material;
    const int typ=read_int(reg,"type");
    switch (typ)
    {
      case 0: // we have a nD region
        material=read_string(reg,"material");
        break;
      case 1: // we have a nD-1 region (contact... drain, source, ...)
        material="Contact";
        // attribute "part 0", no idea for what
        break;
      case 2: material="Interface"+std::to_string(read_int(reg,"bulk 0"))+std::to_string(read_int(reg,"bulk 1"));
        break;
        return;
    }

    region[name].regnr=regnr;
    region[name].name=name;
    region[name].material=material;
    const int n=reg.getNumObjs();
    for (int i=0; i<n; i++)
    {
      if (reg.getObjnameByIdx(i)=="elements_0")
      {
        const DataSet &ds=reg.openDataSet("elements_0");
        region[name].nelements=read_int(ds,"number of elements");
        read_elements(region[name],ds);
        return;
      }
      if (reg.getObjnameByIdx(i)=="part_0")
      {
        const Group &part=reg.openGroup("part_0");
        region[name].nelements=read_int(part,"number of elements");
        read_elements(region[name],part.openDataSet("elements"));
        return;
      }
    }
    mythrow("No elements found on region " << name);
  }

  region_t &find_region(int regnr)
  {
    for (map<string,region_t>::iterator R=region.begin(); R!=region.end(); R++)
      if (R->second.regnr==regnr)
        return R->second;
    mythrow("Region " << regnr << " not found");
  }

  void read_values(dataset_t &dataset,const DataSet &values)
  {
    const DataSpace &dataspace = values.getSpace();
    dataspace.getSimpleExtentNdims();
    hsize_t dims[10];
    int ndims = dataspace.getSimpleExtentDims( dims, NULL);
    if (dataset.nvalues!=dims[0] || ndims!=1)
      mythrow("Dataset " << dataset.name << " should have " << dataset.nvalues << " values, but has " << dims[0] << " with dimension " << ndims);

//    double v[dims[0]];
    double *v = new double [dims[0]];
    values.read( v, PredType::NATIVE_DOUBLE);
    dataset.values.insert(dataset.values.end(),&v[0],&v[dims[0]]);
    delete[] v;
  }
  void read_dataset(const Group &dataset)
  {
          string name = read_string(dataset,"name");
    if (name.find("Stress")!=name.npos)
      return;

          string quantity = read_string(dataset,"quantity");
          int regnr = read_int(dataset,"region");
          int nvalues=read_int(dataset,"number of values");
          double conversion_factor = read_double(dataset,"conversion factor");
          if (read_int(dataset,"location type")!=0)
    {
      return;
    }
          if (read_int(dataset,"structure type") != 0)
      mythrow("Dataset " << name << " structure type not 0");
          if (read_int(dataset,"value type") != 2)
            mythrow("Dataset " << name << " value type not 2");

    // In this dataset we have a group
    // tag_group_0???
    // units: take only the unit name
    // Dataset values: the actual values

    string unit;
    int n=dataset.getNumObjs(),i;
    for (i=0; i<n; i++)
    {
      cerr << dataset.getObjnameByIdx(i) << endl;
      if (dataset.getObjnameByIdx(i)=="unit")
      {
        const Group &u=dataset.openGroup("unit");
        unit=read_string(u,"name");
        break;
      }
    }
    if (i==n)
    {
      unit=read_string(dataset,"unit:name");
    }

    region_t &region=find_region(regnr);
    region.dataset[name].name=name;
    region.dataset[name].quantity=quantity;
    region.dataset[name].nvalues=nvalues;
    region.dataset[name].conversion_factor=conversion_factor;
    region.dataset[name].unit=unit;

    read_values(region.dataset[name],dataset.openDataSet("values"));
  }

  void read_attribs0(const Group &state)
  {
    int i;
    ndatasets=read_int(state,"number of datasets");

    i=read_int(state,"number of plots");
    if (i!=0)
      mythrow("Numberofplots not equal 0");
    i=read_int(state,"number of string streams");
    if (i!=0)
      mythrow("Numberofstringstreams not equal 0");
    i=read_int(state,"number of xy plots");
    if (i!=0)
      mythrow("Numberofxyplots not equal 0");

    for (i=0; i<ndatasets; i++)
    {
      char a[100];
      sprintf(a,"dataset_%d",i);
      cerr << "Attribute: " << a << endl;
      read_dataset(state.openGroup(a));
    }
  }

  void read_geometry(const Group &geometry)
  {
    int typ,dum;
    char name[100];
    typ=read_int(geometry,"type");
    switch (typ)
    {
      case 0: break; // found a geometry with type 0
      case 1: break;
      default : mythrow(__LINE__ << ": Unknown type " << typ << " in geometry");
    }

    dim=read_int(geometry,"dimension");

    nvertices=read_int(geometry,"number of vertices");

    nregions=read_int(geometry,"number of regions");

    dum=read_int(geometry,"number of states");
    if (dum!=0 && dum!=1)
      mythrow("Number of states not 0 or 1");

    for (int i=0; i<nregions; i++)
    {
      sprintf(name,"region_%d",i);
      const Group &reg=geometry.openGroup(name);
      read_region(i,reg);
    }

    const Group &trans=geometry.openGroup("transformation");
    read_transformation(trans);
    const DataSet &vert=geometry.openDataSet("vertex");
    read_vertex(vert);
    // What do the Units???

    read_attribs0(geometry.openGroup("state_0"));
  }

  void read_collection(const Group &collection)
  {
    int i;

    i=read_int(collection,"number of geometries");
    if (i!=1)
      mythrow("Not only one geometry");
    i=read_int(collection,"number of plots");
    if (i!=0)
      fprintf(stderr,"We have plots, skip them\n");

    read_geometry(collection.openGroup("geometry_0"));
  }
};
