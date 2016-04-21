/*
 * DrLithoVTK.cpp
 *
 *  Created on: 2016-04-18
 *      Author: filipov
 */

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <limits>
#include <sstream>
#include "DrLithoVTK.h"

namespace DrL {


	DrLithoVTK::DrLithoVTK(const string& FileName) {
		ReadFile(FileName);
	}
	
	DrLithoVTK::~DrLithoVTK() {
		Vertices.clear();
		Elements.clear();
	}

	void DrLithoVTK::ReadFile(const string& FileName) {

		if (FileName.substr(FileName.find_last_of(".") + 1) == "vtu") {
			ReadVTU(FileName);
		} else if (FileName.substr(FileName.find_last_of(".") + 1) == "vtk") {
			ReadVTK(FileName);
		}
	}

	void DrLithoVTK::ReadVTK(const string& FileName) {

		fstream mask_file(FileName.c_str());

		string line;
		for (int i=0;i<5;i++) getline(mask_file,line);
		double temp_vertex;
		NumVertices=atoi(&line[line.find(" ")+1]);
		Vertices.resize(NumVertices);
		for (int i=0;i<NumVertices;i++) {
			for (int j=0;j<3;j++) {
				mask_file >> temp_vertex;
				Vertices[i].push_back(temp_vertex);
			}
		}
		
		for (int i=0;i<2;i++) getline(mask_file,line);
		NumElements=atoi(&line[line.find(" ")+1]);
		Elements.resize(NumElements);

		int ElemSize;
		int temp_element;
		mask_file >> ElemSize;
		D=3;
		
		for (int i=0;i<NumElements;i++) {
			for (int j=0;j<ElemSize;j++) {
				mask_file >> temp_element;
				Elements[i].push_back(temp_element);
			}
			mask_file >> temp_element;
		}
		mask_file.close();
	}

	void DrLithoVTK::ReadVTU(const string& FileName) {

		fstream mask_file(FileName.c_str());

		string line;
		for (int i=0;i<4;i++) getline(mask_file,line);
		double temp_vertex;
		NumVertices=atoi(&line[line.find("NumberOfPoints")+16]);
		Vertices.resize(NumVertices);
		NumElements=atoi(&line[line.find("NumberOfCells")+15]);
		Elements.resize(NumElements);

		for (int i=0;i<2;i++) getline(mask_file,line);
		for (int i=0;i<NumVertices;i++) {
			for (int j=0;j<3;j++) {
				mask_file >> temp_vertex;
				Vertices[i].push_back(temp_vertex);
			}
		}
		
		for (int i=0;i<6;i++) getline(mask_file,line);

		int ElemSize=(FileName.find("filled")<FileName.npos)?4:3;
		int temp_element;
		D=3;
		
		for (int i=0;i<NumElements;i++) {
			for (int j=0;j<ElemSize;j++) {
				mask_file >> temp_element;
				Elements[i].push_back(temp_element);
			}
		}
		mask_file.close();
	}

	void DrLithoVTK::PrintMaskVTK(const string& FileName) {

		std::ostringstream oss;
		oss << FileName.substr(0, FileName.find("resist")) << "mask";
		if (FileName.find("filled")<FileName.npos) oss << "_volume";
		else oss << "_boundary";
		oss << ".vtk";
		std::ofstream f(oss.str().c_str());
		f << "# vtk DataFile Version 3.0\n";
		f << "vtk output\n";
		f << "ASCII\n";
		f << "DATASET UNSTRUCTURED_GRID\n";
		f << "POINTS " << NumVertices << " float\n";
		for (int i=0;i<NumVertices;i++) {
			for (int j=0;j<3;j++) {
				f << Vertices[i][j] << " ";
			}
			f << "\n";
		}
		f << "CELLS " << NumElements << " " << 5*NumElements << "\n";
		for (int i=0;i<NumElements;i++) {
			f << Elements[i].size() << " ";
			for (int j=0;j<Elements[i].size();j++)
				f << Elements[i][j] << " ";
			f << "\n";
		}
		f << "CELL_TYPES " << NumElements << "\n";
		for (int i=0;i<NumElements;i++) {
			f << "10\n";
		}
		f << "CELL_DATA " << NumElements << "\n";
		f << "SCALARS Material int 1\n";
		f << "LOOKUP_TABLE default\n";
		for (int i=0;i<NumElements;i++) {
			f << "1\n";
		}
		f.close();
	}


	double DrLithoVTK::Max(const int& dir) {
		double maximum=std::numeric_limits<int>::min();
		for (int i=0;i<NumVertices;i++) {
			maximum=std::max(maximum,Vertices[i][dir]);
		}
		return maximum;
	}

	double DrLithoVTK::Min(const int& dir) {
		double minimum=std::numeric_limits<int>::max();
		for (int i=0;i<NumVertices;i++) {
			minimum=std::min(minimum,Vertices[i][dir]);
		}
		return minimum;
	}

	DrLithoWafer::DrLithoWafer(DrLithoVTK & Mask) {
		D=Mask.D;

		double extents[D][2]; // 0 is min, 1 is max
		for (int i=0;i<D;i++) {
			extents[i][0]=2*Mask.Min(i)-Mask.Max(i); //Min extents
			extents[i][1]=2*Mask.Max(i)-Mask.Min(i); //Max extents
		}
		// In the z-direction:
		extents[D-1][1]=Mask.Min(D-1);
			
		for (int i=0;i<2;i++) {
			for (int j=0;j<2;j++) {
				for (int k=0;k<2;k++) {
					Vertices.push_back(extents[0][i]);
					Vertices.push_back(extents[1][j]);
					Vertices.push_back(extents[2][k]);
				}
			}
		}
		NumVertices=Vertices.size()/3;

		NumElements=3*D-4;
		int temp_elements[] = {0,1,2,4,5,1,7,4,3,7,1,2,6,4,2,7,4,1,7,2};
		for (int i=0;i<4*NumElements;i++) {
			Elements.push_back(temp_elements[i]);
		}
	}

	DrLithoWafer::~DrLithoWafer() {
		Vertices.clear();
		Elements.clear();
	}

	void DrLithoWafer::PrintVTK(const string& FileName) {

		std::ostringstream oss;
		oss << FileName.substr(0, FileName.find("resist")) << "wafer_volume.vtk";
		PrintWaferVTK(oss.str().c_str());
	}

	void DrLithoWafer::PrintWaferVTK(const char * FileName) {

		std::ofstream f(FileName);

		f << "# vtk DataFile Version 3.0\n";
		f << "vtk output\n";
		f << "ASCII\n";
		f << "DATASET UNSTRUCTURED_GRID\n";
		f << "POINTS " << NumVertices << " float\n";
		for (int i=0;i<NumVertices;i++) {
			for (int j=0;j<3;j++) {
				f << Vertices[3*i+j] << " ";
			}
			f << "\n";
		}
		f << "CELLS " << NumElements << " " << 5*NumElements << "\n";
		for (int i=0;i<NumElements;i++) {
			f << "4 ";
			for (int j=0;j<D+1;j++)
				f << Elements[4*i+j] << " ";
			f << "\n";
		}
		f << "CELL_TYPES " << NumElements << "\n";
		for (int i=0;i<NumElements;i++) {
			f << "10\n";
		}
		f << "CELL_DATA " << NumElements << "\n";
		f << "SCALARS Material int 1\n";
		f << "LOOKUP_TABLE default\n";
		for (int i=0;i<NumElements;i++) {
			f << "1\n";
		}
		f.close();
	}	
}
