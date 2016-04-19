/*
 * DrLithoVTK.h
 *
 *  Created on: 2016-04-18
 *      Author: filipov
 */

#ifndef DRLITHOVTK_H_
#define DRLITHOVTK_H_

#include <vector>
#include <iostream>

using namespace std;

namespace DrL {

	class DrLithoVTK {

	public:
		int D;
		vector<vector<double> > Vertices;
		vector<vector<int> > Elements;
		int NumVertices;
		int NumElements;

		DrLithoVTK(const string& FileName);
		virtual ~DrLithoVTK();

		void ReadFile(const string& Filename);
		void ReadVTK(const string& Filename);
		void ReadVTU(const string& Filename);
		void PrintMaskVTK(const string& FileName);
		double Max(const int& dir);
		double Min(const int& dir);
	};

	class DrLithoWafer {

	public:
		int D;
		vector<double> Vertices;
		vector<int> Elements;
		int NumVertices;
		int NumElements;

		void PrintVTK(const string& FileName);
		void PrintWaferVTK(const char* FileName);
		DrLithoWafer(DrLithoVTK & mask);
		virtual ~DrLithoWafer();
	};
}

#endif /* DRLITHOVTK_H_ */
