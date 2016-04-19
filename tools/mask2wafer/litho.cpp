/*
 * litho.cpp
 *
 *  Created on: 2016-04-18
 *      Author: filipov
 */


#include "DrLithoVTK.h"
#include "DrLithoVTK.cpp"

using namespace std;

int main(int argc, char * argv[])
{
	// Populate the data structure
	DrL::DrLithoVTK mask(argv[1]);

	// Print ViennaTS-readable mask VTK
	mask.PrintMaskVTK(argv[1]);

	// Pass mask data to generate wafer
	DrL::DrLithoWafer wafer(mask);

	// Print out the wafer data to VTK
	wafer.PrintVTK(argv[1]);

	return 0;
}
