#ifndef PARTITION_H_
#define PARTITION_H_

/* =========================================================================
   Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.

                            -----------------
                 ViennaTS - The Vienna Topography Simulator
                            -----------------

   Contact:         viennats@iue.tuwien.ac.at

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


#include "levelset.hpp"

///Holds different models for grid point/node access strategies.
namespace partition {

	typedef lvlset::boundary_type BoundaryType;
	const BoundaryType NONE=lvlset::INFINITE_BOUNDARY;
	const BoundaryType PERIODIC=lvlset::PERIODIC_BOUNDARY;
	const BoundaryType SYMMETRIC=lvlset::SYMMETRIC_BOUNDARY;

	enum DataStructureType {NEIGHBOR_LINKS_ARRAYS, FULL_GRID, UP_DOWN_LINKED_TREE};
	enum SplittingType {SPATIAL_MEDIAN, OBJECT_MEDIAN, SURFACE_AREA_HEURISTIC};
	typedef double SurfaceAreaHeuristicLambdaType;
}



#endif /*PARTITION_H_*/
