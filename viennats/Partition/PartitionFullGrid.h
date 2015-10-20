#ifndef PARTITIONNONE_H_
#define PARTITIONNONE_H_
/* =========================================================================
   Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.

                            -----------------
                 ViennaTS - The Vienna Topography Simulator
                            -----------------

   Contact:         viennats@iue.tuwien.ac.at

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


#include <cassert>
#include <vector>
#include <cmath>
#include <limits>
#include <fstream>

namespace partition {

	template <class C> class FullGrid {
	
		typedef typename C::IntCoordType 			IntCoordType;
		typedef typename C::IntLinkType 			IntLinkType;
		typedef typename C::CellRefType 			CellRefType;
		static const int D=C::Dimension;
		
		IntCoordType Min_[D];
		IntCoordType Ext_[D];
		IntCoordType Max_[D];
		BoundaryType Boundaries[D];
		
		IntLinkType Increments[D];
		
		void ClearPartition();
		
		std::vector<CellRefType> VoxelList;
		
	public:
		
		class subbox {
			
			friend class FullGrid;
			
			const CellRefType *Voxel;
			IntCoordType Coords[D];
			
		public:	
			
			IntCoordType Extension(int dim) const {
				return 1;
			}
			
			IntCoordType Max(int dim) const {
				return Coords[dim]+1;
			}
					
			IntCoordType Min(int dim) const {
				return Coords[dim];
			}
			
			CellRefType Cell() const {
				return *Voxel;
			}
			
			bool ContainsCell() const {
				return ((Cell()!=C::UndefinedCellRef()));
			}

		};

		IntCoordType Extension(int i) const {
			return Ext_[i];
		}

		IntCoordType Min(int i) const {
			return Min_[i];
		}

		IntCoordType Max(int i) const {
			return Max_[i];
		}
				
		IntLinkType AreaSize(int Direction) const {
			Direction--;
			if (Direction<0) Direction+=D;
			IntLinkType result=Extension(Direction);
			for (int i=0;i<D-2;i++) {
				Direction--;
				if (Direction<0) Direction+=D;				
				result*=Extension(Direction);
			}
			return result;
		}
		
		template <class CellCoordsType, class BncType> void Setup(const CellRefType&, const CellRefType&, const CellCoordsType&, const BncType &,const SplittingType& ,const SurfaceAreaHeuristicLambdaType& );
		
		IntLinkType NumberOfLeaves() const {
			return VoxelList.size();
		}
		
		IntLinkType NumberOfNodes() const {
			return VoxelList.size();
		}
			
		/*IntLinkType UndefinedLink() const {
			return C::UndefinedLink();
		}
		
		CellRefType UndefinedCellRef() const {
			return C::UndefinedCellRef();
		}*/
		
		IntLinkType NumberOfActiveVoxels() const {
			IntLinkType c=0;
			for (IntLinkType i=0;i<VoxelList.size();i++) if (VoxelList[i]!=C::UndefinedCellRef()) c++; 
			return c;
		}
		
		IntLinkType TotalLeavesSurfaceArea() const {
			return NumberOfLeaves()*2*D;
		}
		
		IntLinkType TotalActiveCellsSurfaceArea() const {
			return 2*D*NumberOfActiveVoxels();
		}
		
		IntLinkType BoundingBoxSurfaceArea() const {
			IntLinkType c=0;
			for (int i=0;i<D;++i) {
				if (Boundaries[i]==NONE)c+=AreaSize(i);
			}
			return 2*c;
		}
		
		double AverageNumberOfTraversedLeaves() const {
			return double(TotalLeavesSurfaceArea())/double(BoundingBoxSurfaceArea());
		}
		
		template <class V> subbox Access(V& pos , int Direction, bool DirectionSign) const;
		template <class V> bool GoToNeighborBox(subbox &, V& pos, int Direction, bool DirectionSign) const;
		
		void PrintStatistics(const std::string& FileName) const;

		double get_memory() const {
			return static_cast<double>(sizeof(CellRefType))*VoxelList.size();
		}
				
	};

	
	template <class C> void FullGrid<C>::ClearPartition() {
		VoxelList.clear();
	}

	template <class C> template <class CoordsType, class BncType> void FullGrid<C>::Setup(
            const CellRefType& start,
            const CellRefType& end,
            const CoordsType& CellCoords,
            const BncType& boundaries,
            const SplittingType& splitting_type,
            const SurfaceAreaHeuristicLambdaType& surface_area_heuristic_lambda) {

		
		///delete old partition
		ClearPartition();
		
		for (int i=0;i<D;i++) Boundaries[i]=boundaries[i];

		//determine extensions of domain
		for (int j=0;j<D;j++) {
			Min_[j]=static_cast<typename C::IntCoordType>(CellCoords[0][j]);
			Max_[j]=static_cast<typename C::IntCoordType>(CellCoords[0][j]);
		}
		
		for (typename C::IntLinkType i=1;i<(end-start);i++) {
			for (int j=0;j<D;j++) {
				if (Min_[j]>CellCoords[i][j]) Min_[j]=static_cast<typename C::IntCoordType>(CellCoords[i][j]);
				if (Max_[j]<CellCoords[i][j]) Max_[j]=static_cast<typename C::IntCoordType>(CellCoords[i][j]);
			}
		}
		for (int j=0;j<D;j++) {
			Max_[j]++;
			Ext_[j]=Max_[j]-Min_[j];
		}
		
		
		Increments[0]=1;
		IntLinkType a=Extension(0);
		for (int k=1;k<D;k++) {
			Increments[k]=a;
			a*=Extension(k);
		}
		
		VoxelList.resize(a, C::UndefinedCellRef());
		
		for (IntLinkType i=0;i<end-start;++i) {
			
			IntLinkType c=CellCoords[i][D-1]-Min(D-1);
			
			for (int m=D-2;m>=0;m--) {
				c*=Extension(m);
				c+=CellCoords[i][m]-Min(m);
			}
			VoxelList[c]=i;
		}
		
	}
	
	
	template <class C> template <class V> typename FullGrid<C>::subbox FullGrid<C>::Access(V& pos, int Direction, bool DirectionSign) const {
		
		subbox sb;
		
		for (int i=0;i<D;i++) {
			if (i==Direction) {
				if (DirectionSign) {
					sb.Coords[i]=0;
					pos[i]=0;
				} else {
					sb.Coords[i]=Extension(i)-1;
					pos[i]=1;
				}
			} else {
				sb.Coords[i]=std::min(static_cast<IntLinkType>(pos[i]),static_cast<IntLinkType>(Extension(i)-IntCoordType(1)));
				pos[i]-=sb.Coords[i];
			}
		}

		IntLinkType l=sb.Coords[D-1];
		for (int i=D-2;i>=0;i--) {
			l*=Extension(i);
			l+=sb.Coords[i];
		}
		sb.Voxel=&VoxelList[l];
					
		return sb;
	}
	
	
	
	template <class C> template <class V> bool FullGrid<C>::GoToNeighborBox(subbox &sb, V& pos, int Direction, bool DirectionSign) const {

		if (DirectionSign) {
			if (sb.Coords[Direction]==Extension(Direction)-1) {
				if (Boundaries[Direction]==PERIODIC) {
					sb.Coords[Direction]=0;
					sb.Voxel-=(Increments[Direction]*(Extension(Direction)-1));
				} else {
					return true;
				}
			} else {
				sb.Coords[Direction]++;
				sb.Voxel+=Increments[Direction];
			}
			pos[Direction]=0;
		} else {
			if (sb.Coords[Direction]==0) {
				if (Boundaries[Direction]==PERIODIC) {
					sb.Coords[Direction]=Extension(Direction)-1;
					sb.Voxel+=(Increments[Direction]*(Extension(Direction)-1));
				} else {
					return true;
				}
			} else {
				sb.Coords[Direction]--;
				sb.Voxel-=Increments[Direction];
			}
			pos[Direction]=1;
		}
		return false;
	}
	
	template <class C> void FullGrid<C>::PrintStatistics(const std::string& FileName) const {
		std::ofstream f;
		if(!std::ifstream(FileName.c_str())) {
			f.open(FileName.c_str());
			f << "Number of active Voxels"			<<";";
			f << "Number of Leaves"					<<";";
			f << "Number of Nodes"					<<";";
			f << "TotalLeavesSurfaceArea"			<<";";
			f << "TotalActiveCellsSurfaceArea"		<<";";
			f << "BoundingBoxSurfaceArea"			<<";";
			f << "Average Num. of traversed leaves" <<std::endl;		
		} else {
			f.open(FileName.c_str(),std::ios_base::app);
		}
		f << NumberOfActiveVoxels()				<<";";
		f << NumberOfLeaves()					<<";";
		f << NumberOfNodes()					<<";";
		f << TotalLeavesSurfaceArea()			<<";";
		f << TotalActiveCellsSurfaceArea()		<<";";
		f << BoundingBoxSurfaceArea()			<<";";
		f << AverageNumberOfTraversedLeaves()	<<std::endl;		
		f.close();
	}
	
	
}

#endif /*PARTITIONNONE_H_*/
