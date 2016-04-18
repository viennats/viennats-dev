#ifndef PARTITION4_UpDownLinkTreeH_
#define PARTITION4_UpDownLinkTreeH_

/* =========================================================================
   Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.

                            -----------------
                 ViennaTS - The Vienna Topography Simulator
                            -----------------

   Contact:         viennats@iue.tuwien.ac.at

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


#include <algorithm>
#include <cassert>
#include <stack>
#include <vector>
#include <cmath>
#include <limits>
#include <fstream>
#include "Partition.h"

namespace partition {

	template <class C> class UpDownLinkTree {
	
		typedef typename C::IntCoordType 			IntCoordType;
		typedef typename C::IntLinkType 			IntLinkType;
		typedef typename C::CellRefType 			CellRefType;
		static const int D=C::Dimension;
	
		class node {
			
		public:
			
			int SplitDirection;
			IntLinkType Parent;
			
			union SCUnion {
				struct SC {
					IntCoordType SplitCoordinate;
					IntLinkType Childs[2];
				} SC;
				struct Cells2 {
					CellRefType Cell;
					IntCoordType Min_[D];
					IntCoordType Ext_[D];
				} Cells2;
				
			} SCUnion;
			
			node() : SplitDirection(-1) {}
			
			bool IsLeaf() const {
				return (SplitDirection==-1);
			}
			
			IntCoordType Extension(int dim) const {
				return SCUnion.Cells2.Ext_[dim];
			}
			
			IntCoordType Max(int dim) const {
				return SCUnion.Cells2.Min_[dim]+SCUnion.Cells2.Ext_[dim];
			}
					
			IntCoordType Min(int dim) const {
				return SCUnion.Cells2.Min_[dim];
			}
			
			IntCoordType& Min(int i)  {
				return SCUnion.Cells2.Min_[i];
			}

			IntCoordType& Extension(int i)  {
				return SCUnion.Cells2.Ext_[i];
			}
			
			bool ContainsCell() const {
				return (IsLeaf() && (SCUnion.Cells2.Cell!=C::UndefinedCellRef()));
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
			


		};
		
		typedef typename std::vector<node > node_lst;
		
		class tmp_node {
			IntCoordType Min_[D];
			IntCoordType Ext_[D];
		public:
			
			IntLinkType ParentNode;
			int ChildNumber;
			IntLinkType StartIndex, EndIndex;
			
			tmp_node(	IntLinkType n,
						int c,
						IntLinkType s,
						IntLinkType e,
						IntCoordType* min,
						IntCoordType* ext) :
							ParentNode(n),ChildNumber(c), StartIndex(s), EndIndex(e) {
							
				for (int k=0;k<D;k++) {
					Min_[k]=min[k];
					Ext_[k]=ext[k];
				}
			};
							
			bool IsVoxel() const {
				for (int i=0;i<D;++i) if (Ext_[i]!=static_cast<IntCoordType>(1)) return false;
				return true;
			}
			
			IntCoordType Extension(int i) const {
				return Ext_[i];
			}

			IntCoordType Min(int i) const {
				return Min_[i];
			}

			IntCoordType& Extension(int i)  {
				return Ext_[i];
			}

			IntCoordType& Min(int i)  {
				return Min_[i];
			}
			
			int MaxExtension() const {
				int dir=0;
				for (int i=1;i<D;++i) if (Extension(i)>Extension(dir)) dir=i;
				return dir;
			}



		};
		
		
		IntCoordType Min_[D];
		IntCoordType Ext_[D];
		IntCoordType Max_[D];
		BoundaryType Boundaries[D];
		node_lst NodeList;
		
		void ClearPartition();



		template <class CoordsType> class CompareCoordinate {
			const int dir;
			const CoordsType& coords;
			
		public:
			CompareCoordinate(int d, const CoordsType& cs) :dir(d), coords(cs) {}
			bool operator()(IntCoordType c0, IntCoordType c1) const {
				return (coords[c0][dir]<coords[c1][dir]);
			}
		};
		
		template <class CoordsType> class IsSmallerCoordinate {
			const int dir;
			const IntCoordType coord;
			const CoordsType& coords;
		public:
			IsSmallerCoordinate(int d, IntCoordType co, const CoordsType& cs) :dir(d), coord(co), coords(cs) {}
			bool operator()(IntCoordType c0) const {
				return (coords[c0][dir]<coord);
			}
		};
		
		
	public:
		
		class subbox {
			
			friend class UpDownLinkTree;
			
			//std::vector<IntCoordType> min[D], max[D];
			//std::vector<IntLinkType> ParentNodes;
			
			const node* NodeRef;
			
		public:	
			
			IntCoordType Extension(int dim) const {
				return NodeRef->Extension(dim);
			}
			
			IntCoordType Max(int dim) const {
				return NodeRef->Max(dim);
			}
					
			IntCoordType Min(int dim) const {
				return NodeRef->Min(dim);
			}
			
			CellRefType Cell() const {
				return NodeRef->SCUnion.Cells2.Cell;
			}
			
			bool ContainsCell() const {
				return ((NodeRef->IsLeaf()) && (Cell()!=C::UndefinedCellRef()));
			}

		};

		/*IntLinkType UndefinedLink() const {
			return C::UndefinedLink();
		}
		
		CellRefType UndefinedCellRef() const {
			return C::UndefinedCellRef();
		}*/
		
		
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
			IntLinkType counter=0;
			for (typename C::IntLinkType i=0;i<NodeList.size();i++) {
				if (NodeList[i].IsLeaf()) counter++;
			}
			
			return counter;
		}
		
		IntLinkType NumberOfNodes() const {
			return NodeList.size();
		}
		
		IntLinkType TotalLeavesSurfaceArea() const {
			IntLinkType c=0;
			for (typename C::IntLinkType i=0;i<NodeList.size();i++) {
				if (NodeList[i].IsLeaf()) {
					for (int j=0;j<D;j++) c+=NodeList[i].AreaSize(j);
				}
			}
			return 2*c;
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
		
		IntLinkType NumberOfActiveVoxels() const {
			IntLinkType counter=0;
			for (typename C::IntLinkType i=0;i<NodeList.size();i++) {
				if (NodeList[i].ContainsCell()) counter++;
			}
			return counter;
		}
		
		double AverageNumberOfTraversedLeaves() const {
			return double(TotalLeavesSurfaceArea())/double(BoundingBoxSurfaceArea());
		}
		
		template <class V> subbox Access(V& pos , int Direction, bool DirectionSign) const;
		template <class V> bool GoToNeighborBox(subbox &, V& pos, int Direction, bool DirectionSign) const;
		
		void PrintStatistics(const std::string& FileName) const;

		double get_memory() const {
			return static_cast<double>(sizeof(node))*NodeList.size();
		}
				
	};

	
	template <class C> void UpDownLinkTree<C>::ClearPartition() {
		NodeList.clear();
		
	}
	
	
	
	template <class C> template <class CoordsType, class BncType> void UpDownLinkTree<C>::Setup(
			const CellRefType& start,
            const CellRefType& end,
            const CoordsType& CellCoords,
            const BncType& boundaries,
            const SplittingType& splitting_type,
            const SurfaceAreaHeuristicLambdaType& surface_area_heuristic_lambda) {

		
		///delete old partition
		ClearPartition();
		
		for (int i=0;i<D;i++) Boundaries[i]=boundaries[i];

		std::stack<tmp_node > TempNodeList;
		
		///setup list with all cells
		std::vector<IntLinkType> cell_lst;
		cell_lst.reserve(end-start);
		for (IntLinkType i=0;i<end-start;++i) cell_lst.push_back(i);
		
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
		
		//initialize with first subbox
		
		TempNodeList.push(tmp_node(	C::UndefinedLink(),0,0,cell_lst.size(),Min_,Ext_));
		
		
		
		//split subbox until finished
		while (!TempNodeList.empty()) {
			
			tmp_node& cl_it=TempNodeList.top();
			
			if (cl_it.ParentNode!=C::UndefinedLink()) {
				NodeList[cl_it.ParentNode].SCUnion.SC.Childs[cl_it.ChildNumber]=NodeList.size();
			}
			
			typename C::IntLinkType ParentNumber=NodeList.size();
			NodeList.push_back(node());
			node& cl=NodeList.back();
			
			
			cl.Parent=cl_it.ParentNode;
			
			//if node is leaf
			if ((cl_it.StartIndex==cl_it.EndIndex) || (cl_it.IsVoxel())) {
				if (cl_it.StartIndex!=cl_it.EndIndex) {
					cl.SCUnion.Cells2.Cell=start+cell_lst[cl_it.StartIndex];
				} else {
					cl.SCUnion.Cells2.Cell=C::UndefinedCellRef();
				}
				for (int k=0;k<D;++k) {
					cl.Min(k)=cl_it.Min(k);
					cl.Extension(k)=cl_it.Extension(k);
				}
				
				TempNodeList.pop();
				continue;
			}
			
			//calculate splitting plane

			int& SplitDirection=cl.SplitDirection;
			IntCoordType& SplitIndex=cl.SCUnion.SC.SplitCoordinate;
			
			switch(splitting_type) {
                case SPATIAL_MEDIAN: {
					SplitDirection=cl_it.MaxExtension();
					SplitIndex=cl_it.Extension(SplitDirection)/2+cl_it.Min(SplitDirection);
					break;
				}
				case OBJECT_MEDIAN: {
					SplitDirection=cl_it.MaxExtension();
					std::nth_element(&(cell_lst[cl_it.StartIndex]), &(cell_lst[(cl_it.EndIndex+cl_it.StartIndex)/2]), &(cell_lst[cl_it.EndIndex]), CompareCoordinate<CoordsType>(SplitDirection, CellCoords));
					SplitIndex=static_cast<typename C::IntCoordType>(CellCoords[cell_lst[(cl_it.EndIndex+cl_it.StartIndex)/2]][SplitDirection]);
					if (SplitIndex==cl_it.Min(SplitDirection)) SplitIndex++;
					break;
				}
				default: {
					double best_cost=std::numeric_limits<double>::max();
					
					IntLinkType N=cl_it.EndIndex-cl_it.StartIndex;

					for (int dir=0;dir<D;dir++) {
						IntCoordType ext=cl_it.Extension(dir);
						IntCoordType min=cl_it.Min(dir);
						
						std::vector<IntLinkType> v(ext, IntLinkType(0));
						
						double reciproc_sum=0.;
						for (int i=0;i<D;i++) if (i!=dir) reciproc_sum+=(double(1.)/double(cl_it.Extension(i)));
						
						for (IntLinkType l=cl_it.StartIndex;l<cl_it.EndIndex;l++) v[CellCoords[cell_lst[l]][dir]-min]++;
						
						IntLinkType NA=v[0];
						for (IntCoordType l=1;l<ext;l++) {
							IntLinkType NB=N-NA;
							double cost=double(N)/double(ext)+reciproc_sum*((double(l)/double(ext))*NA+(double(ext-l)/double(ext))*NB);
							if ((NA==0) || (NB==0)) cost*=surface_area_heuristic_lambda;
							if (cost<best_cost) {
								best_cost=cost;
								SplitDirection=dir;
								SplitIndex=min+l;
							}
							NA+=v[l];
						}
					}
					
				}
			}
				
			//partitioning
			IntLinkType SplitCellIndex=(
					std::partition(&(cell_lst[cl_it.StartIndex]), &(cell_lst[cl_it.EndIndex]), IsSmallerCoordinate<CoordsType>(SplitDirection, SplitIndex,CellCoords))
					-&(cell_lst[0]));
	       
			IntLinkType StartCellIndex=cl_it.StartIndex;
			
			IntCoordType Min0[D];
			for (int g=0;g<D;++g) Min0[g]=cl_it.Min(g);
			
			IntCoordType Ext0[D];
			for (int g=0;g<D;++g) Ext0[g]=cl_it.Extension(g);
			Ext0[SplitDirection]=(SplitIndex-Min0[SplitDirection]);
			
			cl_it.StartIndex=SplitCellIndex;
			cl_it.Min(SplitDirection)=SplitIndex;
			cl_it.Extension(SplitDirection)-=Ext0[SplitDirection];
			cl_it.ParentNode=NodeList.size()-1;
			cl_it.ChildNumber=1;
			
			TempNodeList.push(
				tmp_node(ParentNumber, 0,StartCellIndex,SplitCellIndex, Min0, Ext0)
			);
		}
	}
	
	
	template <class C> template <class V> typename UpDownLinkTree<C>::subbox UpDownLinkTree<C>::Access(V& pos, int Direction, bool DirectionSign) const {
		
		subbox sb;
		
		sb.NodeRef=&NodeList[0];
				
		while (!sb.NodeRef->IsLeaf()) {
			if (sb.NodeRef->SplitDirection==Direction) {
				if (DirectionSign) {
					sb.NodeRef=&NodeList[sb.NodeRef->SCUnion.SC.Childs[0]];
				} else {
					sb.NodeRef=&NodeList[sb.NodeRef->SCUnion.SC.Childs[1]];
				}
			} else {	
				if (pos[sb.NodeRef->SplitDirection]<sb.NodeRef->SCUnion.SC.SplitCoordinate) {
					sb.NodeRef=&NodeList[sb.NodeRef->SCUnion.SC.Childs[0]];
				} else {
					sb.NodeRef=&NodeList[sb.NodeRef->SCUnion.SC.Childs[1]];
				} 
			}	
		}
		
		int dim=Direction;
		for (int i=0;i<D-1;i++) { 
			dim++;
			if (dim==D) dim=0;
			pos[dim]+=(Min(dim)-sb.Min(dim));
		}
		pos[Direction]=(DirectionSign)?0:sb.Extension(Direction);
				
		return sb;
	}
	
	
	
	template <class C> template <class V> bool UpDownLinkTree<C>::GoToNeighborBox(subbox &sb, V& pos, int Direction, bool DirectionSign) const {

		if (DirectionSign) {
			if (sb.Max(Direction)==Max(Direction)) {
				if (Boundaries[Direction]==PERIODIC) {
					int dim=Direction;
					for (int i=0;i<D-1;i++) { 
						dim++;
						if (dim==D) dim=0;
						pos[dim]+=sb.Min(dim)-Min(dim);
					}
					sb=Access(pos, Direction, true);
					return false;
				} else {
					return true;
				}
			}
		} else {
			if (sb.Min(Direction)==Min(Direction)) {
				if (Boundaries[Direction]==PERIODIC) {
					int dim=Direction;
					for (int i=0;i<D-1;i++) { 
						dim++;
						if (dim==D) dim=0;
						pos[dim]+=sb.Min(dim)-Min(dim);
					}
					sb=Access(pos, Direction, false);
					return false;
				} else {
					return true;
				}
			}
		}
		
		//transform to absolute coords
		int dim=Direction;
		for (int i=0;i<D-1;i++) { 
			dim++;
			if (dim==D) dim=0;
			pos[dim]+=sb.Min(dim);
		}
		
		//Traverse up and go to other child
		
		while (true) {
			
			IntLinkType tmp=(sb.NodeRef-(&NodeList[0]));
			sb.NodeRef=&NodeList[sb.NodeRef->Parent];
			
			if (sb.NodeRef->SplitDirection==Direction) {
				if (sb.NodeRef->SCUnion.SC.Childs[DirectionSign]!=tmp) break;
			}
		}
		
		//go to other child
		sb.NodeRef=&NodeList[sb.NodeRef->SCUnion.SC.Childs[DirectionSign]];
		
		//traverse down
		while (!sb.NodeRef->IsLeaf()) {
			
			if (sb.NodeRef->SplitDirection==Direction) {
				sb.NodeRef=&NodeList[sb.NodeRef->SCUnion.SC.Childs[!DirectionSign]];
			} else {	
				if (pos[sb.NodeRef->SplitDirection]<sb.NodeRef->SCUnion.SC.SplitCoordinate) {
					sb.NodeRef=&NodeList[sb.NodeRef->SCUnion.SC.Childs[0]];
				} else {
					sb.NodeRef=&NodeList[sb.NodeRef->SCUnion.SC.Childs[1]];
				} 
			}	
					
		}
			
		//transform to relative coordinates
		pos[dim]-=sb.Min(dim);
		for (int i=0;i<D-2;i++) { 
			if (dim==0) dim=D;
			dim--;
			pos[dim]-=sb.Min(dim);
		}
		pos[Direction]=(DirectionSign)?0:sb.Extension(Direction);
		
		return false;
	}
	
	template <class C> void UpDownLinkTree<C>::PrintStatistics(const std::string& FileName) const {
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



#endif /*PARTITION4_UpDownLinkTreeH_*/
