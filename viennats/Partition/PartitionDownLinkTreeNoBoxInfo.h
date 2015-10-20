#ifndef PARTITION2_DownLinkTreeNoBoxInfoH_
#define PARTITION2_DownLinkTreeNoBoxInfoH_

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

	template <class C> class DownLinkTreeNoBoxInfo {
	
		typedef typename C::IntCoordType 			IntCoordType;
		typedef typename C::IntLinkType 			IntLinkType;
		typedef typename C::CellRefType 			CellRefType;
		static const int D=C::Dimension;
	
		class node {
		public:	
			union {
				IntLinkType Childs[2];
				CellRefType Cell;	
			};
			
			IntCoordType SplitCoordinate;
			int SplitDirection;
			
			node() : SplitDirection(-1) {}
			
			bool IsLeaf() const {
				return (SplitDirection==-1);
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
			
			friend class DownLinkTreeNoBoxInfo;
			
			std::vector<IntCoordType> min[D], max[D];
			std::vector<IntLinkType> ParentNodes;
			
			const node* NodeArray;
			
		public:	
			
			IntCoordType Extension(int dim) const {
				return max[dim].back()-min[dim].back();
			}
			
			IntCoordType Max(int dim) const {
				return max[dim].back();
			}
					
			IntCoordType Min(int dim) const {
				return min[dim].back();
			}
			
			CellRefType Cell() const {
				return NodeArray[ParentNodes.back()].Cell;
			}
			
			bool ContainsCell() const {
				return ((NodeArray[ParentNodes.back()].IsLeaf()) && (Cell()!=C::UndefinedCellRef()));
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
		
		template <class CellCoordsType> void Setup(const CellRefType&, const CellRefType&, const CellCoordsType&, const BoundaryType*);
		
		IntLinkType CountLeaves() const {
			IntLinkType counter=0;
			for (typename C::IntLinkType i=0;i<NodeList.size();i++) {
				if (NodeList[i].SplitDirection==-1) counter++;
			}
			
			return counter;
		}
		
		IntLinkType CountNodes() const {
							
			return NodeList.size();
		}
		
		IntLinkType UndefinedLink() const {
			return C::UndefinedLink();
		}
		
		CellRefType UndefinedCellRef() const {
			return C::UndefinedCellRef();
		}
		
		template <class V> subbox Access(V& pos , int Direction, bool DirectionSign) const;
		template <class V> bool GoToNeighborBox(subbox &, V& pos, int Direction, bool DirectionSign) const;
				
	};

	
	template <class C> void DownLinkTreeNoBoxInfo<C>::ClearPartition() {
		NodeList.clear();
		
	}
	
	
	
	

	template <class C> template <class CoordsType> void DownLinkTreeNoBoxInfo<C>::Setup(
			const CellRefType& start, 
			const CellRefType& end, 
			const CoordsType& CellCoords,
			const BoundaryType* boundaries) {
		
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
				NodeList[cl_it.ParentNode].Childs[cl_it.ChildNumber]=NodeList.size();
			}
			
			typename C::IntLinkType ParentNumber=NodeList.size();
			NodeList.push_back(node());
			node& cl=NodeList.back();
			
			//if node is leaf
			if ((cl_it.StartIndex==cl_it.EndIndex) || (cl_it.IsVoxel())) {
				if (cl_it.StartIndex!=cl_it.EndIndex) {
					cl.Cell=start+cell_lst[cl_it.StartIndex];
				} else {
					cl.Cell=C::UndefinedCellRef();
				}
				TempNodeList.pop();
				continue;
			}
			
			//calculate splitting plane

			int& SplitDirection=cl.SplitDirection;
			IntCoordType& SplitIndex=cl.SplitCoordinate;
			
			switch(C::PartitionMode) {
				case SpatialMedian: {
					SplitDirection=cl_it.MaxExtension();
					SplitIndex=cl_it.Extension(SplitDirection)/2+cl_it.Min(SplitDirection);
					break;
				}
				case ObjectMedian: {
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
							if ((NA==0) || (NB==0)) cost*=C::SurfaceAreaHeuristicLambda;
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
			cl_it.ParentNode=ParentNumber;
			cl_it.ChildNumber=1;
			
			TempNodeList.push(
				tmp_node(ParentNumber, 0,StartCellIndex,SplitCellIndex, Min0, Ext0)
			);

			
		}
	}
	
	
	template <class C> template <class V> typename DownLinkTreeNoBoxInfo<C>::subbox DownLinkTreeNoBoxInfo<C>::Access(V& pos, int Direction, bool DirectionSign) const {
		
		subbox sb;
		
		for (int i=0;i<D;i++) {
			sb.min[i].push_back(Min(i));
			sb.max[i].push_back(Max(i));
		}
		sb.NodeArray=&NodeList[0];
		sb.ParentNodes.push_back(0);
				
		while (!NodeList[sb.ParentNodes.back()].IsLeaf()) {
			const node& n=NodeList[sb.ParentNodes.back()];
			
			if (n.SplitDirection==Direction) {
				if (DirectionSign) {
					sb.ParentNodes.push_back(n.Childs[0]);
					sb.max[n.SplitDirection].push_back(n.SplitCoordinate);
				} else {
					sb.ParentNodes.push_back(n.Childs[1]);
					sb.min[n.SplitDirection].push_back(n.SplitCoordinate);
				}
			} else {	
				if (pos[n.SplitDirection]<n.SplitCoordinate) {
					sb.ParentNodes.push_back(n.Childs[0]);
					sb.max[n.SplitDirection].push_back(n.SplitCoordinate);
				} else {
					sb.ParentNodes.push_back(n.Childs[1]);
					sb.min[n.SplitDirection].push_back(n.SplitCoordinate);
				} 
			}	
		}
		
		int dim=Direction;
		for (int i=0;i<D-1;i++) { 
			dim++;
			if (dim==D) dim=0;
			pos[dim]+=Min(dim);
			pos[dim]-=sb.Min(dim);
		}
		pos[Direction]=(DirectionSign)?0:sb.Extension(Direction);
				
		return sb;
	}
	
	
	
	template <class C> template <class V> bool DownLinkTreeNoBoxInfo<C>::GoToNeighborBox(subbox &sb, V& pos, int Direction, bool DirectionSign) const {

		if (DirectionSign) {
			if (sb.max[Direction].back()==Max(Direction)) {
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
			if (sb.min[Direction].back()==Min(Direction)) {
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
			
			IntLinkType tmp=sb.ParentNodes.back();
			sb.ParentNodes.pop_back();
			const node& n=NodeList[sb.ParentNodes.back()];
			
			if (n.Childs[0]==tmp) {
				sb.max[n.SplitDirection].pop_back();
				if ((DirectionSign) && (n.SplitDirection==Direction)) {
					sb.ParentNodes.push_back(n.Childs[1]);
					sb.min[n.SplitDirection].push_back(n.SplitCoordinate);
					break;
				}
			} else {
				sb.min[n.SplitDirection].pop_back();
				if ((!DirectionSign) && (n.SplitDirection==Direction)) {
					sb.ParentNodes.push_back(n.Childs[0]);
					sb.max[n.SplitDirection].push_back(n.SplitCoordinate);
					break;
				}
			}
			
		}
		
		//traverse down
		while (!NodeList[sb.ParentNodes.back()].IsLeaf()) {
			const node& n=NodeList.at(sb.ParentNodes.back());
			
			if (n.SplitDirection==Direction) {
				if (DirectionSign) {
					sb.ParentNodes.push_back(n.Childs[0]);
					sb.max[n.SplitDirection].push_back(n.SplitCoordinate);
				} else {
					sb.ParentNodes.push_back(n.Childs[1]);
					sb.min[n.SplitDirection].push_back(n.SplitCoordinate);
				}
			} else {	
				if (pos[n.SplitDirection]<n.SplitCoordinate) {
					sb.ParentNodes.push_back(n.Childs[0]);
					sb.max[n.SplitDirection].push_back(n.SplitCoordinate);
				} else {
					sb.ParentNodes.push_back(n.Childs[1]);
					sb.min[n.SplitDirection].push_back(n.SplitCoordinate);
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
	
	
}



#endif /*PARTITION2_DownLinkTreeNoBoxInfoH_*/
