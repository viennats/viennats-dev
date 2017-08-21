#ifndef PARTITION_NeighborLinksArraysH_
#define PARTITION_NeighborLinksArraysH_

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
#include <iomanip>
#include "Partition.h"
//#include "../Time.h"

namespace partition {

	template <class C> class NeighborLinksArrays {

		typedef typename C::IntCoordType 		IntCoordType;
		typedef typename C::IntLinkType 		IntLinkType;
		typedef typename C::CellRefType 		CellRefType;
		static const int D=C::Dimension;



		class leaf {
			IntCoordType Min_[D];
			IntCoordType Ext_[D];

		public:

			IntCoordType Extension(int dim) const {
				return Ext_[dim];
			}

			IntCoordType Max(int dim) const {
				return Min_[dim]+Ext_[dim];
			}

			IntCoordType Min(int dim) const {
				return Min_[dim];
			}

			IntCoordType& Min(int i)  {
				return Min_[i];
			}

			IntCoordType& Extension(int i)  {
				return Ext_[i];
			}

			bool IsVoxel() const {
				for (int i=0;i<D;++i) if (Extension(i)!=static_cast<IntCoordType>(1)) return false;
				return true;
			}

			int MaxExtension() const {
				int dir=0;
				for (int i=1;i<D;++i) if (Extension(i)>Extension(dir)) dir=i;
				return dir;
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


			IntLinkType NeighbourStartIndex[D];

			CellRefType Cell;

			leaf(const IntCoordType* min_,const IntCoordType* ext_) : Cell(C::UndefinedCellRef()) {
				for (int k=0;k<D;k++) Min_[k]=min_[k];
				for (int k=0;k<D;k++) Ext_[k]=ext_[k];
			}
		};

		typedef typename std::vector<leaf> leaf_lst;

		class tmp_node {
		public:
			leaf Subbox;
			IntLinkType StartIndex, EndIndex;

			tmp_node(const leaf& cl,IntLinkType start,IntLinkType end):Subbox(cl),StartIndex(start), EndIndex(end) {};
		};


		//double Time;
		IntCoordType Min_[D];
		IntCoordType Ext_[D];
		IntCoordType Max_[D];
		BoundaryType Boundaries[D];
		leaf_lst LeafList;


		std::vector<IntLinkType> NeighbourConnections[2*D];
		std::vector<IntLinkType> AccessList[2*D];

		void AddToSubboxList(const leaf&,  std::vector<IntLinkType>*);

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

		/*const leaf& Leaf(IntLinkType i) const {
			return LeafList[i];
		}*/

	public:

		class subbox {

			friend class NeighborLinksArrays;

			const leaf* Subbox;
		public:

			subbox(const leaf& s):Subbox(&s) {}
			subbox() {}

			const subbox& operator=(const subbox& b) {
				Subbox=b.Subbox;
				return *this;
			}

			IntCoordType Extension(int dim) const {
				return Subbox->Extension(dim);
			}

			IntCoordType Max(int dim) const {
				return Subbox->Min(dim)+Subbox->Extension(dim);
			}

			IntCoordType Min(int dim) const {
				return Subbox->Min(dim);
			}

			CellRefType Cell() const {
				return Subbox->Cell;
			}

			bool ContainsCell() const {
				return (Cell()!=C::UndefinedCellRef());
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

		double get_memory() const {
			double mem=0;
			for (int i=0;i<2*D;++i) {
				mem+=static_cast<double>(sizeof(IntLinkType))*NeighbourConnections[i].size();
				mem+=static_cast<double>(sizeof(IntLinkType))*AccessList[i].size();
			}
			mem+=static_cast<double>(sizeof(leaf))*LeafList.size();
			return mem;
		}


		template <class CellCoordsType, class BncType> void Setup(const CellRefType&, const CellRefType&, const CellCoordsType&, const BncType &, const SplittingType& ,const SurfaceAreaHeuristicLambdaType& );

		IntLinkType NumberOfLeaves() const {return LeafList.size();}

		IntLinkType NumberOfNodes() const {return LeafList.size();}

		IntLinkType TotalLeavesSurfaceArea() const {
			IntLinkType c=0;
			for (typename leaf_lst::const_iterator it=LeafList.begin();it!=LeafList.end();++it) {
				for (int i=0;i<D;++i) {
					c+=it->AreaSize(i);
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
			IntLinkType c=0;
			for (typename leaf_lst::const_iterator it=LeafList.begin();it!=LeafList.end();++it) {
				if (it->Cell!=C::UndefinedCellRef()) c++;
			}
			return c;
		}

		double AverageNumberOfTraversedLeaves() const {
			return double(TotalLeavesSurfaceArea())/double(BoundingBoxSurfaceArea());
		}


		void WriteToFileDX(const std::string&) const;

		void PrintStatistics(const std::string&) const;

		template <class V> subbox Access(V& pos , int Direction, bool DirectionSign) const;
		template <class V> bool GoToNeighborBox(subbox &, V& pos, int Direction, bool DirectionSign) const;

		bool Test() const;


	};

	template <class C> void NeighborLinksArrays<C>::AddToSubboxList(const leaf& c, std::vector<IntLinkType>* tmp_access) {

		IntLinkType cl_num=LeafList.size();
		LeafList.push_back(c);
		leaf& cl=LeafList.back();

		for (int dim=0;dim<D;dim++) {
			cl.NeighbourStartIndex[dim]=NeighbourConnections[dim].size();
			NeighbourConnections[dim].resize(NeighbourConnections[dim].size()+cl.AreaSize(dim),C::UndefinedLink());
			NeighbourConnections[dim+D].resize(NeighbourConnections[dim+D].size()+cl.AreaSize(dim),C::UndefinedLink());

			for (IntLinkType pos_cl=cl.NeighbourStartIndex[dim];pos_cl<NeighbourConnections[dim].size();pos_cl++) {
				IntCoordType coords[D-1];
				int dimx=dim+1;
				{
					if (dimx==D) dimx=0;
					coords[D-2]=pos_cl-cl.NeighbourStartIndex[dim];

					for (int s=0;s<D-2;s++) {
						coords[s]=coords[D-2]%cl.Extension(dimx)+cl.Min(dimx);
						coords[D-2]/=cl.Extension(dimx);
						dimx++;
						if (dimx>=D) dimx=0;
					}
					coords[D-2]+=cl.Min(dimx);
				}

				//dimx=dim+D-1 (mod D)

				IntLinkType pos_border=(coords[D-2]-Min(dimx));
				for (int s=D-3;s>=0;s--) {
					if (dimx==0) dimx=D;
					dimx--;
					pos_border*=Extension(dimx);
					pos_border+=(coords[s]-Min(dimx));
				}

				//dimx=dim+D-1-(D-2)=dim+1 (mod D)

				if (cl.Min(dim)==Min(dim)) {
					tmp_access[dim+D][pos_border]=cl_num;
					tmp_access[dim][pos_border]=cl_num;
				} else {

					NeighbourConnections[dim+D][pos_cl]=tmp_access[dim][pos_border];
					tmp_access[dim][pos_border]=cl_num;

					leaf& cl2=LeafList[NeighbourConnections[dim+D][pos_cl]];

					dimx-=2;

					//dimx=dim-1 (mod D)

					if (dimx<0) dimx+=D;
					IntLinkType pos_cl2=(coords[D-2]-cl2.Min(dimx));
					for (int s=D-3;s>=0;s--) {
						if (dimx==0) dimx=D;
						dimx--;
						pos_cl2*=cl2.Extension(dimx);
						pos_cl2+=(coords[s]-cl2.Min(dimx));
					}

					pos_cl2+=cl2.NeighbourStartIndex[dim];
					NeighbourConnections[dim][pos_cl2]=cl_num;
				}

			}
		}
	}



	template <class C> void NeighborLinksArrays<C>::ClearPartition() {
		for (int i=0;i<2*D;i++) {
			NeighbourConnections[i].clear();
			AccessList[i].clear();
		}
		LeafList.clear();
	}


	template <class C> template <class CoordsType, class BncType> void NeighborLinksArrays<C>::Setup(
			const CellRefType& start,
			const CellRefType& end,
			const CoordsType& CellCoords,
			const BncType& boundaries,
			const SplittingType& splitting_type,
			const SurfaceAreaHeuristicLambdaType& surface_area_heuristic_lambda) {

		//Time=my::time::GetTime();

		///delete old partition
		ClearPartition();

		for (int i=0;i<D;i++) Boundaries[i]=boundaries[i];

		std::stack<tmp_node> TempSubboxList;

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

		//setup borders
		for (int g=0;g<D;g++) {
			AccessList[g].resize(AreaSize(g), C::UndefinedLink());
			AccessList[g+D].resize(AreaSize(g),C::UndefinedLink());
		}

		//initialize with first leaf
		TempSubboxList.push(tmp_node(leaf(&Min_[0],&Ext_[0]),0,cell_lst.size()));

		int SplitDirection(0);
		IntCoordType SplitIndex(0);


		//split leaf until finished
		while (!TempSubboxList.empty()) {

			tmp_node& cl_it=TempSubboxList.top();
			leaf& cl=cl_it.Subbox;

			if ((cl_it.StartIndex==cl_it.EndIndex) || (cl.IsVoxel())) {
				if (cl_it.StartIndex!=cl_it.EndIndex) cl.Cell=start+cell_lst[cl_it.StartIndex];
				AddToSubboxList(cl, AccessList);
				TempSubboxList.pop();
				continue;
			}


			//calculate splitting plane

			switch(splitting_type) {
				case SPATIAL_MEDIAN: {
					SplitDirection=cl.MaxExtension();
					SplitIndex=cl.Extension(SplitDirection)/2+cl.Min(SplitDirection);
					break;
				}
				case OBJECT_MEDIAN: {
					SplitDirection=cl.MaxExtension();//MaxIndex(cl.Max-cl.Min);
					std::nth_element(&(cell_lst[cl_it.StartIndex]), &(cell_lst[(cl_it.EndIndex+cl_it.StartIndex)/2]), &(cell_lst[cl_it.EndIndex]), CompareCoordinate<CoordsType>(SplitDirection, CellCoords));
					SplitIndex=static_cast<typename C::IntCoordType>(CellCoords[cell_lst[(cl_it.EndIndex+cl_it.StartIndex)/2]][SplitDirection]);
					if (SplitIndex==cl.Min(SplitDirection)) SplitIndex++;
					break;
				}
				default: {
					double best_cost=std::numeric_limits<double>::max();

					IntLinkType N=cl_it.EndIndex-cl_it.StartIndex;

					for (int dir=0;dir<D;dir++) {
						IntCoordType ext=cl.Extension(dir);
						IntCoordType min=cl.Min(dir);

						std::vector<IntLinkType> v(ext, IntLinkType(0));

						double reciproc_sum=0.;
						for (int i=0;i<D;i++) if (i!=dir) reciproc_sum+=(double(1.)/double(cl.Extension(i)));

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
			for (int g=0;g<D;++g) Min0[g]=cl.Min(g);

			IntCoordType Ext0[D];
			for (int g=0;g<D;++g) Ext0[g]=cl.Extension(g);
			Ext0[SplitDirection]=(SplitIndex-Min0[SplitDirection]);

			cl_it.StartIndex=SplitCellIndex;
			cl.Min(SplitDirection)=SplitIndex;
			cl.Extension(SplitDirection)-=Ext0[SplitDirection];

			TempSubboxList.push(tmp_node(leaf(Min0, Ext0),StartCellIndex,SplitCellIndex));

		}

		//if periodic boundaries, link opposite boundaries

		for (int dim=0;dim<D;dim++) {

			if (Boundaries[dim]==PERIODIC) {

				for (IntLinkType pos_border=0;pos_border < this->AreaSize(dim);++pos_border) {

					leaf& cl_a=LeafList[AccessList[dim][pos_border]];
					leaf& cl_b=LeafList[AccessList[dim+D][pos_border]];

					IntCoordType coords[D-1];
					int dimx=dim+1;
					{
						if (dimx==D) dimx=0;
						coords[D-2]=pos_border;

						for (int s=0;s<D-2;s++) {
							coords[s]=coords[D-2]%Extension(dimx)+Min(dimx);
							coords[D-2]/=Extension(dimx);
							dimx++;
							if (dimx==D) dimx=0;
						}
						coords[D-2]+=Min(dimx);
					}

					IntLinkType pos_a=(coords[D-2]-cl_a.Min(dimx));
					IntLinkType pos_b=(coords[D-2]-cl_b.Min(dimx));

					for (int s=D-3;s>=0;s--) {
						if (dimx==0) dimx=D;
						dimx--;
						pos_a*=cl_a.Extension(dimx);
						pos_a+=(coords[s]-cl_a.Min(dimx));
						pos_b*=cl_b.Extension(dimx);
						pos_b+=(coords[s]-cl_b.Min(dimx));
					}
					pos_a+=cl_a.NeighbourStartIndex[dim];
					pos_b+=cl_b.NeighbourStartIndex[dim];

					//assert(NeighbourConnections[dim][pos_a]==UndefinedLink());
					//assert(NeighbourConnections[dim+D][pos_b]==UndefinedLink());

					NeighbourConnections[dim][pos_a]=AccessList[dim+D][pos_border];
					NeighbourConnections[dim+D][pos_b]=AccessList[dim][pos_border];

				}
			}
		}

		//Time=my::time::GetTime()-Time;
	}


	template <class C> void NeighborLinksArrays<C>::WriteToFileDX(const std::string& FileName) const {
		const int num_edges=D*(1<<(D-1));


		std::ofstream f;
		f.open(FileName.c_str());

		f<< "object 1 class array type float rank 1 shape " << D << " items "<< (LeafList.size()*(1<<D)) <<" data follows" << std::endl;

		for (typename leaf_lst::const_iterator it_cl=LeafList.begin();it_cl!=LeafList.end();it_cl++) {
			const leaf* cl=&(*it_cl);

			for (int i=0;i<(1<<D);i++) {
				for (int j=0;j<D;j++) f << ((((i>>j) & 1)==0)?cl->Min(j):(cl->Extension(j)+cl->Min(j))) << " ";
				f<<std::endl;
			}
		}

		f << "object 2 class array type int rank 1 shape 2 items "<< (LeafList.size()*num_edges) <<" data follows" << std::endl;
		IntLinkType count=0;
		for (typename leaf_lst::const_iterator it_cl=LeafList.begin();it_cl!=LeafList.end();it_cl++) {
			for (int i=0;i<((1<<D)-1);i++) {
				for (int j=i+1;j<(1<<D);j++) {
					int tmp=i^j;
					for (int k=0;k<D;k++) {
						if (tmp==(1<<k)) f << (i+count) << " " << (j+count) << std::endl;
					}
				}
			}
			count+=(1<<D);
		}
		f << "attribute \"element type\" string \"lines\"" << std::endl;
		f << "attribute \"ref\" string \"positions\"" << std::endl;

		f << "object \"data\" class field" << std::endl;
		f << "component \"positions\" value 1" << std::endl;
		f << "component \"connections\" value 2" << std::endl;
		f << "end" << std::endl;

		f.close();
	}

	template <class C> template <class V> typename NeighborLinksArrays<C>::subbox NeighborLinksArrays<C>::Access(V& pos, int Direction, bool DirectionSign) const {
		int dim=Direction;
		if (dim==0) dim=D;
		dim--;
		IntLinkType idx=(pos[dim]<0)?IntLinkType(0):((pos[dim]<Extension(dim))?static_cast<IntLinkType>(pos[dim]):static_cast<IntLinkType>(Extension(dim)-IntCoordType(1)));
		for (int i=0;i<D-2;i++) {
			if (dim==0) dim=D;
			dim--;
			idx*=Extension(dim);
			idx+=(pos[dim]<0)?IntLinkType(0):((pos[dim]<Extension(dim))?static_cast<IntLinkType>(pos[dim]):static_cast<IntLinkType>(Extension(dim)-IntCoordType(1)));
		}

		subbox sb(LeafList[AccessList[(DirectionSign)?(Direction+D):Direction][idx]]);

		pos[dim]+=(Min(dim)-sb.Min(dim));
		for (int i=0;i<D-2;i++) {
			dim++;
			if (dim==D) dim=0;
			pos[dim]+=(Min(dim)-sb.Min(dim));
		}

		pos[Direction]=(DirectionSign)?0:sb.Extension(Direction);

		return sb;
	}

	template <class C> template <class V> bool NeighborLinksArrays<C>::GoToNeighborBox(subbox &sb, V& pos, int Direction, bool DirectionSign) const {

		const IntLinkType *links=&(NeighbourConnections[(DirectionSign)?Direction:(Direction+D)][sb.Subbox->NeighbourStartIndex[Direction]]);

		if (links[0]==C::UndefinedLink()) return true;

		int dim=Direction;
		if (dim==0) dim=D;
		dim--;
		IntLinkType idx=(pos[dim]<0)?IntLinkType(0):((pos[dim]<sb.Extension(dim))?static_cast<IntLinkType>(pos[dim]):static_cast<IntLinkType>(sb.Extension(dim)-IntCoordType(1)));
		pos[dim]+=sb.Min(dim);
		for (int i=0;i<D-2;i++) {
			if (dim==0) dim=D;
			dim--;
			idx*=sb.Extension(dim);
			idx+=(pos[dim]<0)?IntLinkType(0):((pos[dim]<sb.Extension(dim))?static_cast<IntLinkType>(pos[dim]):static_cast<IntLinkType>(sb.Extension(dim)-IntCoordType(1)));
			pos[dim]+=sb.Min(dim);
		}

		sb.Subbox=&LeafList[links[idx]];

		pos[dim]-=sb.Min(dim);
		for (int i=0;i<D-2;i++) {
			dim++;
			if (dim==D) dim=0;
			pos[dim]-=sb.Min(dim);
		}

		pos[Direction]=(DirectionSign)?0:sb.Extension(Direction);

		return false;

	}


	template <class C>
	bool NeighborLinksArrays<C>::Test() const {
		for (IntLinkType it_cl=0;it_cl!=LeafList.size();it_cl++) {

			const leaf& cl=LeafList[it_cl];

			for (int dir=0;dir<D;dir++) {

				for (IntLinkType itArea=0;itArea<cl.AreaSize(dir);itArea++) {

					IntLinkType tmp=itArea;

					IntCoordType v[D];

					for (int l=(dir+1)%D; l!=dir; l=(l+1)%D) {
						v[l]=tmp%cl.Extension(l);
						tmp/=cl.Extension(l);
					}

					IntLinkType l0=NeighbourConnections[dir  ][itArea+cl.NeighbourStartIndex[dir]];
					IntLinkType l1=NeighbourConnections[dir+D][itArea+cl.NeighbourStartIndex[dir]];

					if (l0!=C::UndefinedLink()) {
						const leaf& cl0=LeafList[l0];

						if (!((cl.Max(dir)==cl0.Min(dir)) || ((Boundaries[dir]==PERIODIC) && (cl.Max(dir)==Max(dir)) && (cl0.Min(dir)==Min(dir))))) return false;
						for (int l=(dir+1)%D; l!=dir; l=(l+1)%D) if(!((cl.Max(l)>cl0.Min(l)) && (cl0.Max(l)>cl.Min(l)))) return false;

						IntLinkType x0=0;
						for (int l=(dir+D-1)%D; l!=dir; l=(l+D-1)%D) {
							x0*=cl0.Extension(l);
							x0+=(v[l]+cl.Min(l)-cl0.Min(l));
						}
						if (NeighbourConnections[dir+D][x0+cl0.NeighbourStartIndex[dir]]!=it_cl) return false;
					} else {
						if (Boundaries[dir]==PERIODIC) return false;
						if(cl.Max(dir)!=Max(dir)) return false;
						IntLinkType x0=0;
						for (int l=(dir+D-1)%D; l!=dir; l=(l+D-1)%D) {
							x0*=Extension(l);
							x0+=(v[l]+cl.Min(l)-Min(l));
						}
						if(AccessList[dir][x0]!=it_cl) return false;
					}

					if (l1!=C::UndefinedLink()) {
						const leaf& cl1=LeafList[l1];
						if (!((cl.Min(dir)==cl1.Max(dir)) || ((Boundaries[dir]==PERIODIC) && (cl.Min(dir)==Min(dir)) && (cl1.Max(dir)==Max(dir))))) return false;
						for (int l=(dir+1)%D; l!=dir; l=(l+1)%D) if(!((cl.Max(l)>cl1.Min(l)) && (cl1.Max(l)>cl.Min(l)))) return false;

						IntLinkType x1=0;
						for (int l=(dir+D-1)%D; l!=dir; l=(l+D-1)%D) {
							x1*=cl1.Extension(l);
							x1+=v[l]+cl.Min(l)-cl1.Min(l);
						}
						if (NeighbourConnections[dir ][x1+cl1.NeighbourStartIndex[dir]]!=it_cl) return false;
					} else {
						if (Boundaries[dir]==PERIODIC) return false;
						if(cl.Min(dir)!=Min(dir)) return false;
						IntLinkType x1=0;
						for (int l=(dir+D-1)%D; l!=dir; l=(l+D-1)%D) {
							x1*=Extension(l);
							x1+=(v[l]+cl.Min(l)-Min(l));
						}
						if(AccessList[dir+D][x1]!=it_cl) return false;
					}
				}
			}
		}

		return true;
	}


	template <class C> void NeighborLinksArrays<C>::PrintStatistics(const std::string& FileName) const {
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






#endif /*PARTITION_NeighborLinksArraysH_*/
