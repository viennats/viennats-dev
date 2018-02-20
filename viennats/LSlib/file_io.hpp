#ifndef FILE_IO_HPP_
#define FILE_IO_HPP_

#include <iostream>
#include <cmath>
#include <climits>
#include <cstdint>
#include <fstream>

#include "grid.hpp"
//#include "kernel.hpp"

#define FILE_NAME "./exported_LVST.lvl"
#define LVST_FILE_VERSION_NUMBER 1
#define BITS_PER_RUNTYPE 2
#define BITS_PER_DISTANCE 4
//#define BITS_OVERFOW(x) CHAR_BIT%x!=0

#if BITS_PER_DISTANCE > CHAR_BIT
_Pragma ("GCC error  \"You are trying to use more BITS_PER_DISTANCE than one char can hold. This will result in an invalid levelset.\"")
#endif // BITS LONGER_THAN_CHAR


#define GRID_IN_FILE_HEADER
#define FORCE_BYTESIZE
#define BYTES_START_INDEX 4
#define BYTES_RUNTYPE 4
#define BYTES_RUNBREAK 4

namespace lvlset {
//class DefaultLevelSetTraitsType;
class GridTraitsType;
class LevelSetTraitsType;
template <class GridTraitsType, class LevelSetTraitsType>
class levelset;

	union {
		uint16_t shortVar;    // binary  number of length 16 Bits
		uint8_t  charVar[2];  // 2 binary numbers, each 8 Bits
	} test_endianness;

	bool bigEndian(){
		test_endianness.shortVar = 0x8000; // MSB of 16
		return test_endianness.charVar[0] != 0;
	}



	template <typename GTT, typename LSTT> void importLevelSet(const levelset<GTT, LSTT>& l, std::string path){

		//GridTraitsType<D> GridProperties(grid_min, grid_max, bnc, p.GridDelta);
		//lvlset::grid_type<GridTraitsType<D> > grid(GridProperties);

	}
/*
	template <class GTT, class D> template<typename A, typename B, typename C> GTT<D> fillLVSTData(const A& stIn, const A& rnTy, const B& rnBr, const C& dist, std::string path){
		std::ifstream fin("./exportedLVLSet_g.lvl");
		char buff[10] = {};
		fin.read(buff, 10);
		const int dim = buff[5]-48;
		typedef typename GridTraitsType<D>::index_type index_type;
		index_type * g_min = new index_type[dim];
		index_type * g_max = new index_type[dim];
		lvlset::boundary_type* bc = new lvlset::boundary_type [dim];
		std::cout << "D: " << D << ", Dim: " << dim << std::endl;
		for(int i=dim;i--;){
			fin.read((char *)&g_min[i], 4);
			fin.read((char *)&g_max[i], 4);
			fin.read((char *)&bc[i], 1);
		}
		double gridDelta; //TODO
		fin.read((char *)&gridDelta, sizeof(double));
	}
*/
	template<typename A, typename B> void fillGridProperties(A& grid_min, A& grid_max, B& boundary_conditions, double& gD, std::string path){

		std::ifstream fin(path);
		char buff[10] = {};
		char byte;

		uint32_t uInt;
		//int32_t sInt;
		fin.read(buff, 10);
		const int dimension = buff[5]-48;

		int bits_per_run = buff[8]-48;
		int bits_per_distance = buff[9]-48;
		std::cout << "Dimensions: " << dimension << std::endl << "Bits per runtype:" << bits_per_run << std::endl << "Bits per distance:" << bits_per_distance << std::endl;
		double gridDelta; //TODO
		fin.read((char *)&gridDelta, sizeof(double));
		std::cout << "Delta: " << gridDelta << std::endl;
		gD = gridDelta;

		for(int i=dimension;i--;){
			long pos;
			uint32_t num_st_indices, num_run_types, num_run_breaks;
			int32_t bytesPerStIndex, bytesPerRnType, bytesPerRnBreak;
			int32_t gd_min, gd_max;
			boundary_type bcondition = 0;
			//reading in the HRLE header
		 	fin.read(&byte, 1);
			fin.read((char *)&num_st_indices, 4);
			fin.read((char *)&num_run_types, 4);
			fin.read((char *)&num_run_breaks, 4);
			fin.read((char *)&gd_min, 4);
			fin.read((char *)&gd_max, 4);
			fin.read((char *)&bcondition, 1);
			grid_min[i] = gd_min;
			grid_max[i] = gd_max;
			boundary_conditions[i] = bcondition;

			bytesPerStIndex = (byte & 0x3) +1;
			bytesPerRnType = (byte >> 2 & 0x3) +1;
			bytesPerRnBreak = (byte >> 4 & 0x3) +1;

			std::cout << bytesPerStIndex << " byte(s) per start index." << std::endl
				  << bytesPerRnType << " byte(s) per runtype." << std::endl
				  << bytesPerRnBreak << " byte(s) per runbreak." << std::endl
				  << "Grid min: " << grid_min[i] << std::endl
				  << "Grid max: " << grid_max[i] << std::endl
				  << "Boundary Condition: " << bcondition << std::endl;

			pos = fin.tellg();
			fin.seekg(pos+num_st_indices*bytesPerStIndex);

			uint32_t count = 0, values_read = 0;
			//reading runtypes
			for(int y = 0; y < std::ceil(num_run_types/4.0); y++){
				fin.read(&byte, 1);
				for(int z = 4; z--;){
					if(values_read == num_run_types) break;
					uInt = byte >> z*bits_per_run & 0x3;
					if(uInt == 0) count++;
					values_read++;
				}
			}

			pos = fin.tellg();
			fin.seekg(pos+count*bytesPerRnType+bytesPerRnBreak*num_run_breaks);
			std::cout << fin.tellg();
		}
		fin.close();
	}

	char bitsToChar(std::bitset<CHAR_BIT> &s){
	    char c = 0;
	    for(unsigned int i = 0; i < CHAR_BIT; i++){
		if(s[i] == 1){
		    //unsigned char m = 0x1 << i; // is for string and bitset reversed
		    c |= 0x1 << i;
		}
	    }
	    return c;
	}

	char bitsToChar(std::string s){
	    char c = 0;
	    for(unsigned int i = 0; i < CHAR_BIT; i++){
		if(s[i] == '1'){
		    //unsigned char m = (int)pow(2, CHAR_BIT-1) >> i;// is for string and bitset reversed
		    c |= (int)std::pow(2, CHAR_BIT-1) >> i;
		}
	    }
	    return c;
	}



}
#endif /*FILE_IO_HPP_*/
