#ifndef FILE_IO_HPP_
#define FILE_IO_HPP_

#include <iostream>
#include <cmath>
#include <climits>
#include <cstdint>
#include <fstream>

#include "grid.hpp"
//#include "kernel.hpp"



	namespace lvlset {
	//class DefaultLevelSetTraitsType;
	class GridTraitsType;
	class LevelSetTraitsType;
	template <class GridTraitsType, class LevelSetTraitsType>
	class levelset;




	template <class GridTraitsType, class LevelSetTraitsType>
	void importLevelSetData(
		std::vector<unsigned int>* startIndices,
		std::vector<unsigned int>* runTypes,
		std::vector<int>* runBreaks,
		std::vector<double>& distances,
		std::string path){
	    std::ifstream fin(path);
	    if(!fin.is_open()) { std::cout << "ERROR: Couldn't open the file: " << path << "\n";return;}
	    std::cout << "Reading in LevelSet from " << path << std::endl;
	    char buff[11] = {};
	    char byte;

			//type & constant definitions
			typedef typename levelset<GridTraitsType, LevelSetTraitsType>::size_type size_type;
			typedef typename levelset<GridTraitsType, LevelSetTraitsType>::index_type index_type;
			typedef typename levelset<GridTraitsType, LevelSetTraitsType>::value_type value_type;
			const size_type POS_PT  = levelset<GridTraitsType, LevelSetTraitsType>::POS_PT;
			const size_type NEG_PT = levelset<GridTraitsType, LevelSetTraitsType>::NEG_PT;
			const size_type UNDEF_PT = levelset<GridTraitsType, LevelSetTraitsType>::UNDEF_PT;

	    uint32_t uInt;
	    fin.read(buff, 11);
			//std::cout << buff << std::endl << "LVSTx" << std::endl<< strncmp((char *)&buff, "LVSTx", 5) << std::endl;
	    if(strncmp((char *)&buff, "LVSTx", 5)) {std::cout << "Error: File is not a levelset file." << std::endl; return;}
	    const int dimension = buff[5]-48;
	    if(LVST_FILE_VERSION_NUMBER !=  buff[6]-48) std::cout << "WARNING: File version does not match!" << std::endl;
	    if(bigEndian() != buff[7]-48) std::cout << "WARNING: File was written in a different byte order than it is being read. Results may be incorrect!" << std::endl;
	    int bits_per_run = buff[8]-48;
	    int bits_per_distance = buff[9]-48;
	    std::cout << "Dimensions: " << dimension << std::endl << "Bits per runtype:" << bits_per_run << std::endl << "Bits per distance:" << bits_per_distance << std::endl;
	    std::cout << "Bytes per min/max: " << (buff[10] >> 4) << ", Bytes per delta: " << (buff[10] & 0xF) << std::endl;
	    fin.seekg(int(fin.tellg()) + ((buff[10]>>4)*2+1)*dimension + (buff[10] & 0xF));
			/*for(int i=dimension;i--;){
	      int32_t grid_min, grid_max;
	      boundary_type bcondition = 0;
	      fin.read((char *)&grid_min, buff[10] >> 4);
	      fin.read((char *)&grid_max, buff[10] >> 4);
	      fin.read((char *)&bcondition, 1);
	      std::cout << "Dimension: " << i << std::endl << "Grid min: " << grid_min << std::endl
	          << "Grid max: " << grid_max << std::endl
	          << "Boundary Condition: " << bcondition << std::endl;
	    }

	    double gridDelta; //TODO
	    fin.read((char *)&gridDelta, buff[10] & 0xF);
	    std::cout << "Delta: " << gridDelta << std::endl;*/

	    for(int i=dimension;i--;){
	      //std::vector<size_type>& startIndices = sub_levelsets[0].start_indices[i];
	      //std::vector<size_type>& runTypes = sub_levelsets[0].runtypes[i];
	      //std::vector<index_type>& runBreaks = sub_levelsets[0].runbreaks[i];
	      uint32_t num_st_indices, num_run_types, num_run_breaks;
	      int32_t bytesPerStIndex, bytesPerRnType, bytesPerRnBreak;
	      //reading in the HRLE header
	      fin.read(&byte, 1);
	      fin.read((char *)&num_st_indices, 4);
	      fin.read((char *)&num_run_types, 4);
	      fin.read((char *)&num_run_breaks, 4);

	      bytesPerStIndex = (byte & 0x3) +1;
	      bytesPerRnType = (byte >> 2 & 0x3) +1;
	      bytesPerRnBreak = (byte >> 4 & 0x3) +1;

	      std::cout << bytesPerStIndex << " byte(s) per start index." << std::endl
	          << bytesPerRnType << " byte(s) per runtype." << std::endl
	          << bytesPerRnBreak << " byte(s) per runbreak." << std::endl;

	      uint32_t values_read = 0;
	      //reading start indices
	      //if(startIndices.size() > 0) startIndices[i].clear();
	      for(unsigned int x = 0; x < num_st_indices; x++){
	        uInt = 0;
	        fin.read((char*) &uInt, bytesPerStIndex);
	        startIndices[i].push_back(uInt);
	        values_read++;
	      }
	      std::cout << "Dimension " << i << std::endl << "\t" << values_read << " of " << num_st_indices << " start indices read." << std::endl;

	      uint32_t count = 0;
	      values_read = 0;
	      //reading runtypes
	      //if(runTypes.size() > 0) runTypes[i].clear();
	      for(int y = 0; y < std::ceil(num_run_types/4.0); y++){
	        fin.read(&byte, 1);
	        for(int z = 4; z--;){
	          if(values_read == num_run_types) break;
	          uInt = byte >> z*bits_per_run & 0x3;
	          if(uInt == 1) runTypes[i].push_back(POS_PT);
	          else if(uInt == 3) runTypes[i].push_back(levelset<GridTraitsType, LevelSetTraitsType>::NEG_PT);
	          else if(uInt == 2) runTypes[i].push_back(levelset<GridTraitsType, LevelSetTraitsType>::UNDEF_PT);
	          else if(uInt == 0) runTypes[i].push_back(101), count++;
	          values_read++;
	        }
	      }
	      std::cout << "\t" << values_read << " of " << num_run_types << " runtypes read." << " Defined runtypes: " << count << std::endl;
	      //reading defined runtypes
	      for(unsigned int j=0; j<runTypes[i].size(); j++){
	        if(runTypes[i][j] == 101){
	          uInt = 0;
	          fin.read((char *) &uInt, bytesPerRnType);
	          runTypes[i][j] = uInt;
	        }
	      }

	      values_read = 0;

	      //Allocate the exact amount of bytes you need for a runbreak
	      //Example: We have 2 bytes and the runbreak is -7, which corresponds to FFF9 (Two's complement).
	      //Now if we read it into a 4 byte integer we will have 0000 FFF9, which is not interpreted as -7, but as a positive number.
	      //-7 as a 4 byte integer is: FFFF FFF9.

	      //int8_t *sInt = new int8_t[bytesPerRnBreak];
	      int8_t *sInt = (int8_t *) malloc(bytesPerRnBreak*sizeof(int8_t));
	      //reading runbreaks
	      //if(runBreaks.size() > 0) runBreaks[i].clear();
	      for(unsigned int z = 0; z < num_run_breaks; z++){
	        fin.read((char *) sInt, bytesPerRnBreak);
	        runBreaks[i].push_back(*sInt);
	        values_read++;
	      }
	      //delete[] sInt;
	      free(sInt);
	      std::cout << "\t" << values_read << " of " << num_run_breaks << " runbreaks." << std::endl;
	    }

	    uint32_t num_distances = 0, values_read = 0;
	    //reading the number of distances to read
	    fin.read((char *)&num_distances, 4);
	    //std::vector<value_type> & distances = sub_levelsets[0].distances;
	    //
	    double value = (std::pow(2, bits_per_distance)-1)/2;
	    int count = CHAR_BIT/bits_per_distance;
	    int num_reads = std::ceil(num_distances/count);
	    char mask = 0xFF >> (CHAR_BIT - bits_per_distance);
	    //reading distances
	    //if(distances.size() > 0) distances.clear();
			uint8_t uInt8;
	    for(int i = 0; i < num_reads; i++){
	      fin.read(&byte, 1);
	      for(unsigned int z = count; z--;){
	        if(values_read == num_distances) break; //if distances are odd, skip padding bits
	        uInt8 = byte >> z*bits_per_distance & mask;
	        distances.push_back(uInt8 / value - 1);
	        values_read++;
	      }
	    }
	    std::cout << values_read << " of " << num_distances << " distances read." << std::endl;

	    if(fin.fail()) std::cout << "ERROR: Couldn't read from file: " << path << std::endl;
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
