#ifndef FILE_IO_HPP_
#define FILE_IO_HPP_

#include <iostream>
#include <cmath>
#include <climits>
#include <cstdint>
#include <fstream>

#include "grid.hpp"
//#include "kernel.hpp"


#define MAJOR_FILE_LVLST_VERSION 0
#define MINOR_FILE_LVLST_VERSION 1
#define BITS_PER_RUN 2
#define BITS_PER_DISTANCE 4
#define BYTES_PER_INDEX 4


namespace lvlset {
//class DefaultLevelSetTraitsType;
//class GridTraitsType;
//class LevelSetTraitsType;
//template <class GridTraitsType, class LevelSetTraitsType>
//class levelset;

	union {
		uint16_t shortVar;    // binary  number of length 16 Bits
		uint8_t  charVar[2];  // 2 binary numbers, each 8 Bits
	} test_endianness;

	bool bigEndian(){
		test_endianness.shortVar = 0x8000; // MSB of 16
		return test_endianness.charVar[0] != 0;
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