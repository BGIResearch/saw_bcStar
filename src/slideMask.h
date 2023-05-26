/* Copyright (C) BGI-Reasearch - All Rights Reserved
* Unauthorized copying of this file, via any medium is strictly prohibited
* Proprietary and confidential
* Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
*/
#ifndef SLIDE_MASK_H
#define SLIDE_MASK_H

#include <unordered_map>
#include <map>
#include <vector>
#include <string>
#include "common.h"
#include <iostream>
#include <fstream>
#include <sstream>
//#include <boost/serialization/string.hpp>
//#include <boost/serialization/unordered_map.hpp>
#include "util.h"

using namespace std;

typedef struct slideRange{
	uint32_t colStart;
	uint32_t colEnd;
	uint32_t rowStart;
	uint32_t rowEnd;
}slideRange;

typedef struct Position
{
//	friend class boost::serialization::access;
	uint8_t fov_c;
	uint8_t fov_r;
	uint32_t x;
	uint32_t y;
	
	template <typename Archive>
	void serialize(Archive &ar, const unsigned int version){
		ar & fov_c ;
		ar & fov_r ;
		ar & x ;
		ar & y ;
	}
}Position;

typedef struct Position1{
//	friend class boost::serialization::access;
#ifdef OPT_MEMORY_PLAN_B
  	uint64_t x:19;
  	uint64_t count:25;
  	uint64_t y:19;
  	uint64_t fullTag:1;
//	uint32_t x;
//	uint32_t y;
#else
	uint32_t x;
	uint32_t y;
#endif
	template <typename Archive>
	void serialize(Archive &ar, const unsigned int version){
		ar & x;
		ar & y;
	}
}Position1;
typedef struct Position2{
  // rx: 0x12345  ry: 0x67890
  // x: 0x12346789  y: 0x50
  uint32_t x;
  uint8_t y;
  
  template <typename Archive>
  void serialize(Archive &ar, const unsigned int version){
	ar & x;
	ar & y;
  }
  uint32_t getX(){
    return ((x>>16)<<4) | (y>>4);
  }
  uint32_t getY(){
    return ((x & 0xffff)<<4) | (y & 0xf);
  }
  void set(uint32_t rx,uint32_t ry){
	x=(((rx & 0xfffff)>>4)<<16) | ((ry & 0xfffff)>>4);
	y=((rx & 0xf)<<4) | (ry & 0xf);
  }
  Position2(uint32_t rx,uint32_t ry){
    x=(((rx & 0xfffff)>>4)<<16) | ((ry & 0xfffff)>>4);
    y=((rx & 0xf)<<4) | (ry & 0xf);
  }
  Position2(){
    x=0;
    y=0;
  }
}Position2;
//added by gc
//used for serializing
//class Position2{
//public:
//    Position2(){
//        x=-1;
//        y=-1;
//    }
//    Position2(Position1 position1){
//        x=position1.x;
//        y=position1.y;
//    }
//
//    Position2(int i,int i1){
//        x=i;
//        y=i1;
//    }
//
//    uint32_t x;
//    uint32_t y;
//
//    template<class B> void serialize(B &buf) const{
//        buf<<x<<y;
//    }
//
//    template<class B> void parse(B &buf){
//        buf>>x>>y;
//    }
//};

class SlideMask {
public:
	SlideMask(string maskFile, bool isSeq500, int turnFovDegree, int fovGap = FOV_GAP);
	~SlideMask();
	void makeIdx(string maskFile, int turnFovDegree);
	pair<int32_t, int32_t> getIdx(int idx, int fovCol,  int fovRow);
	bool isInIntermediate(int idx);
	slideRange getSlideRange(int minCol, int maxCol, int minRow, int maxRow);
	uint16 getSlidePitch(std::string platform);
private:
	pair<int32_t, int32_t> getShift(int block_c, int block_r);
	//pair<int, int> getPosition(pair<int, int>& shit, int dnb_c, int dnb_r, int turnFovDegree);
	//void addIdx(int block_c, int block_r, int turnFovDegree);
public:
//	unordered_map<long, pair<int, int>> dnbIdx;
	pair<int, int>* dnbIdx;
	pair<int, int> maskOrder[10][10];
	pair<int, int> dnbPos;
	map<int, int> rowOrder;
	map<int, int> colOrder;
	int maxdnbCol = 0;
	int maxdnbRow = 0;
	int totaldnbs = 0;
	int intOffset = 0;
	int idx = 0;
	long  idxCol;
	long  idxRow;
	int turnFovDegree;
	int fovGap;
	bool isSeq500;
private:
	map<std::string, uint16> slidePitchMap = {
		{"SEQ500", 900},
		{"T1", 715},
		{"T5", 600},
		{"T10", 500}
	};
};

#endif
